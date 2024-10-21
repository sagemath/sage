# sage_setup: distribution = sagemath-modules
# sage.doctest: needs sage.groups
"""
Representations of a semigroup

AUTHORS:

- Travis Scrimshaw (2015-11-21): initial version
- Siddharth Singh  (2020-03-21): signed representation
- Travis Scrimshaw (2024-02-17): tensor products
"""

##############################################################################
#       Copyright (C) 2015-2024 Travis Scrimshaw <tcscrims at gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  https://www.gnu.org/licenses/
##############################################################################

from sage.misc.abstract_method import abstract_method
from sage.misc.cachefunc import cached_method
from sage.structure.element import Element
from sage.combinat.free_module import CombinatorialFreeModule, CombinatorialFreeModule_Tensor
from sage.modules.with_basis.subquotient import SubmoduleWithBasis
from sage.categories.modules import Modules
from sage.matrix.constructor import matrix
from sage.modules.free_module_element import vector
from sage.modules.with_basis.subquotient import SubmoduleWithBasis, QuotientModuleWithBasis

class Representation_abstract:
    """
    Abstract base class for representations of semigroups.

    INPUT:

    - ``semigroup`` -- a semigroup
    - ``side`` -- (default: ``'left'``) whether this is a
      ``'left'`` or ``'right'`` representation
    - ``algebra`` -- (default: ``semigroup.algebra(self.base_ring())``)
      the semigroup algebra

    .. NOTE::

        This class should come before :class:`CombinatorialFreeModule` in the
        MRO in order for tensor products to use the correct class.
    """
    def __init__(self, semigroup, side, algebra=None):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: G = FreeGroup(3)
            sage: T = G.trivial_representation()
            sage: TestSuite(T).run()

        Verify that the dynamic mixin classes work::

            sage: SGA = SymmetricGroupAlgebra(QQ, 3)
            sage: S21 = SGA.specht_module([2,1])
            sage: T = tensor([S21, S21])
            sage: s21 = S21.frobenius_image()
            sage: T.frobenius_image() == s21.kronecker_product(s21)
            True
        """
        self._semigroup = semigroup
        if algebra is None:
            algebra = semigroup.algebra(self.base_ring())
        self._semigroup_algebra = algebra
        self._side = side
        if side not in ["left", "right", "twosided"]:
            raise ValueError("the side must be either 'left', 'right', or 'twosided'")
        self._left_repr = bool(side == "left" or side == "twosided")
        self._right_repr = bool(side == "right" or side == "twosided")

        if hasattr(self._semigroup_algebra, "_representation_mixin_class"):
            mixin = self._semigroup_algebra._representation_mixin_class
            # No need to do anything if it is already in the MRO
            if mixin not in self.__class__.__mro__:
                from sage.structure.dynamic_class import dynamic_class
                cat = self.category()
                # perhaps the category has not been initialized yet
                if not isinstance(self, cat.parent_class):
                    self.__class__ = dynamic_class(f"{type(self).__name__}_with_mixin",
                                                   (type(self), mixin),
                                                   doccls=type(self))
                else:
                    base = self.__class__.__base__  # strip off the category dynamic class
                    # recreate the dynamic class with adding the mixin
                    self.__class__ = dynamic_class(f"{base.__name__}_with_category",
                                                   (base, mixin, cat.parent_class),
                                                   doccls=base)

    def semigroup(self):
        """
        Return the semigroup whose representation ``self`` is.

        EXAMPLES::

            sage: G = SymmetricGroup(4)
            sage: M = CombinatorialFreeModule(QQ, ['v'])
            sage: on_basis = lambda g,m: M.term(m, g.sign())
            sage: R = G.representation(M, on_basis)
            sage: R.semigroup()
            Symmetric group of order 4! as a permutation group
        """
        return self._semigroup

    def semigroup_algebra(self):
        """
        Return the semigroup algebra whose representation ``self`` is.

        EXAMPLES::

            sage: G = SymmetricGroup(4)
            sage: M = CombinatorialFreeModule(QQ, ['v'])
            sage: on_basis = lambda g,m: M.term(m, g.sign())
            sage: R = G.representation(M, on_basis)
            sage: R.semigroup_algebra()
            Symmetric group algebra of order 4 over Rational Field
        """
        return self._semigroup_algebra

    def side(self):
        """
        Return whether ``self`` is a left, right, or two-sided representation.

        OUTPUT: the string ``'left'``, ``'right'``, or ``'twosided'``

        EXAMPLES::

            sage: G = groups.permutation.Dihedral(4)
            sage: R = G.regular_representation()
            sage: R.side()
            'left'
            sage: S = G.regular_representation(side='right')
            sage: S.side()
            'right'
            sage: R = G.sign_representation()
            sage: R.side()
            'twosided'
            sage: R = G.trivial_representation()
            sage: R.side()
            'twosided'
        """
        return self._side

    def invariant_module(self, S=None, **kwargs):
        r"""
        Return the submodule of ``self`` invariant under the action of ``S``.

        For a semigroup `S` acting on a module `M`, the invariant
        submodule is given by

        .. MATH::

            M^S = \{m \in M : s \cdot m = m \forall s \in S\}.

        INPUT:

        - ``S`` -- a finitely-generated semigroup (default: the semigroup
          this is a representation of)
        - ``action`` -- a function (default: :obj:`operator.mul`)
        - ``side`` -- ``'left'`` or ``'right'`` (default: :meth:`side()`);
          which side of ``self`` the elements of ``S`` acts

        .. NOTE::

            Two sided actions are considered as left actions for the
            invariant module.

        OUTPUT: :class:`~sage.modules.with_basis.invariant.FiniteDimensionalInvariantModule`

        EXAMPLES::

            sage: S3 = SymmetricGroup(3)
            sage: M = S3.regular_representation()
            sage: I = M.invariant_module()
            sage: [I.lift(b) for b in I.basis()]
            [() + (2,3) + (1,2) + (1,2,3) + (1,3,2) + (1,3)]

        We build the `D_4`-invariant representation inside of the regular
        representation of `S_4`::

            sage: D4 = groups.permutation.Dihedral(4)
            sage: S4 = SymmetricGroup(4)
            sage: R = S4.regular_representation()
            sage: I = R.invariant_module(D4)
            sage: [I.lift(b) for b in I.basis()]
            [() + (2,4) + (1,2)(3,4) + (1,2,3,4) + (1,3) + (1,3)(2,4) + (1,4,3,2) + (1,4)(2,3),
             (3,4) + (2,3,4) + (1,2) + (1,2,4) + (1,3,2) + (1,3,2,4) + (1,4,3) + (1,4,2,3),
             (2,3) + (2,4,3) + (1,2,3) + (1,2,4,3) + (1,3,4,2) + (1,3,4) + (1,4,2) + (1,4)]
        """
        if S is None:
            S = self.semigroup()
        side = kwargs.pop('side', self.side())
        if side == "twosided":
            side = "left"

        return super().invariant_module(S, side=side, **kwargs)

    def twisted_invariant_module(self, chi, G=None, **kwargs):
        r"""
        Create the isotypic component of the action of ``G`` on
        ``self`` with irreducible character given by ``chi``.

        .. SEEALSO::

            - :class:`~sage.modules.with_basis.invariant.FiniteDimensionalTwistedInvariantModule`

        INPUT:

        - ``chi`` -- list/tuple of character values or an instance
          of :class:`~sage.groups.class_function.ClassFunction_gap`
        - ``G`` -- a finitely-generated semigroup (default: the semigroup
          this is a representation of)

        This also accepts the first argument to be the group.

        OUTPUT: :class:`~sage.modules.with_basis.invariant.FiniteDimensionalTwistedInvariantModule`

        EXAMPLES::

            sage: G = SymmetricGroup(3)
            sage: R = G.regular_representation(QQ)
            sage: T = R.twisted_invariant_module([2,0,-1])
            sage: T.basis()
            Finite family {0: B[0], 1: B[1], 2: B[2], 3: B[3]}
            sage: [T.lift(b) for b in T.basis()]
            [() - (1,2,3), -(1,2,3) + (1,3,2), (2,3) - (1,2), -(1,2) + (1,3)]

        We check the different inputs work::

            sage: R.twisted_invariant_module([2,0,-1], G) is T
            True
            sage: R.twisted_invariant_module(G, [2,0,-1]) is T
            True
        """
        from sage.categories.groups import Groups
        if G is None:
            G = self.semigroup()
        elif chi in Groups():
            G, chi = chi, G
        side = kwargs.pop('side', self.side())
        if side == "twosided":
            side = "left"

        return super().twisted_invariant_module(G, chi, side=side, **kwargs)

    def representation_matrix(self, g, side=None, sparse=False):
        r"""
        Return the matrix representation of ``g`` acting on ``self``.

        EXAMPLES::

            sage: S3 = SymmetricGroup(3)
            sage: g = S3.an_element(); g
            (2,3)
            sage: L = S3.regular_representation(side='left')
            sage: R = S3.regular_representation(side='right')
            sage: R.representation_matrix(g)
            [0 0 0 1 0 0]
            [0 0 0 0 0 1]
            [0 0 0 0 1 0]
            [1 0 0 0 0 0]
            [0 0 1 0 0 0]
            [0 1 0 0 0 0]
            sage: L.representation_matrix(g)
            [0 0 0 1 0 0]
            [0 0 0 0 1 0]
            [0 0 0 0 0 1]
            [1 0 0 0 0 0]
            [0 1 0 0 0 0]
            [0 0 1 0 0 0]
            sage: A = S3.algebra(ZZ)
            sage: R.representation_matrix(sum(A.basis()), side='right')
            [1 1 1 1 1 1]
            [1 1 1 1 1 1]
            [1 1 1 1 1 1]
            [1 1 1 1 1 1]
            [1 1 1 1 1 1]
            [1 1 1 1 1 1]

        We verify tensor products agree::

            sage: T = tensor([L, R])
            sage: for g in S3:
            ....:     gL = L.representation_matrix(g, side='left')
            ....:     gR = R.representation_matrix(g, side='left')
            ....:     gT = T.representation_matrix(g, side='left')
            ....:     assert gL.tensor_product(gR) == gT

        Some examples with Specht modules::

            sage: SM = Partition([3,1,1]).specht_module(QQ)
            sage: SM.representation_matrix(Permutation([2,1,3,5,4]))
            [-1  0  1  0  1  0]
            [ 0  0  0 -1 -1 -1]
            [ 0  0  0  0  1  0]
            [ 0 -1 -1  0  0  1]
            [ 0  0  1  0  0  0]
            [ 0  0  0  0  0 -1]

            sage: SGA = SymmetricGroupAlgebra(QQ, 5)
            sage: SM = SGA.specht_module([(0,0), (0,1), (0,2), (1,0), (2,0)])
            sage: SM.representation_matrix(Permutation([2,1,3,5,4]))
            [-1  0  0  0  0  0]
            [ 0  0  1  0  0  0]
            [ 0  1  0  0  0  0]
            [ 1  0 -1  0 -1  0]
            [-1 -1  0 -1  0  0]
            [ 0  1  1  0  0 -1]
            sage: SM.representation_matrix(SGA([3,1,5,2,4]))
            [ 0  0  0  0  1  0]
            [-1  0  0  0  0  0]
            [ 0  0  0 -1  0  0]
            [ 1  0 -1  0 -1  0]
            [ 0  0  0  1  1  1]
            [-1 -1  0 -1  0  0]

            sage: SGA = SymmetricGroupAlgebra(QQ, 4)
            sage: SM = SGA.specht_module([3, 1])
            sage: all(SM.representation_matrix(g) * b.to_vector() == (g * b).to_vector()
            ....:     for b in SM.basis() for g in SGA.group())
            True

            sage: SGA = SymmetricGroupAlgebra(QQ, SymmetricGroup(4))
            sage: SM = SGA.specht_module([3, 1])
            sage: all(SM.representation_matrix(g) * b.to_vector() == (g * b).to_vector()
            ....:     for b in SM.basis() for g in SGA.group())
            True
        """
        if self.dimension() == float('inf'):
            raise NotImplementedError("only implemented for finite dimensional modules")

        B = self.basis()
        order = self.get_order()
        inv_order = {b: i for i, b in enumerate(order)}
        ret = matrix.zero(self.base_ring(), len(order), sparse=sparse)
        if side is None:
            if self._side == "twosided":
                side = "left"
            else:
                side = self._side
        use_left = side == "left"
        for i, k in enumerate(order):
            if use_left:
                temp = g * B[k]
            else:
                temp = B[k] * g
            for m, c in temp.monomial_coefficients(copy=False).items():
                if not use_left:
                    ret[i, inv_order[m]] = c
                else:
                    ret[inv_order[m], i] = c
        return ret

    @cached_method
    def character(self):
        r"""
        Return the character of ``self``.

        EXAMPLES::

            sage: SGA = SymmetricGroupAlgebra(QQ, 5)
            sage: SM = SGA.specht_module([3,2])
            sage: SM.character()
            (5, 1, 1, -1, 1, -1, 0)
            sage: matrix(SGA.specht_module(la).character() for la in Partitions(5))
            [ 1  1  1  1  1  1  1]
            [ 4  2  0  1 -1  0 -1]
            [ 5  1  1 -1  1 -1  0]
            [ 6  0 -2  0  0  0  1]
            [ 5 -1  1 -1 -1  1  0]
            [ 4 -2  0  1  1  0 -1]
            [ 1 -1  1  1 -1 -1  1]

            sage: SGA = SymmetricGroupAlgebra(QQ, SymmetricGroup(5))
            sage: SM = SGA.specht_module([3,2])
            sage: SM.character()
            Character of Symmetric group of order 5! as a permutation group
            sage: SM.character().values()
            [5, 1, 1, -1, 1, -1, 0]
            sage: matrix(SGA.specht_module(la).character().values() for la in reversed(Partitions(5)))
            [ 1 -1  1  1 -1 -1  1]
            [ 4 -2  0  1  1  0 -1]
            [ 5 -1  1 -1 -1  1  0]
            [ 6  0 -2  0  0  0  1]
            [ 5  1  1 -1  1 -1  0]
            [ 4  2  0  1 -1  0 -1]
            [ 1  1  1  1  1  1  1]
            sage: SGA.group().character_table()
            [ 1 -1  1  1 -1 -1  1]
            [ 4 -2  0  1  1  0 -1]
            [ 5 -1  1 -1 -1  1  0]
            [ 6  0 -2  0  0  0  1]
            [ 5  1  1 -1  1 -1  0]
            [ 4  2  0  1 -1  0 -1]
            [ 1  1  1  1  1  1  1]
        """
        G = self._semigroup
        B = self.basis()
        chi = [sum((g * B[k])[k] for k in B.keys())
               for g in G.conjugacy_classes_representatives()]
        try:
            return G.character(chi)
        except AttributeError:
            return vector(chi, immutable=True)

    @cached_method
    def brauer_character(self):
        r"""
        Return the Brauer character of ``self``.

        EXAMPLES::

            sage: SGA = SymmetricGroupAlgebra(GF(2), 5)
            sage: SM = SGA.specht_module([3, 2])
            sage: SM.brauer_character()
            (5, -1, 0)
            sage: SM.simple_module().brauer_character()
            (4, -2, -1)

            sage: T = SM.subrepresentation([])
            sage: T.brauer_character()
            (0, 0, 0)

            sage: W = CoxeterGroup(['D', 4], implementation='permutation')
            sage: R = W.reflection_representation(GF(2))
            sage: R.brauer_character()
            (4, 1)
            sage: T = R.subrepresentation([])
            sage: T.brauer_character()
            (0, 0)
        """
        G = self._semigroup
        p = self.base_ring().characteristic()

        if self.dimension() == 0:
            from sage.rings.rational_field import QQ
            ccrep = [g for g in G.conjugacy_classes_representatives()
                     if not p.divides(g.order())]
            return vector(QQ, [QQ.zero()] * len(ccrep))

        from sage.rings.number_field.number_field import CyclotomicField
        chi = []
        for g in G.conjugacy_classes_representatives():
            if p.divides(g.order()):
                # ignore the non-p-regular elements
                continue
            evals = self.representation_matrix(g).eigenvalues()
            K = evals[0].parent()
            val = 0
            orders = {la: la.multiplicative_order() for la in evals if la != K.one()}
            zetas = {o: CyclotomicField(o).gen() for o in orders.values()}
            prims = {o: K.zeta(o) for o in orders.values()}
            for la in evals:
                if la == K.one():
                    val += 1
                    continue
                o = la.multiplicative_order()
                zeta = zetas[o]
                prim = prims[o]
                for deg in range(o):
                    if prim ** deg == la:
                        val += zeta ** deg
                        break
            chi.append(val)

        return vector(chi, immutable=True)

    def exterior_power(self, degree=None):
        r"""
        Return the exterior power of ``self``.

        INPUT:

        - ``degree`` -- (optional) if given, then only consider the
          given degree

        EXAMPLES::

            sage: DC3 = groups.permutation.DiCyclic(3)
            sage: L = DC3.regular_representation(QQ, side='left')
            sage: E5 = L.exterior_power(5)
            sage: E5
            Exterior power representation of Left Regular Representation of
             Dicyclic group of order 12 as a permutation group over Rational Field
             in degree 5
            sage: L.exterior_power()
            Exterior algebra representation of Left Regular Representation of
             Dicyclic group of order 12 as a permutation group over Rational Field
        """
        if degree is None or degree == 0:
            return Representation_ExteriorAlgebra(self, degree)
        return Representation_Exterior(self, degree)

    def symmetric_power(self, degree=None):
        r"""
        Return the symmetric power of ``self`` in degree ``degree``.

        EXAMPLES::

            sage: W = CoxeterGroup(['H', 3])
            sage: R = W.reflection_representation()
            sage: S3R = R.symmetric_power(3)
            sage: S3R
            Symmetric power representation of Reflection representation of
            Finite Coxeter group over ... with Coxeter matrix:
            [1 3 2]
            [3 1 5]
            [2 5 1] in degree 3
        """
        return Representation_Symmetric(self, degree)

    def schur_functor(self, la):
        r"""
        Return the :class:`Schur functor
        <sage.modules.with_basis.representation.SchurFunctorRepresentation>`
        with shape ``la`` applied to ``self``.

        EXAMPLES::

            sage: W = CoxeterGroup(['H', 3])
            sage: R = W.reflection_representation()
            sage: S111 = R.schur_functor([1,1,1])
            sage: S111.dimension()
            1
            sage: S3 = R.schur_functor([3])
            sage: S3.dimension()
            10
        """
        return SchurFunctorRepresentation(self, la)

    @abstract_method(optional=True)
    def _semigroup_action(self, g, vec, vec_on_left):
        """
        Return the action of the semigroup element ``g`` on the
        vector ``vec`` of ``self``.

        If this is not defined, the representation element must
        override ``_acted_upon_``.

        EXAMPLES::

            sage: DC3 = groups.permutation.DiCyclic(3)
            sage: T = DC3.trivial_representation()
            sage: T._semigroup_action(DC3.an_element(), T.basis()['v'], True)
            B['v']
        """

    @cached_method
    def is_irreducible(self) -> bool:
        r"""
        Return if ``self`` is an irreducible module or not.

        A representation `M` is *irreducible* (also known as simple)
        if the only subrepresentations of `M` are the trivial module
        `\{0\}` and `M` itself.

        EXAMPLES::

            sage: DC3 = groups.permutation.DiCyclic(3)
            sage: L = DC3.regular_representation(GF(3), side='left')
            sage: L.is_irreducible()
            False
            sage: E3 = L.exterior_power(3)
            sage: E3.is_irreducible()
            False
            sage: E12 = L.exterior_power(12)
            sage: E12.is_irreducible()
            True

            sage: SGA = SymmetricGroupAlgebra(GF(2), 5)
            sage: SGA.specht_module([3, 2]).is_irreducible()
            False
            sage: SGA.specht_module([4, 1]).is_irreducible()
            True
        """
        return bool(self.find_subrepresentation() is None)

    def find_subrepresentation(self):
        r"""
        Return a nontrivial (not ``self`` or the trivial module) submodule
        of ``self`` or ``None`` if ``self`` is irreducible.

        EXAMPLES::

            sage: SGA = SymmetricGroupAlgebra(GF(2), 5)
            sage: SM = SGA.specht_module([3, 2])
            sage: U = SM.find_subrepresentation()
            sage: [SM(b) for b in U.basis()]
            [S[[1, 2, 3], [4, 5]] + S[[1, 2, 4], [3, 5]]
             + S[[1, 2, 5], [3, 4]] + S[[1, 3, 4], [2, 5]]]
            sage: B = U.basis()[0].lift()
            sage: all(g * B == B for g in SGA.gens())
            True
            sage: SGA.specht_module([4, 1]).find_subrepresentation() is None
            True
            sage: SGA.specht_module([5]).find_subrepresentation() is None
            True
        """
        if self.dimension() <= 1:  # trivially irreducible
            return None

        if not self.is_finite():
            raise NotImplementedError("only implemented for finite modules")

        if self._semigroup.cardinality() < 2:
            return self.subrepresentation([next(iter(self.basis()))], is_closed=True)

        R = self.base_ring()
        gens = [self.representation_matrix(g) for g in self._semigroup.gens()]
        try:
            # add the identity to the generators (if it exists)
            gens.append(self.representation_matrix(self._semigroup.one()))
        except (AttributeError, NotImplementedError, TypeError, ValueError):
            pass
        gens_transpose = tuple([g.transpose() for g in gens])
        amb_dim = self.dimension()

        def check_submodule(xi, G):
            N = xi.right_kernel_matrix()
            SM = N[:1]
            dim = 0
            while dim < SM.nrows() < amb_dim:
                dim = SM.nrows()
                added = matrix([g * vec for g in G for vec in SM.rows()])
                SM = SM.stack(added)
                SM.echelonize()
                SM = SM[:SM.rank()]
            if SM.nrows() < amb_dim:
                return SM
            return None

        def generate_elements():
            cur = sum(R.random_element() * g for g in gens)
            while True:
                if cur.is_zero():
                    cur = sum(R.random_element() * g for g in gens)
                yield cur
                cur *= sum(R.random_element() * g for g in gens)

        for theta in generate_elements():
            chi = theta.charpoly()
            factors = sorted(chi.factor(), key=lambda c: c[0].degree())
            for f, _ in factors:
                xi = f(theta)
                if not xi.is_singular():  # no nullspaces
                    continue
                SM = check_submodule(xi, gens)
                if SM is not None:
                    return self.subrepresentation([self.from_vector(v) for v in SM.rows()],
                                                  is_closed=True)
                SM = check_submodule(xi.transpose(), gens_transpose)
                if SM is not None:
                    # We instead want the submodule given by the orthogonal complement
                    return self.subrepresentation([self.from_vector(v) for v in SM.right_kernel_matrix().rows()],
                                                  is_closed=True)
                if xi.right_kernel_matrix().nrows() == f.degree():  # good factor
                    return None  # irreducible

    def subrepresentation(self, gens, check=True, already_echelonized=False,
                          *args, is_closed=False, **opts):
        """
        Construct a subrepresentation of ``self`` generated by ``gens``.

        INPUT:

        - ``gens`` -- the generators of the submodule
        - ``check`` -- ignored
        - ``already_echelonized`` -- (default: ``False``) whether
          the elements of ``gens`` are already in (not necessarily
          reduced) echelon form
        - ``is_closed`` -- (keyword only; default: ``False``) whether ``gens``
          already spans the subspace closed under the semigroup action

        EXAMPLES::

            sage: SGA = SymmetricGroupAlgebra(GF(2), 5)
            sage: SM = SGA.specht_module([3, 2])
            sage: B = next(iter(SM.basis()))
            sage: U = SM.subrepresentation([B])
            sage: U.dimension()
            5
        """
        if not is_closed and gens:
            R = self.base_ring()
            repr_mats = [self.representation_matrix(g) for g in self._semigroup.gens()]
            amb_dim = self.dimension()
            SM = matrix([v._vector_() for v in gens])
            SM.echelonize()
            SM = SM[:SM.rank()]
            dim = 0
            while dim < SM.nrows() < amb_dim:
                dim = SM.nrows()
                added = matrix([g * vec for g in repr_mats for vec in SM.rows()])
                SM = SM.stack(added)
                SM.echelonize()
                SM = SM[:SM.rank()]
            gens = [self.from_vector(v) for v in SM.rows()]
            # it might not be echelonized w.r.t. the module's basis ordering
            already_echelonized = False
        return self.submodule(gens, *args, submodule_class=Subrepresentation, check=check,
                              already_echelonized=already_echelonized, **opts)

    def quotient_representation(self, subrepr, already_echelonized=False, **kwds):
        r"""
        Construct a quotient representation of ``self`` by ``subrepr``.

        EXAMPLES::

            sage: SGA = SymmetricGroupAlgebra(GF(2), 5)
            sage: SM = SGA.specht_module([3, 2])
            sage: v = sum(list(SM.basis())[1:])
            sage: Q = SM.quotient_representation([v]); Q
            Quotient representation with basis {[[1, 3, 5], [2, 4]],
             [[1, 3, 4], [2, 5]], [[1, 2, 4], [3, 5]], [[1, 2, 3], [4, 5]]}
             of Specht module of [3, 2] over Finite Field of size 2
            sage: Q.is_irreducible()
            True
        """
        if not isinstance(subrepr, Subrepresentation):
            subrepr = self.subrepresentation(subrepr, unitriangular=True,
                                             already_echelonized=already_echelonized)
        return QuotientRepresentation(subrepr, **kwds)

    @cached_method
    def _composition_series_data(self):
        r"""
        Return a composition series and the simple quotients (i.e.,
        the composition factors) of ``self``.

        EXAMPLES:

        The algorithm used here uses random elements, so we set the seed
        for testing purposes::

            sage: set_random_seed(0)
            sage: G = groups.permutation.Cyclic(6)
            sage: R = G.regular_representation(GF(3))
            sage: CS, CF = R._composition_series_data()
            sage: [[R(b) for b in F.basis()] for F in CS]
            [[(),
              (1,2,3,4,5,6),
              (1,3,5)(2,4,6),
              (1,4)(2,5)(3,6),
              (1,5,3)(2,6,4),
              (1,6,5,4,3,2)],
             [() + 2*(1,6,5,4,3,2),
              (1,2,3,4,5,6) + 2*(1,6,5,4,3,2),
              (1,3,5)(2,4,6) + 2*(1,6,5,4,3,2),
              (1,4)(2,5)(3,6) + 2*(1,6,5,4,3,2),
              (1,5,3)(2,6,4) + 2*(1,6,5,4,3,2)],
             [() + (1,5,3)(2,6,4) + (1,6,5,4,3,2),
              (1,2,3,4,5,6) + 2*(1,5,3)(2,6,4),
              (1,3,5)(2,4,6) + 2*(1,6,5,4,3,2),
              (1,4)(2,5)(3,6) + (1,5,3)(2,6,4) + (1,6,5,4,3,2)],
             [() + 2*(1,4)(2,5)(3,6),
              (1,2,3,4,5,6) + 2*(1,5,3)(2,6,4),
              (1,3,5)(2,4,6) + 2*(1,6,5,4,3,2)],
             [() + 2*(1,3,5)(2,4,6) + 2*(1,4)(2,5)(3,6) + (1,6,5,4,3,2),
              (1,2,3,4,5,6) + (1,3,5)(2,4,6) + 2*(1,5,3)(2,6,4) + 2*(1,6,5,4,3,2)],
             [() + 2*(1,2,3,4,5,6) + (1,3,5)(2,4,6) + 2*(1,4)(2,5)(3,6)
              + (1,5,3)(2,6,4) + 2*(1,6,5,4,3,2)],
             []]
            sage: [F.dimension() for F in CF]
            [1, 1, 1, 1, 1, 1]
        """
        from sage.data_structures.blas_dict import linear_combination
        series = [self]
        cur = 0
        # The natural condition is ``while cur < len(series)``. However, the
        #   loop will always terminate from the break statement. So we skip
        #   this test for speed.
        while True:
            V = series[cur]
            if cur == len(series) - 1:
                W = V.find_subrepresentation()
                if W is None:  # V is irreducible
                    break
                # Construct W as a subrepresentation of ``self`` for consistency
                Wp = self.subrepresentation([self(b) for b in W._basis],
                                            already_echelonized=True, is_closed=True)
                series.append(Wp)
            else:
                W = V.subrepresentation([V.retract(b) for b in series[cur+1]._basis],
                                        already_echelonized=True, is_closed=True)

            Q = V.quotient_representation(W)
            S = Q.find_subrepresentation()
            while S is not None:  # found a nontrivial subrepresentation
                # Lift S -> Q -> V -> self
                # S_basis = [b.lift().lift().lift() for b in S.basis()]
                # Fast version not creating transient elements
                if V is self:
                    S_basis = [self.element_class(self, b._monomial_coefficients) for b in S._basis]
                else:
                    S_basis = [self.element_class(self, linear_combination((V._basis[i]._monomial_coefficients, coeff)
                                                                           for i, coeff in b._monomial_coefficients.items()))
                               for b in S._basis]
                # Lift the basis of W', which is W as a subrepresentation of ``self``
                # This is equivalent to [b.lift() for b in series[cur+1].basis()]
                Wp_basis = list(series[cur+1]._basis)
                Wp = self.subrepresentation(S_basis + Wp_basis, is_closed=True)
                series.insert(cur+1, Wp)
                if V is self:
                    W = Wp
                else:
                    W = V.subrepresentation([V.retract(b) for b in series[cur+1]._basis],
                                            already_echelonized=True, is_closed=True)
                Q = V.quotient_representation(W)
                S = Q.find_subrepresentation()

            cur += 1

        # Special case when ``self`` is irreducible
        if len(series) == 1:
            simples = (self,)
            series.append(self.subrepresentation([], is_closed=True))
            return (tuple(series), simples)

        # Convert series to a tower and include the 0 dimensional
        #   representation at the end.
        ret = [self, series[1]]
        prev = series[1]
        lift = prev.lift
        retract = lift.section()
        for W in series[2:]:
            prev = prev.subrepresentation([retract(b) for b in W._basis],
                                          already_echelonized=True, is_closed=True)
            ret.append(prev)

            # Construct the lift map prev -> self
            data = {i: lift(prev._basis[i]) for i in prev._basis.keys()}
            lift = prev.module_morphism(data.__getitem__,
                                        codomain=self,
                                        triangular='lower',
                                        unitriangular=False,
                                        key=W._support_key,
                                        inverse_on_support='compute')
            retract = lift.section()

        ret.append(ret[-1].subrepresentation([], is_closed=True))

        # Construct the simples
        simples = [ret[i].quotient_representation(ret[i+1])
                   for i in range(len(ret)-2)]
        simples.append(ret[-2])

        return (tuple(ret), tuple(simples))

    def composition_series(self):
        r"""
        Return a composition series of ``self``.

        EXAMPLES:

        The algorithm used here uses random elements, so we set the seed
        for testing purposes::

            sage: set_random_seed(0)
            sage: G = groups.permutation.Dihedral(5)
            sage: CFM = CombinatorialFreeModule(GF(2), [1, 2, 3, 4, 5],
            ....:                               bracket=False, prefix='e')
            sage: CFM.an_element()
            e3
            sage: R = G.representation(CFM, lambda g, i: CFM.basis()[g(i)], side='right')
            sage: CS = R.composition_series()
            sage: len(CS)
            3
            sage: [[R(b) for b in F.basis()] for F in CS]
            [[e1, e2, e3, e4, e5], [e1 + e5, e2 + e5, e3 + e5, e4 + e5], []]
            sage: [F.brauer_character() for F in CS]
            [(5, 0, 0), (4, -1, -1), (0, 0, 0)]
            sage: [F.brauer_character() for F in R.composition_factors()]
            [(1, 1, 1), (4, -1, -1)]
            sage: Reg = G.regular_representation(GF(2))
            sage: simple_brauer_chars = set([F.brauer_character()
            ....:                            for F in Reg.composition_factors()])
            sage: sorted(simple_brauer_chars)
            [(1, 1, 1), (4, -1, -1)]
        """
        return self._composition_series_data()[0]

    def composition_factors(self):
        r"""
        Return the composition factors of ``self``.

        Given a composition series `V = V_0 \subseteq V_1 \subseteq \cdots
        \subseteq V_{\ell} = 0`, the composition factor `S_i` is defined as
        `V_i / V_{i+1}`.

        EXAMPLES:

        The algorithm used here uses random elements, so we set the seed
        for testing purposes::

            sage: set_random_seed(0)
            sage: SGA = SymmetricGroupAlgebra(GF(3), 6)
            sage: SM = SGA.specht_module([4, 1, 1])
            sage: CF = SM.composition_factors()
            sage: CF
            (Quotient representation with basis
              {[[1, 2, 5, 6], [3], [4]], [[1, 2, 4, 6], [3], [5]], [[1, 2, 3, 6], [4], [5]],
               [[1, 2, 4, 5], [3], [6]], [[1, 2, 3, 5], [4], [6]], [[1, 2, 3, 4], [5], [6]]}
              of Specht module of [4, 1, 1] over Finite Field of size 3,
             Subrepresentation with basis {0, 1, 2, 3} of Specht module
              of [4, 1, 1] over Finite Field of size 3)
            sage: x = SGA.an_element()
            sage: v = CF[1].an_element(); v
            2*B[0] + 2*B[1]
            sage: x * v
            B[1] + B[2]

        We reproduce the decomposition matrix for `S_5` over `\GF{2}`::

            sage: SGA = SymmetricGroupAlgebra(GF(2), 5)
            sage: simples = [SGA.simple_module(la).brauer_character()
            ....:            for la in SGA.simple_module_parameterization()]
            sage: D = []
            sage: for la in Partitions(5):
            ....:     SM = SGA.specht_module(la)
            ....:     data = [CF.brauer_character() for CF in SM.composition_factors()]
            ....:     D.append([data.count(bc) for bc in simples])
            sage: matrix(D)
            [1 0 0]
            [0 1 0]
            [1 0 1]
            [2 0 1]
            [1 0 1]
            [0 1 0]
            [1 0 0]
        """
        return self._composition_series_data()[1]

    class Element(CombinatorialFreeModule.Element):
        def _acted_upon_(self, scalar, self_on_left=False):
            """
            Return the action of ``scalar`` on ``self``.

            EXAMPLES::

                sage: G = groups.misc.WeylGroup(['B',2], prefix='s')
                sage: R = G.regular_representation()
                sage: s1,s2 = G.gens()
                sage: x = R.an_element(); x
                2*s2*s1*s2 + s1*s2 + 3*s2 + 1
                sage: 2 * x
                4*s2*s1*s2 + 2*s1*s2 + 6*s2 + 2
                sage: s1 * x
                2*s2*s1*s2*s1 + 3*s1*s2 + s1 + s2
                sage: s2 * x
                s2*s1*s2 + 2*s1*s2 + s2 + 3

                sage: G = groups.misc.WeylGroup(['B',2], prefix='s')
                sage: R = G.regular_representation(side='right')
                sage: s1,s2 = G.gens()
                sage: x = R.an_element(); x
                2*s2*s1*s2 + s1*s2 + 3*s2 + 1
                sage: x * s1
                2*s2*s1*s2*s1 + s1*s2*s1 + 3*s2*s1 + s1
                sage: x * s2
                2*s2*s1 + s1 + s2 + 3

                sage: G = groups.misc.WeylGroup(['B',2], prefix='s')
                sage: R = G.regular_representation()
                sage: R.base_ring()
                Integer Ring
                sage: A = G.algebra(ZZ)
                sage: s1,s2 = A.algebra_generators()
                sage: x = R.an_element(); x
                2*s2*s1*s2 + s1*s2 + 3*s2 + 1
                sage: s1 * x
                2*s2*s1*s2*s1 + 3*s1*s2 + s1 + s2
                sage: s2 * x
                s2*s1*s2 + 2*s1*s2 + s2 + 3
                sage: (2*s1 - s2) * x
                4*s2*s1*s2*s1 - s2*s1*s2 + 4*s1*s2 + 2*s1 + s2 - 3
                sage: (3*s1 + s2) * R.zero()
                0

                sage: A = G.algebra(QQ)
                sage: s1,s2 = A.algebra_generators()
                sage: a = 1/2 * s1
                sage: a * x
                Traceback (most recent call last):
                ...
                TypeError: unsupported operand parent(s) for *:
                 'Algebra of Weyl Group of type ['B', 2] ... over Rational Field'
                 and 'Left Regular Representation of Weyl Group of type ['B', 2] ... over Integer Ring'

            Check that things that coerce into the group (algebra) also have
            an action::

                sage: D4 = groups.permutation.Dihedral(4)
                sage: S4 = SymmetricGroup(4)
                sage: S4.has_coerce_map_from(D4)
                True
                sage: R = S4.regular_representation()
                sage: D4.an_element() * R.an_element()
                2*(2,4) + 3*(1,2,3,4) + (1,3) + (1,4,2,3)
            """
            if isinstance(scalar, Element):
                P = self.parent()
                sP = scalar.parent()
                if sP is P._semigroup:
                    if not self:
                        return self
                    return P._semigroup_action(scalar, self, self_on_left)

                if sP is P._semigroup_algebra:
                    if not self:
                        return self
                    return P.linear_combination(((P._semigroup_action(ms, self, self_on_left), cs)
                                                 for ms, cs in scalar), not self_on_left)

                if P._semigroup.has_coerce_map_from(sP):
                    scalar = P._semigroup(scalar)
                    return self._acted_upon_(scalar, self_on_left)

                # Check for scalars first before general coercion to the semigroup algebra.
                # This will result in a faster action for the scalars.
                ret = CombinatorialFreeModule.Element._acted_upon_(self, scalar, self_on_left)
                if ret is not None:
                    return ret

                if P._semigroup_algebra.has_coerce_map_from(sP):
                    scalar = P._semigroup_algebra(scalar)
                    return self._acted_upon_(scalar, self_on_left)

                return None

            return CombinatorialFreeModule.Element._acted_upon_(self, scalar, self_on_left)


class Representation(Representation_abstract, CombinatorialFreeModule):
    r"""
    Representation of a semigroup.

    INPUT:

    - ``semigroup`` -- a semigroup
    - ``module`` -- a module with a basis
    - ``on_basis`` -- function which takes as input ``g``, ``m``, where
      ``g`` is an element of the semigroup and ``m`` is an element of the
      indexing set for the basis, and returns the result of ``g`` acting
      on ``m``
    - ``side`` -- (default: ``'left'``) whether this is a
      ``'left'`` or ``'right'`` representation

    EXAMPLES:

    We construct the sign representation of a symmetric group::

        sage: G = SymmetricGroup(4)
        sage: M = CombinatorialFreeModule(QQ, ['v'])
        sage: on_basis = lambda g,m: M.term(m, g.sign())
        sage: R = G.representation(M, on_basis)
        sage: x = R.an_element(); x
        2*B['v']
        sage: c,s = G.gens()
        sage: c,s
        ((1,2,3,4), (1,2))
        sage: c * x
        -2*B['v']
        sage: s * x
        -2*B['v']
        sage: c * s * x
        2*B['v']
        sage: (c * s) * x
        2*B['v']

    This extends naturally to the corresponding group algebra::

        sage: A = G.algebra(QQ)
        sage: s,c = A.algebra_generators()
        sage: c,s
        ((1,2,3,4), (1,2))
        sage: c * x
        -2*B['v']
        sage: s * x
        -2*B['v']
        sage: c * s * x
        2*B['v']
        sage: (c * s) * x
        2*B['v']
        sage: (c + s) * x
        -4*B['v']

    REFERENCES:

    - :wikipedia:`Group_representation`
    """
    def __init__(self, semigroup, module, on_basis, side='left', **kwargs):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: G = SymmetricGroup(4)
            sage: M = CombinatorialFreeModule(QQ, ['v'])
            sage: def on_basis(g, m):
            ....:     return M.term(m, g.sign())
            sage: R = G.representation(M, on_basis)
            sage: R._test_representation()

            sage: G = CyclicPermutationGroup(3)
            sage: M = algebras.Exterior(QQ, 'x', 3)
            sage: def on_basis(g, m):  # cyclically permute generators
            ....:     return M.prod([M.monomial(FrozenBitset([g(j+1)-1])) for j in m])
            sage: from sage.categories.algebras import Algebras
            sage: R = G.representation(M, on_basis, category=Algebras(QQ).WithBasis().FiniteDimensional())
            sage: r = R.an_element(); r
            1 + 2*x0 + x0*x1 + 3*x1
            sage: r*r
            1 + 4*x0 + 2*x0*x1 + 6*x1
            sage: x0, x1, x2 = M.gens()
            sage: s = R(x0*x1)
            sage: g = G.an_element()
            sage: g*s
            x1*x2
            sage: g*R(x1*x2)
            -x0*x2
            sage: g*r
            1 + 2*x1 + x1*x2 + 3*x2
            sage: g^2*r
            1 + 3*x0 - x0*x2 + 2*x2

            sage: G = SymmetricGroup(4)
            sage: A = SymmetricGroup(4).algebra(QQ)
            sage: def action(g, x):
            ....:     return A.monomial(g*x)
            sage: category = Algebras(QQ).WithBasis().FiniteDimensional()
            sage: R = G.representation(A, action, 'left', category=category)
            sage: r = R.an_element(); r
            () + (2,3,4) + 2*(1,3)(2,4) + 3*(1,4)(2,3)
            sage: r^2
            14*() + 2*(2,3,4) + (2,4,3) + 12*(1,2)(3,4) + 3*(1,2,4) + 2*(1,3,2) + 4*(1,3)(2,4) + 5*(1,4,3) + 6*(1,4)(2,3)
            sage: g = G.an_element(); g
            (2,3,4)
            sage: g*r
            (2,3,4) + (2,4,3) + 2*(1,3,2) + 3*(1,4,3)
        """
        try:
            self.product_on_basis = module.product_on_basis
        except AttributeError:
            pass

        category = kwargs.pop('category', Modules(module.base_ring()).WithBasis())

        self._on_basis = on_basis
        self._module = module
        if side == "twosided":
            raise ValueError("the defined action must be either left or right")

        indices = module.basis().keys()

        if 'FiniteDimensional' in module.category().axioms():
            category = category.FiniteDimensional()

        CombinatorialFreeModule.__init__(self, module.base_ring(), indices, category=category,
                                         **module.print_options())
        Representation_abstract.__init__(self, semigroup, side)

    def _test_representation(self, **options):
        """
        Check (on some elements) that ``self`` is a representation of the
        given semigroup.

        EXAMPLES::

            sage: G = groups.permutation.Dihedral(4)
            sage: R = G.regular_representation()
            sage: R._test_representation()

            sage: G = CoxeterGroup(['A',4,1], base_ring=ZZ)
            sage: M = CombinatorialFreeModule(QQ, ['v'])
            sage: def on_basis(g, m):
            ....:     return M.term(m, (-1)**g.length())
            sage: R = G.representation(M, on_basis, side='right')
            sage: R._test_representation(max_runs=500)
        """
        from sage.misc.functional import sqrt
        tester = self._tester(**options)
        S = tester.some_elements()
        L = []
        max_len = int(sqrt(tester._max_runs)) + 1
        for i, x in enumerate(self._semigroup):
            L.append(x)
            if i >= max_len:
                break
        for x in L:
            for y in L:
                for elt in S:
                    if self._left_repr:
                        tester.assertEqual(x*(y*elt), (x*y)*elt)
                    else:
                        tester.assertEqual((elt*y)*x, elt*(y*x))

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: P = Permutations(4)
            sage: M = CombinatorialFreeModule(QQ, ['v'])
            sage: on_basis = lambda g,m: M.term(m, g.sign())
            sage: P.representation(M, on_basis)
            Representation of Standard permutations of 4 indexed by {'v'}
             over Rational Field
        """
        return "Representation of {} indexed by {} over {}".format(
            self._semigroup, self.basis().keys(), self.base_ring())

    def _repr_term(self, b):
        """
        Return a string representation of a basis index ``b`` of ``self``.

        EXAMPLES::

            sage: SGA = SymmetricGroupAlgebra(QQ, 3)
            sage: R = SGA.regular_representation()
            sage: all(R._repr_term(b) == SGA._repr_term(b) for b in SGA.basis().keys())
            True
        """
        return self._module._repr_term(b)

    def _latex_term(self, b):
        """
        Return a LaTeX representation of a basis index ``b`` of ``self``.

        EXAMPLES::

            sage: SGA = SymmetricGroupAlgebra(QQ, 3)
            sage: R = SGA.regular_representation()
            sage: all(R._latex_term(b) == SGA._latex_term(b) for b in SGA.basis().keys())
            True
        """
        return self._module._latex_term(b)

    def _element_constructor_(self, x):
        """
        Construct an element of ``self`` from ``x``.

        EXAMPLES::

            sage: G = groups.permutation.Dihedral(4)
            sage: A = G.algebra(ZZ)
            sage: R = A.regular_representation()
            sage: x = A.an_element(); x
            () + (1,3) + 2*(1,3)(2,4) + 3*(1,4,3,2)
            sage: R(x)
            () + (1,3) + 2*(1,3)(2,4) + 3*(1,4,3,2)
        """
        if isinstance(x, Element) and x.parent() is self._module:
            return self._from_dict(x.monomial_coefficients(copy=False), remove_zeros=False)
        return super()._element_constructor_(x)

    def product_by_coercion(self, left, right):
        r"""
        Return the product of ``left`` and ``right`` by passing to
        ``self._module`` and then building a new element of ``self``.

        EXAMPLES::

            sage: G = groups.permutation.KleinFour()
            sage: E = algebras.Exterior(QQ,'e',4)
            sage: on_basis = lambda g,m: E.monomial(m) # the trivial representation
            sage: R = G.representation(E, on_basis)
            sage: r = R.an_element(); r
            1 + 2*e0 + 3*e1 + e1*e2
            sage: g = G.an_element();
            sage: g * r == r  # indirect doctest
            True
            sage: r * r  # indirect doctest
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for *:
             'Representation of The Klein 4 group of order 4, as a permutation
             group indexed by Subsets of {0,1,...,3} over Rational Field' and
             'Representation of The Klein 4 group of order 4, as a permutation
             group indexed by Subsets of {0,1,...,3} over Rational Field'

            sage: from sage.categories.algebras import Algebras
            sage: category = Algebras(QQ).FiniteDimensional().WithBasis()
            sage: T = G.representation(E, on_basis, category=category)
            sage: t = T.an_element(); t
            1 + 2*e0 + 3*e1 + e1*e2
            sage: g * t == t  # indirect doctest
            True
            sage: t * t  # indirect doctest
            1 + 4*e0 + 4*e0*e1*e2 + 6*e1 + 2*e1*e2
        """
        M = self._module

        # Multiply in self._module
        p = M._from_dict(left._monomial_coefficients, False, False) * M._from_dict(right._monomial_coefficients, False, False)

        # Convert from a term in self._module to a term in self
        return self._from_dict(p.monomial_coefficients(copy=False), False, False)

    def _semigroup_action(self, g, vec, vec_on_left):
        """
        Return the action of the semigroup element ``g`` on the
        vector ``vec`` of ``self``.

        EXAMPLES::

            sage: G = groups.permutation.KleinFour()
            sage: E = algebras.Exterior(QQ,'e',4)
            sage: on_basis = lambda g,m: E.monomial(m) # the trivial representation
            sage: R = G.representation(E, on_basis)
            sage: R._semigroup_action(G.an_element(), R.an_element(), True)
            1 + 2*e0 + 3*e1 + e1*e2
        """
        if self._left_repr == vec_on_left:
            g = ~g
        return self.linear_combination(((self._on_basis(g, m), c)
                                       for m, c in vec._monomial_coefficients.items()), not vec_on_left)


class Subrepresentation(Representation_abstract, SubmoduleWithBasis):
    r"""
    A subrepresentation.

    Let `R` be a representation of an algebraic object `X`. A
    subrepresentation is a submodule of `R` that is closed under
    the action of `X`.
    """
    # Use the same normalization as the base class
    __classcall_private__ = SubmoduleWithBasis.__classcall_private__

    def __init__(self, basis, support_order, ambient, *args, **opts):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: G = groups.permutation.Dihedral(4)
            sage: R = G.regular_representation(QQ)
            sage: S = R.subrepresentation([sum(R.basis())])
            sage: TestSuite(S).run()
        """
        SubmoduleWithBasis.__init__(self, basis, support_order, ambient, *args, **opts)
        Representation_abstract.__init__(self, ambient.semigroup(), ambient.side(), ambient.semigroup_algebra())

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: G = groups.permutation.Dihedral(4)
            sage: R = G.regular_representation()
            sage: R.subrepresentation([sum(R.basis())], is_closed=True)
            Subrepresentation with basis {0} of Left Regular Representation of
             Dihedral group of order 8 as a permutation group over Integer Ring
        """
        return "Subrepresentation with basis {} of {}".format(self.basis().keys(), self._ambient)

    class Element(SubmoduleWithBasis.Element):
        def _acted_upon_(self, scalar, self_on_left=True):
            """
            Return the action of ``scalar`` on ``self``.

            EXAMPLES::

                sage: G = groups.permutation.Dihedral(4)
                sage: CFM = CombinatorialFreeModule(GF(2), [1, 2, 3, 4],
                ....:                               bracket=False, prefix='e')
                sage: R = G.representation(CFM, lambda g, i: CFM.basis()[g(i)], side='right')
                sage: e1, e2, e3, e4 = R.basis()
                sage: S = R.subrepresentation([e1 + e3, e2 + e4], is_closed=True)
                sage: x = G.an_element(); x
                (1,3)
                sage: v = sum(S.basis()); v
                B[0] + B[1]
                sage: x * v
                B[0] + B[1]
                sage: [x * b for b in S.basis()]
                [B[0], B[1]]
                sage: [[g * b for g in G] for b in S.basis()]
                [[B[0], B[0], B[1], B[1], B[0], B[0], B[1], B[1]],
                 [B[1], B[1], B[0], B[0], B[1], B[1], B[0], B[0]]]
                sage: 2 * v
                0

                sage: Q = R.quotient_representation(S)
                sage: [[g * b for g in G] for b in Q.basis()]
                [[B[3], B[3], B[4], B[4], B[3], B[3], B[4], B[4]],
                 [B[4], B[4], B[3], B[3], B[4], B[4], B[3], B[3]]]
            """
            ret = super()._acted_upon_(scalar, self_on_left)
            if ret is not None:
                return ret
            P = self.parent()
            if scalar in P._semigroup or scalar in P._semigroup_algebra:
                return P.retract(self.lift()._acted_upon_(scalar, self_on_left))
            return None


class QuotientRepresentation(Representation_abstract, QuotientModuleWithBasis):
    """
    The quotient of a representation by another representation, which
    admits a natural structure of a representation.
    """
    # Use the same normalization as the base class
    __classcall_private__ = QuotientModuleWithBasis.__classcall_private__

    def __init__(self, *args, **kwds):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: G = groups.permutation.Dihedral(4)
            sage: R = G.regular_representation(QQ)
            sage: S = R.subrepresentation([sum(R.basis())], is_closed=True)
            sage: Q = R.quotient_representation(S)
            sage: TestSuite(Q).run()
        """
        QuotientModuleWithBasis.__init__(self, *args, **kwds)
        amb = self.ambient()
        Representation_abstract.__init__(self, amb.semigroup(), amb.side(), amb.semigroup_algebra())

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: G = groups.permutation.Dihedral(4)
            sage: R = G.regular_representation(QQ)
            sage: S = R.subrepresentation([sum(R.basis())], is_closed=True)
            sage: R.quotient_representation(S)
            Quotient representation with basis
             {(1,3)(2,4), (1,4,3,2), (1,2,3,4), (2,4), (1,3), (1,4)(2,3), (1,2)(3,4)}
             of Left Regular Representation of Dihedral group of order 8
             as a permutation group over Rational Field
        """
        return "Quotient representation with basis {} of {}".format(self.basis().keys(), self._ambient)

    Element = Subrepresentation.Element


class Representation_Tensor(Representation_abstract, CombinatorialFreeModule_Tensor):
    r"""
    Tensor product of representations.
    """
    @staticmethod
    def __classcall_private__(cls, reps, **options):
        r"""
        Normalize input to ensure a unique representation.

        EXAMPLES::

            sage: S3 = SymmetricGroup(3)
            sage: L = S3.regular_representation(side='left')
            sage: S = S3.sign_representation()
            sage: R = S3.regular_representation(side='right')
            sage: tensor([tensor([L, S]), R]) == tensor([L, S, R])
            True
            sage: tensor([L, tensor([S, R])]) == tensor([L, S, R])
            True

        Check that the tensor product with more general modules
        can be constructed::

            sage: C = CombinatorialFreeModule(ZZ, ['a','b'])
            sage: T = tensor([C, R])
            sage: type(T)
            <class 'sage.combinat.free_module.CombinatorialFreeModule_Tensor_with_category'>
            sage: T = tensor([R, C])
            sage: type(T)
            <class 'sage.combinat.free_module.CombinatorialFreeModule_Tensor_with_category'>
        """
        assert len(reps) > 0
        assert isinstance(reps[0], Representation_abstract)
        S = reps[0].semigroup()
        if not all(isinstance(module, Representation_abstract)
                   and module.semigroup() == S for module in reps):
            return CombinatorialFreeModule_Tensor(reps, **options)
        R = reps[0].base_ring()
        if not all(module in Modules(R).WithBasis() for module in reps):
            raise ValueError("not all representations over the same base ring")
        # flatten the list of modules so that tensor(A, tensor(B,C)) gets rewritten into tensor(A, B, C)
        reps = sum((module._sets if isinstance(module, Representation_Tensor) else (module,) for module in reps), ())
        if all('FiniteDimensional' in M.category().axioms() for M in reps):
            options['category'] = options['category'].FiniteDimensional()
        return super().__classcall__(cls, reps, **options)

    def __init__(self, reps, **options):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: G = groups.permutation.Alternating(5)
            sage: L = G.regular_representation(side='left')
            sage: S = G.sign_representation()
            sage: T = tensor([L, S, L])
            sage: TestSuite(T).run()
        """
        sides = set(M.side() for M in reps)
        if "left" and "right" in sides:
            side = reps[0].side()  # make a choice as this is not fundamentally important
        else:
            if len(sides) == 2:  # mix of one side and twosided
                sides.remove("twosided")
            side, = sides  # get the unique side remaining
        CombinatorialFreeModule_Tensor.__init__(self, reps, **options)
        Representation_abstract.__init__(self, reps[0].semigroup(), side)

    def _semigroup_action(self, g, vec, vec_on_left):
        """
        Return the action of the semigroup element ``g`` on the
        vector ``vec`` of ``self``.

        EXAMPLES::

            sage: S3 = SymmetricGroup(3)
            sage: L = S3.regular_representation(side='left')
            sage: R = S3.regular_representation(side='right')
            sage: T = tensor([R, L])
            sage: g = S3.an_element(); g
            (2,3)
            sage: v = T.an_element(); v
            2*() # () + 3*() # (1,2,3) + 2*() # (1,3,2)
            sage: g * v
            2*(2,3) # (2,3) + 3*(2,3) # (1,2) + 2*(2,3) # (1,3)
            sage: T._semigroup_action(g, v, True)
            2*(2,3) # (2,3) + 3*(2,3) # (1,2) + 2*(2,3) # (1,3)
        """
        bases = [M.basis() for M in self._sets]
        if vec_on_left:
            return self.linear_combination((self._tensor_of_elements([B[k] * g for B, k in zip(bases, b)]), c)
                                           for b, c in vec._monomial_coefficients.items())
        return self.linear_combination((self._tensor_of_elements([g * B[k] for B, k in zip(bases, b)]), c)
                                       for b, c in vec._monomial_coefficients.items())

    class Element(Representation_abstract.Element):
        pass


Representation_abstract.Tensor = Representation_Tensor


class Representation_Exterior(Representation_abstract, CombinatorialFreeModule):
    r"""
    The exterior power representation (in a fixed degree).
    """
    def __init__(self, rep, degree=None, category=None, **options):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: G = groups.matrix.GL(3, 2)
            sage: R = G.regular_representation(side='right')
            sage: E2 = R.exterior_power(2)
            sage: E2.category()
            Category of finite dimensional modules with basis over Integer Ring
            sage: TestSuite(E2).run()

            sage: G = groups.matrix.GL(2, 3)
            sage: L = G.regular_representation(side='left')
            sage: E48 = L.exterior_power(48)
            sage: TestSuite(E48).run()

            sage: L.exterior_power(-2)
            Traceback (most recent call last):
            ...
            ValueError: the degree must be an integer in [0, 48]
            sage: L.exterior_power(120)
            Traceback (most recent call last):
            ...
            ValueError: the degree must be an integer in [0, 48]
            sage: L.exterior_power(5/6)
            Traceback (most recent call last):
            ...
            ValueError: the degree must be an integer in [0, 48]
        """
        from sage.algebras.clifford_algebra import ExteriorAlgebra
        from sage.algebras.clifford_algebra import CliffordAlgebraIndices
        from sage.rings.integer_ring import ZZ
        self._degree = degree
        self._rep = rep
        R = rep.base_ring()
        dim = rep.dimension()
        if degree is not None and (degree not in ZZ or degree > dim or degree < 0):
            raise ValueError(f"the degree must be an integer in [0, {dim}]")
        self._extalg = ExteriorAlgebra(R, dim)
        self._basis_order = list(rep.basis().keys())
        self._inv_map = {b: i for i, b in enumerate(self._basis_order)}
        ind = CliffordAlgebraIndices(dim, degree)
        R = rep.base_ring()
        category = Modules(R).WithBasis().or_subcategory(category)
        CombinatorialFreeModule.__init__(self, R, ind, category=category, **options)
        Representation_abstract.__init__(self, rep.semigroup(), rep.side(), rep.semigroup_algebra())

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: DC3 = groups.permutation.DiCyclic(3)
            sage: L = DC3.regular_representation(QQ, side='left')
            sage: L.exterior_power(7)
            Exterior power representation of Left Regular Representation of
             Dicyclic group of order 12 as a permutation group over Rational Field
             in degree 7
            sage: L.exterior_power()
            Exterior algebra representation of Left Regular Representation of
             Dicyclic group of order 12 as a permutation group over Rational Field
        """
        if self._degree is None:
            return "Exterior algebra representation of {}".format(repr(self._rep))
        return "Exterior power representation of {} in degree {}".format(repr(self._rep), self._degree)

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: DC3 = groups.permutation.DiCyclic(3)
            sage: L = DC3.regular_representation(QQ, side='left')
            sage: latex(L.exterior_power(4))
            \bigwedge^{4} ...
            sage: latex(L.exterior_power())
            \bigwedge ...
        """
        from sage.misc.latex import latex
        if self._degree is None:
            return "\\bigwedge " + latex(self._rep)
        return "\\bigwedge^{{{}}} ".format(self._degree) + latex(self._rep)

    def _repr_term(self, m):
        r"""
        Return a string representation of the basis element indexed by
        ``m``.

        EXAMPLES::

            sage: DC3 = groups.permutation.DiCyclic(3)
            sage: L = DC3.regular_representation(QQ, side='left')
            sage: E2 = L.exterior_power(2)
            sage: E2._repr_term(E2.an_element().leading_support())
            '()*(5,6,7)'
        """
        if len(m) == 0:
            return '1'
        B = self._rep.basis()
        return '*'.join(repr(B[self._basis_order[i]]) for i in m)

    def _ascii_art_term(self, m):
        r"""
        Return ascii art for the basis element indexed by ``m``.

        EXAMPLES::

            sage: DC3 = groups.permutation.DiCyclic(3)
            sage: L = DC3.regular_representation(QQ, side='left')
            sage: E2 = L.exterior_power(2)
            sage: E2._ascii_art_term(E2.an_element().leading_support())
            ()/\(5,6,7)
            sage: ascii_art(E2.an_element())
            2*()/\(5,6,7) + 2*()/\(5,7,6) + 3*()/\(1,2)(3,4)
        """
        from sage.typeset.ascii_art import ascii_art
        if len(m) == 0:
            return ascii_art('1')
        wedge = '/\\'
        B = self._rep.basis()
        return ascii_art(*[B[self._basis_order[i]] for i in m], sep=wedge)

    def _unicode_art_term(self, m):
        r"""
        Return unicode art for the basis element indexed by ``m``.

        EXAMPLES::

            sage: DC3 = groups.permutation.DiCyclic(3)
            sage: L = DC3.regular_representation(QQ, side='left')
            sage: E2 = L.exterior_power(2)
            sage: E2._unicode_art_term(E2.an_element().leading_support())
            ()∧(5,6,7)
            sage: unicode_art(E2.an_element())
            2*()∧(5,6,7) + 2*()∧(5,7,6) + 3*()∧(1,2)(3,4)
        """
        from sage.typeset.unicode_art import unicode_art
        if len(m) == 0:
            return unicode_art('1')
        import unicodedata
        wedge = unicodedata.lookup('LOGICAL AND')
        B = self._rep.basis()
        return unicode_art(*[B[self._basis_order[i]] for i in m], sep=wedge)

    def _latex_term(self, m):
        r"""
        Return a `\LaTeX` representation of the basis element indexed
        by ``m``.

        EXAMPLES::

            sage: DC3 = groups.permutation.DiCyclic(3)
            sage: L = DC3.regular_representation(QQ, side='left')
            sage: E2 = L.exterior_power(2)
            sage: E2._latex_term(E2.an_element().leading_support())
            '1 \\wedge (5,6,7)'
        """
        if len(m) == 0:
            return '1'
        from sage.misc.latex import latex
        B = self._rep.basis()
        return " \\wedge ".join(latex(B[self._basis_order[i]]) for i in m)

    def _from_repr_to_ext(self, elt):
        r"""
        Return the element ``elt`` from the defining representation
        to the corresponding exterior algebra.

        EXAMPLES::

            sage: G = groups.matrix.GL(2, 2)
            sage: L = G.regular_representation(side='left')
            sage: E = L.exterior_power()
            sage: E._from_repr_to_ext(sum(i*b for i,b in enumerate(L.basis(), start=1)))
            e0 + 2*e1 + 3*e2 + 4*e3 + 5*e4 + 6*e5
        """
        ind = self._indices
        data = {ind([self._inv_map[k]]): c for k, c in elt._monomial_coefficients.items()}
        return self._extalg.element_class(self._extalg, data)

    def _semigroup_action(self, g, vec, vec_on_left):
        r"""
        Return the action of the semigroup element ``g`` on the
        vector ``vec`` of ``self``.

        EXAMPLES::

            sage: DC3 = groups.permutation.DiCyclic(3)
            sage: g = DC3.an_element(); g
            (1,4,2,3)(5,6)
            sage: R = DC3.regular_representation(side='right')
            sage: E2 = R.exterior_power(2)
            sage: vec = E2.an_element(); vec
            2*()*(5,6,7) + 2*()*(5,7,6) + 3*()*(1,2)(3,4)
            sage: E2._semigroup_action(g, vec, True)
            -2*(1,4,2,3)(6,7)*(1,4,2,3)(5,6) + 2*(1,4,2,3)(5,6)*(1,4,2,3)(5,7)
             + 3*(1,4,2,3)(5,6)*(1,3,2,4)(5,6)
            sage: E2._semigroup_action(g, vec, False)
            -2*(1,3,2,4)(6,7)*(1,3,2,4)(5,6) + 2*(1,3,2,4)(5,6)*(1,3,2,4)(5,7)
             - 3*(1,4,2,3)(5,6)*(1,3,2,4)(5,6)
        """
        return self.linear_combination(((self._action_on_basis(g, b, vec_on_left), c)
                                        for b, c in vec._monomial_coefficients.items()), not vec_on_left)

    def _action_on_basis(self, g, b, vec_on_left):
        r"""
        Return the action of ``g`` on the basis element indexed by ``b``.

        EXAMPLES::

            sage: S3 = SymmetricGroup(3)
            sage: g = S3.an_element(); g
            (2,3)
            sage: L = S3.regular_representation(side='left')
            sage: E2 = L.exterior_power(2)
            sage: vec = E2.an_element(); vec
            2*()*(1,3,2) + 2*()*(1,2,3) + 3*()*(2,3)
            sage: g * vec
            2*(2,3)*(1,3) + 2*(2,3)*(1,2) - 3*()*(2,3)
            sage: vec * g
            2*(2,3)*(1,3) + 2*(2,3)*(1,2) - 3*()*(2,3)
            sage: supp = vec.leading_support(); supp
            11
            sage: E2._action_on_basis(g, supp, True)
            (2,3)*(1,3)
            sage: E2._action_on_basis(g, supp, False)
            (2,3)*(1,3)
        """
        B = self._rep.basis()
        if vec_on_left:
            temp = self._extalg.prod(self._from_repr_to_ext(B[self._basis_order[bk]] * g)
                                     for bk in b)
        else:
            temp = self._extalg.prod(self._from_repr_to_ext(g * B[self._basis_order[bk]])
                                     for bk in b)
        return self.element_class(self, temp._monomial_coefficients)


class Representation_ExteriorAlgebra(Representation_Exterior):
    r"""
    The exterior algebra representation.
    """
    def __init__(self, rep, degree=None, category=None, **options):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: G = groups.matrix.GL(3, 2)
            sage: R = G.regular_representation(side='right')
            sage: E0 = R.exterior_power(0)
            sage: E0.category()
            Category of finite dimensional algebras with basis over Integer Ring
            sage: TestSuite(E0).run()

            sage: G = groups.matrix.GL(2, 3)
            sage: L = G.regular_representation(side='left')
            sage: E = L.exterior_power()
            sage: E.category()
            Category of finite dimensional algebras with basis over Integer Ring
            sage: TestSuite(E).run()
        """
        R = rep.base_ring()
        from sage.categories.algebras_with_basis import AlgebrasWithBasis
        category = AlgebrasWithBasis(R).or_subcategory(category)
        Representation_Exterior.__init__(self, rep, degree=degree, category=category, **options)

    @cached_method
    def one_basis(self):
        r"""
        Return the basis element indexing `1` in ``self`` if it exists.

        EXAMPLES::

            sage: S3 = SymmetricGroup(3)
            sage: L = S3.regular_representation(side='left')
            sage: E = L.exterior_power()
            sage: E.one_basis()
            0
            sage: E0 = L.exterior_power(0)
            sage: E0.one_basis()
            0
        """
        return self._indices([])

    def product_on_basis(self, x, y):
        r"""
        Return the product of basis elements indexed by ``x`` and ``y``.

        EXAMPLES::

            sage: S3 = SymmetricGroup(3)
            sage: L = S3.regular_representation(side='left')
            sage: E = L.exterior_power()
            sage: B = list(E.basis())
            sage: B[:7]
            [1, (), (1,3,2), (1,2,3), (2,3), (1,3), (1,2)]
            sage: B[2] * B[4]  # indirect doctest
            (1,3,2)*(2,3)
        """
        B = self._extalg.basis()
        temp = B[x] * B[y]
        return self.element_class(self, temp._monomial_coefficients)


class Representation_Symmetric(Representation_abstract, CombinatorialFreeModule):
    r"""
    The symmetric power representation in a fixed degree.
    """
    def __init__(self, rep, degree, **options):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: G = groups.matrix.GL(3, 2)
            sage: R = G.regular_representation(side='right')
            sage: S2 = R.symmetric_power(2)
            sage: TestSuite(S2).run()
            sage: S0 = R.symmetric_power(0)
            sage: TestSuite(S2).run()

            sage: R.symmetric_power(-2)
            Traceback (most recent call last):
            ...
            ValueError: the degree must be a nonnegative integer
            sage: R.symmetric_power(3/2)
            Traceback (most recent call last):
            ...
            ValueError: the degree must be a nonnegative integer
        """
        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
        from sage.combinat.integer_vector import IntegerVectors
        from sage.rings.integer_ring import ZZ
        self._degree = degree
        self._rep = rep
        R = rep.base_ring()
        dim = rep.dimension()
        if degree not in ZZ or degree < 0:
            raise ValueError(f"the degree must be a nonnegative integer")
        self._symalg = PolynomialRing(R, 'e', dim)
        self._basis_order = list(rep.basis().keys())
        G = self._symalg.gens()
        self._inv_map = {b: G[i] for i, b in enumerate(self._basis_order)}
        ind = IntegerVectors(degree, dim)
        CombinatorialFreeModule.__init__(self, rep.base_ring(), ind, **options)
        Representation_abstract.__init__(self, rep.semigroup(), rep.side(), rep.semigroup_algebra())

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: DC3 = groups.permutation.DiCyclic(3)
            sage: L = DC3.regular_representation(QQ, side='left')
            sage: L.symmetric_power(7)
            Symmetric power representation of Left Regular Representation of
             Dicyclic group of order 12 as a permutation group over Rational Field
             in degree 7
        """
        return "Symmetric power representation of {} in degree {}".format(repr(self._rep), self._degree)

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: DC3 = groups.permutation.DiCyclic(3)
            sage: L = DC3.regular_representation(QQ, side='left')
            sage: latex(L.symmetric_power(4))
            S^{4} ...
        """
        from sage.misc.latex import latex
        return "S^{{{}}} {}".format(self._degree, latex(self._rep))

    def _repr_term(self, m):
        r"""
        Return a string representation of the basis element indexed by
        ``m``.

        EXAMPLES::

            sage: DC3 = groups.permutation.DiCyclic(3)
            sage: L = DC3.regular_representation(QQ, side='left')
            sage: S2L = L.symmetric_power(2)
            sage: S2L.an_element()
            3*()*(5,7,6) + 2*()*(5,6,7) + 2*()^2
            sage: S2L._repr_term(S2L.an_element().trailing_support())
            '()*(5,7,6)'
            sage: S2L._repr_term(S2L.an_element().leading_support())
            '()^2'
            sage: L.symmetric_power(0).an_element()
            2
        """
        if not self._degree:
            return '1'
        B = self._rep.basis()
        return '*'.join(repr(B[self._basis_order[i]]) if e == 1 else repr(B[self._basis_order[i]]) + f'^{e}'
                        for i,e in enumerate(m) if e)

    def _ascii_art_term(self, m):
        r"""
        Return ascii art for the basis element indexed by ``m``.

        EXAMPLES::

            sage: DC3 = groups.permutation.DiCyclic(3)
            sage: L = DC3.regular_representation(QQ, side='left')
            sage: S2L = L.symmetric_power(2)
            sage: S2L._ascii_art_term(S2L.an_element().leading_support())
              2
            ()
            sage: ascii_art(S2L.an_element())
                                              2
            3*()*(5,7,6) + 2*()*(5,6,7) + 2*()
            sage: ascii_art(L.symmetric_power(0).an_element())
            2*1
        """
        from sage.typeset.ascii_art import ascii_art
        if not self._degree:
            return ascii_art('1')
        B = self._rep.basis()
        ret = ascii_art("")
        for i, e in enumerate(m):
            if not e:
                continue
            cur = ascii_art(B[self._basis_order[i]])
            if e > 1:
                cur += ascii_art(e, baseline=-cur.height())
            if ret:
                ret += ascii_art('*')
            ret += cur
        return ret

    def _unicode_art_term(self, m):
        r"""
        Return unicode art for the basis element indexed by ``m``.

        EXAMPLES::

            sage: DC3 = groups.permutation.DiCyclic(3)
            sage: L = DC3.regular_representation(QQ, side='left')
            sage: S2L = L.symmetric_power(2)
            sage: S2L._unicode_art_term(S2L.an_element().leading_support())
              2
            ()
            sage: unicode_art(S2L.an_element())
                                              2
            3*()*(5,7,6) + 2*()*(5,6,7) + 2*()
            sage: unicode_art(L.symmetric_power(0).an_element())
            2*1
        """
        from sage.typeset.unicode_art import unicode_art
        if not self._degree:
            return unicode_art('1')
        B = self._rep.basis()
        ret = unicode_art("")
        for i, e in enumerate(m):
            if not e:
                continue
            cur = unicode_art(B[self._basis_order[i]])
            if e > 1:
                cur += unicode_art(e, baseline=-cur.height())
            if ret:
                ret += unicode_art('*')
            ret += cur
        return ret

    def _latex_term(self, m):
        r"""
        Return a `\LaTeX` representation of the basis element indexed
        by ``m``.

        EXAMPLES::

            sage: DC3 = groups.permutation.DiCyclic(3)
            sage: L = DC3.regular_representation(QQ, side='left')
            sage: S2L = L.symmetric_power(2)
            sage: S2L._latex_term(S2L.an_element().leading_support())
            '1 ^{2}'
            sage: latex(S2L.an_element())
            3 1 (5,7,6) + 2 1 (5,6,7) + 2 1 ^{2}
            sage: latex(L.symmetric_power(0).an_element())
            2
        """
        if not self._degree:
            return '1'
        from sage.misc.latex import latex
        B = self._rep.basis()
        return " ".join(latex(B[self._basis_order[i]]) if e == 1 else latex(B[self._basis_order[i]]) + f"^{{{e}}}"
                        for i, e in enumerate(m) if e)

    def _from_repr_to_sym(self, elt):
        r"""
        Return the element ``elt`` from the defining representation
        to the corresponding exterior algebra.

        EXAMPLES::

            sage: G = groups.matrix.GL(2, 2)
            sage: L = G.regular_representation(side='left')
            sage: S3L = L.symmetric_power(3)
            sage: S3L._from_repr_to_sym(sum(i*b for i,b in enumerate(L.basis(), start=1)))
            e0 + 2*e1 + 3*e2 + 4*e3 + 5*e4 + 6*e5
        """
        return self._symalg.sum(c * self._inv_map[k]
                                for k, c in elt._monomial_coefficients.items())

    def _semigroup_action(self, g, vec, vec_on_left):
        r"""
        Return the action of the semigroup element ``g`` on the
        vector ``vec`` of ``self``.

        EXAMPLES::

            sage: DC3 = groups.permutation.DiCyclic(3)
            sage: g = DC3.an_element(); g
            (1,4,2,3)(5,6)
            sage: R = DC3.regular_representation(side='right')
            sage: S2L = R.symmetric_power(2)
            sage: vec = S2L.an_element(); vec
            3*()*(5,7,6) + 2*()*(5,6,7) + 2*()^2
            sage: S2L._semigroup_action(g, vec, True)
            3*(1,4,2,3)(5,6)*(1,4,2,3)(5,7) + 2*(1,4,2,3)(5,6)^2
             + 2*(1,4,2,3)(6,7)*(1,4,2,3)(5,6)
            sage: S2L._semigroup_action(g, vec, False)
            3*(1,3,2,4)(5,6)*(1,3,2,4)(5,7) + 2*(1,3,2,4)(5,6)^2
             + 2*(1,3,2,4)(6,7)*(1,3,2,4)(5,6)
        """
        return self.linear_combination(((self._action_on_basis(g, b, vec_on_left), c)
                                        for b, c in vec._monomial_coefficients.items()), not vec_on_left)

    def _action_on_basis(self, g, b, vec_on_left):
        r"""
        Return the action of ``g`` on the basis element indexed by ``b``.

        EXAMPLES::

            sage: S3 = SymmetricGroup(3)
            sage: g = S3.an_element(); g
            (2,3)
            sage: L = S3.regular_representation(side='left')
            sage: S2L = L.symmetric_power(2)
            sage: vec = S2L.an_element(); vec
            3*()*(1,2,3) + 2*()*(1,3,2) + 2*()^2
            sage: g * vec
            3*(2,3)*(1,2) + 2*(2,3)*(1,3) + 2*(2,3)^2
            sage: vec * g
            3*(2,3)*(1,2) + 2*(2,3)*(1,3) + 2*(2,3)^2
            sage: supp = vec.leading_support(); supp
            [2, 0, 0, 0, 0, 0]
            sage: S2L._action_on_basis(g, supp, True)
            (2,3)^2
            sage: S2L._action_on_basis(g, supp, False)
            (2,3)^2
        """
        B = self._rep.basis()
        if vec_on_left:
            temp = self._symalg.prod(self._from_repr_to_sym(B[self._basis_order[bk]] * g) ** e
                                     for bk, e in enumerate(b))
        else:
            temp = self._symalg.prod(self._from_repr_to_sym(g * B[self._basis_order[bk]]) ** e
                                     for bk, e in enumerate(b))
        ind = self._indices
        data = {ind(mon.exponents()[0]): c for c, mon in temp}
        return self.element_class(self, data)


class RegularRepresentation(Representation):
    r"""
    The regular representation of a semigroup.

    The left regular representation of a semigroup `S` over a commutative
    ring `R` is the semigroup ring `R[S]` equipped with the left
    `S`-action `x b_y = b_{xy}`, where `(b_z)_{z \in S}` is the natural
    basis of `R[S]` and `x,y \in S`.

    INPUT:

    - ``semigroup`` -- a semigroup
    - ``base_ring`` -- the base ring for the representation
    - ``side`` -- (default: ``'left'``) whether this is a
      ``'left'`` or ``'right'`` representation

    REFERENCES:

    - :wikipedia:`Regular_representation`
    """
    def __init__(self, semigroup, base_ring, side='left'):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: G = groups.permutation.Dihedral(4)
            sage: R = G.regular_representation()
            sage: TestSuite(R).run()
        """
        if side == "left":
            on_basis = self._left_on_basis
        else:
            on_basis = self._right_on_basis
        module = semigroup.algebra(base_ring)
        Representation.__init__(self, semigroup, module, on_basis, side)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: G = groups.permutation.Dihedral(4)
            sage: G.regular_representation()
            Left Regular Representation of Dihedral group of order 8
             as a permutation group over Integer Ring
            sage: G.regular_representation(side='right')
            Right Regular Representation of Dihedral group of order 8
             as a permutation group over Integer Ring
        """
        if self._left_repr:
            base = "Left Regular Representation"
        else:
            base = "Right Regular Representation"
        return base + " of {} over {}".format(self._semigroup, self.base_ring())

    def _left_on_basis(self, g, m):
        """
        Return the left action of ``g`` on ``m``.

        EXAMPLES::

            sage: G = groups.permutation.Dihedral(4)
            sage: R = G.regular_representation()
            sage: R._test_representation()  # indirect doctest
        """
        return self.monomial(g * m)

    def _right_on_basis(self, g, m):
        """
        Return the right action of ``g`` on ``m``.

        EXAMPLES::

            sage: G = groups.permutation.Dihedral(4)
            sage: R = G.regular_representation(side='right')
            sage: R._test_representation()  # indirect doctest
        """
        return self.monomial(m * g)


class TrivialRepresentation(Representation_abstract, CombinatorialFreeModule):
    """
    The trivial representation of a semigroup.

    The trivial representation of a semigroup `S` over a commutative ring
    `R` is the `1`-dimensional `R`-module on which every element of `S`
    acts by the identity.

    This is simultaneously a left and right representation.

    INPUT:

    - ``semigroup`` -- a semigroup
    - ``base_ring`` -- the base ring for the representation

    REFERENCES:

    - :wikipedia:`Trivial_representation`
    """
    def __init__(self, semigroup, base_ring):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: G = groups.permutation.PGL(2, 3)
            sage: V = G.trivial_representation()
            sage: TestSuite(V).run()
        """
        cat = Modules(base_ring).WithBasis().FiniteDimensional()
        from sage.sets.finite_enumerated_set import FiniteEnumeratedSet
        indices = FiniteEnumeratedSet(['v'])
        CombinatorialFreeModule.__init__(self, base_ring, indices, category=cat)
        Representation_abstract.__init__(self, semigroup, "twosided")

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: G = groups.permutation.Dihedral(4)
            sage: G.trivial_representation()
            Trivial representation of Dihedral group of order 8
             as a permutation group over Integer Ring
        """
        return "Trivial representation of {} over {}".format(self._semigroup,
                                                             self.base_ring())

    def _semigroup_action(self, g, vec, vec_on_left):
        r"""
        Return the action of the semigroup element ``g`` on the
        vector ``vec`` of ``self``.

        EXAMPLES::

            sage: SGA = SymmetricGroupAlgebra(QQ, 4)
            sage: V = SGA.trivial_representation()
            sage: x = V.an_element()
            sage: V._semigroup_action(SGA.group().random_element(), x, True) == x
            True
        """
        return vec

    class Element(Representation_abstract.Element):
        def _acted_upon_(self, scalar, self_on_left=False):
            """
            Return the action of ``scalar`` on ``self``.

            EXAMPLES::

                sage: SGA = SymmetricGroupAlgebra(QQ, 3)
                sage: V = SGA.trivial_representation()
                sage: x = V.an_element()
                sage: 2 * x
                4*B['v']
                sage: all(x * b == x for b in SGA.basis())
                True
                sage: all(b * x == x for b in SGA.basis())
                True
                sage: z = V.zero()
                sage: all(b * z == z for b in SGA.basis())
                True

                sage: H = groups.permutation.Dihedral(5)
                sage: G = SymmetricGroup(5)
                sage: G.has_coerce_map_from(H)
                True
                sage: R = G.trivial_representation(QQ)
                sage: H.an_element() * R.an_element()
                2*B['v']

                sage: AG = G.algebra(QQ)
                sage: AG.an_element() * R.an_element()
                14*B['v']

                sage: AH = H.algebra(ZZ)
                sage: AG.has_coerce_map_from(AH)
                True
                sage: AH.an_element() * R.an_element()
                14*B['v']
            """
            if isinstance(scalar, Element):
                P = self.parent()
                if P._semigroup.has_coerce_map_from(scalar.parent()):
                    return self
                if P._semigroup_algebra.has_coerce_map_from(scalar.parent()):
                    if not self:
                        return self
                    scalar = P._semigroup_algebra(scalar)
                    d = self.monomial_coefficients(copy=True)
                    d['v'] *= sum(scalar.coefficients())
                    return P._from_dict(d)
            return CombinatorialFreeModule.Element._acted_upon_(self, scalar, self_on_left)


class SignRepresentation_abstract(Representation_abstract, CombinatorialFreeModule):
    """
    Generic implementation of a sign representation.

    The sign representation of a semigroup `S` over a commutative ring
    `R` is the `1`-dimensional `R`-module on which every element of `S`
    acts by `1` if order of element is even (including 0) or `-1` if
    order of element if odd.

    This is simultaneously a left and right representation.

    INPUT:

    - ``permgroup`` -- a permgroup
    - ``base_ring`` -- the base ring for the representation
    - ``sign_function`` -- a function which returns `1` or `-1` depending
      on the elements sign

    REFERENCES:

    - :wikipedia:`Representation_theory_of_the_symmetric_group`
    """
    def __init__(self, group, base_ring, sign_function=None):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: G = groups.permutation.PGL(2, 3)
            sage: V = G.sign_representation()
            sage: TestSuite(V).run()
        """
        self.sign_function = sign_function
        if sign_function is None:
            try:
                self.sign_function = self._default_sign
            except AttributeError:
                raise TypeError("a sign function must be given")

        cat = Modules(base_ring).WithBasis().FiniteDimensional()

        CombinatorialFreeModule.__init__(self, base_ring, ["v"], category=cat)
        Representation_abstract.__init__(self, group, "twosided")

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: G = groups.permutation.Dihedral(4)
            sage: G.sign_representation()
            Sign representation of Dihedral group of order 8
             as a permutation group over Integer Ring
        """
        return "Sign representation of {} over {}".format(
            self._semigroup, self.base_ring()
        )

    def _semigroup_action(self, g, vec, vec_on_left):
        r"""
        Return the action of the semigroup element ``g`` on the
        vector ``vec`` of ``self``.

        EXAMPLES::

            sage: G = PermutationGroup(gens=[(1,2,3), (1,2)])
            sage: S = G.sign_representation()
            sage: x = S.an_element(); x
            2*B['v']
            sage: s, c = G.gens(); c
            (1,2,3)
            sage: S._semigroup_action(s, x, True)
            -2*B['v']
            sage: S._semigroup_action(s, x, False)
            -2*B['v']
            sage: s * x
            -2*B['v']
            sage: s*x*s
            2*B['v']
            sage: s*x*s*s*c
            -2*B['v']
            sage: A = G.algebra(ZZ)
            sage: s,c = A.algebra_generators()
            sage: c
            (1,2,3)
            sage: s
            (1,2)
            sage: c*x
            2*B['v']
            sage: c*c*x
            2*B['v']
            sage: c*x*s
            -2*B['v']
            sage: c*x*s*s
            2*B['v']
            sage: (c+s)*x
            0
            sage: (c-s)*x
            4*B['v']

            sage: H = groups.permutation.Dihedral(4)
            sage: G = SymmetricGroup(4)
            sage: G.has_coerce_map_from(H)
            True
            sage: R = G.sign_representation()
            sage: H.an_element() * R.an_element()
            -2*B['v']

            sage: AG = G.algebra(ZZ)
            sage: AH = H.algebra(ZZ)
            sage: AG.has_coerce_map_from(AH)
            True
            sage: AH.an_element() * R.an_element()
            -2*B['v']
        """
        return vec if self.sign_function(g) > 0 else -vec


class SignRepresentationPermgroup(SignRepresentation_abstract):
    """
    The sign representation for a permutation group.

    EXAMPLES::

        sage: G = groups.permutation.PGL(2, 3)
        sage: V = G.sign_representation()
        sage: TestSuite(V).run()
    """
    def _default_sign(self, elem):
        """
        Return the sign of the element.

        INPUT:

        - ``elem`` -- the element of the group

        EXAMPLES::

            sage: G = groups.permutation.PGL(2, 3)
            sage: V = G.sign_representation()
            sage: elem = G.an_element()
            sage: elem
            (1,2,4,3)
            sage: V._default_sign(elem)
            -1
        """
        return elem.sign()


class SignRepresentationMatrixGroup(SignRepresentation_abstract):
    """
    The sign representation for a matrix group.

    EXAMPLES::

        sage: G = groups.permutation.PGL(2, 3)
        sage: V = G.sign_representation()
        sage: TestSuite(V).run()
    """
    def _default_sign(self, elem):
        """
        Return the sign of the element.

        INPUT:

        - ``elem`` -- the element of the group

        EXAMPLES::

            sage: G = GL(2, QQ)
            sage: V = G.sign_representation()
            sage: m = G.an_element()
            sage: m
            [1 0]
            [0 1]
            sage: V._default_sign(m)
            1
        """
        return 1 if elem.matrix().det() > 0 else -1


class SignRepresentationCoxeterGroup(SignRepresentation_abstract):
    r"""
    The sign representation for a Coxeter group.

    EXAMPLES::

        sage: G = WeylGroup(["A", 1, 1])
        sage: V = G.sign_representation()
        sage: TestSuite(V).run()

        sage: # optional - gap3
        sage: W = CoxeterGroup(['B', 3], implementation="coxeter3")
        sage: S = W.sign_representation()
        sage: TestSuite(S).run()
    """
    def _default_sign(self, elem):
        """
        Return the sign of the element.

        INPUT:

        - ``elem`` -- the element of the group

        EXAMPLES::

            sage: G = WeylGroup(["A", 1, 1])
            sage: elem = G.an_element()
            sage: V = G.sign_representation()
            sage: V._default_sign(elem)
            1

            sage: # optional - gap3
            sage: W = CoxeterGroup(['B', 3], implementation="coxeter3")
            sage: S = W.sign_representation()
            sage: elem = W.an_element()
            sage: S._default_sign(elem)
            1
        """
        return -1 if elem.length() % 2 else 1


class ReflectionRepresentation(Representation_abstract, CombinatorialFreeModule):
    r"""
    The reflection representation of a Coxeter group.

    This is the canonical faithful representation of a Coxeter group.

    EXAMPLES::

        sage: W = CoxeterGroup(['B', 4])
        sage: R = W.reflection_representation()
        sage: all(g.matrix() == R.representation_matrix(g) for g in W)
        True
    """
    @staticmethod
    def __classcall_private__(cls, W, base_ring=None):
        r"""
        Normalize input to ensure a unique representation.

        EXAMPLES::

            sage: W = CoxeterGroup(['D', 4])
            sage: R1 = W.reflection_representation()
            sage: R2 = W.reflection_representation(ZZ)
            sage: R1 is R2
            True
        """
        if base_ring is None:
            base_ring = W.one().canonical_matrix().base_ring()
        return super().__classcall__(cls, W, base_ring)

    def __init__(self, W, base_ring):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: W = CoxeterGroup(['C', 3])
            sage: R = W.reflection_representation()
            sage: TestSuite(R).run()

            sage: W = CoxeterGroup(['E', 6, 1])
            sage: R = W.reflection_representation()
            sage: TestSuite(R).run()
        """
        self._W = W
        rk = W.coxeter_matrix().rank()
        from sage.sets.finite_enumerated_set import FiniteEnumeratedSet
        indices = FiniteEnumeratedSet(range(rk))
        CombinatorialFreeModule.__init__(self, base_ring, indices, prefix='e', bracket=False)
        Representation_abstract.__init__(self, W, "left")

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: W = CoxeterGroup(['E', 8])
            sage: W.reflection_representation()
            Reflection representation of Finite Coxeter group over Integer Ring with Coxeter matrix:
            [1 2 3 2 2 2 2 2]
            [2 1 2 3 2 2 2 2]
            [3 2 1 3 2 2 2 2]
            [2 3 3 1 3 2 2 2]
            [2 2 2 3 1 3 2 2]
            [2 2 2 2 3 1 3 2]
            [2 2 2 2 2 3 1 3]
            [2 2 2 2 2 2 3 1]
        """
        return "Reflection representation of {}".format(self._W)

    def _semigroup_action(self, g, vec, vec_on_left):
        r"""
        Return the action of the Coxeter group element ``g`` on the
        vector ``vec`` of ``self``.

        EXAMPLES::

            sage: W = CoxeterGroup(['E', 8])
            sage: R = W.reflection_representation()
            sage: g = W.an_element()
            sage: g == ~g
            False
            sage: vec = R.an_element()
            sage: R._semigroup_action(g, vec, True)
            e0 - 2*e1 - 2*e2 - 4*e3 - 4*e4 - 4*e5 - 4*e6 - 4*e7
            sage: R._semigroup_action(g, vec, False)
            2*e0 + 3*e1 + 4*e2 + 5*e3
        """
        if vec_on_left:
            g = ~g
        return self.from_vector(g.canonical_matrix() * vec.to_vector())


class NaturalMatrixRepresentation(Representation):
    r"""
    The natural representation of a matrix semigroup.

    A matrix semigroup is defined by its representation on a (finite
    dimensional) vector space `V`, which is called the *natural
    representation*.

    INPUT:

    - ``matrix_semigroup`` -- a matrix semigroup
    - ``base_ring`` -- (optional) the base ring; the default is the base ring
      of the semigroup
    """
    @staticmethod
    def __classcall_private__(cls, semigroup, base_ring=None):
        r"""
        Normalize input to ensure a unique representation.

        EXAMPLES::

            sage: G = groups.matrix.SO(5, 2)
            sage: N1 = G.natural_representation()
            sage: N2 = G.natural_representation(GF(2))
            sage: N1 is N2
            True
        """
        if base_ring is None:
            base_ring = semigroup.base_ring()
        return super().__classcall__(cls, semigroup, base_ring)

    def __init__(self, semigroup, base_ring):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: G = groups.matrix.SU(2, 2)
            sage: N = G.natural_representation()
            sage: TestSuite(N).run()

            sage: G = groups.matrix.Sp(4, 2)
            sage: N = G.natural_representation()
            sage: TestSuite(N).run()  # long time
        """
        base_ring = semigroup.base_ring()
        n = semigroup.degree()
        module = CombinatorialFreeModule(base_ring, list(range(n)), prefix='e')
        Representation.__init__(self, semigroup, module, None, "left")

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: G = groups.matrix.GL(3, 2)
            sage: G.natural_representation()
            Natural representation of General Linear Group of degree 3
             over Finite Field of size 2
        """
        return "Natural representation of {}".format(self._semigroup)

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: G = groups.matrix.GL(3, 2)
            sage: latex(G.natural_representation())
            \Bold{F}_{2} ^{3}
        """
        from sage.misc.latex import latex
        return latex(self.base_ring()) + "^{{{}}}".format(self.dimension())

    def _semigroup_action(self, g, vec, vec_on_left):
        r"""
        Return the action of the semigroup element ``g`` on the
        vector ``vec`` of ``self``.

        EXAMPLES::

            sage: G = groups.matrix.SL(2, 3)
            sage: N = G.natural_representation()
            sage: N._test_representation()  # indirect doctest
        """
        if self._left_repr == vec_on_left:
            g = ~g
        # TODO: currently, the slowest part is getting the matrix from g
        return self.from_vector(g.matrix() * vec.to_vector())


class SchurFunctorRepresentation(Subrepresentation):
    r"""
    The representation constructed by the Schur functor.

    Let `G` be a semigroup and let `V` be a representation of `G`. The
    *Schur functor* for a partition `\lambda` of size `k` is the functor
    `\mathbb{S}_{\lambda}` that sends `V` to the `G`-subrepresentation of
    `V^{\otimes k}` spanned by `(v_1 \otimes \cdots \otimes v_k) c_{\lambda}`,
    where `c_{\lambda}` is the :meth:`Young symmetrizer
    <sage.combinat.symmetric_group_algebra.SymmetricGroupAlgebra.young_symmetrizer>`
    corresponding to `\lambda`. When `G = GL_n(F)`, the Schur functor image
    `\mathbb{S}_{\lambda} F^n` is the (irreducible when `F` has characteristic
    `0`) highest representation of shape `\lambda`.

    EXAMPLES::

        sage: G = groups.permutation.Dihedral(3)
        sage: V = G.regular_representation()
        sage: S21V = V.schur_functor([2, 1])
        sage: S21V.dimension()
        70
        sage: SemistandardTableaux([2,1], max_entry=V.dimension()).cardinality()
        70
        sage: chi = G.character([S21V.representation_matrix(g).trace()
        ....:                    for g in G.conjugacy_classes_representatives()])
        sage: chi.values()
        [70, 0, -2]
        sage: G.character_table()
        [ 1 -1  1]
        [ 2  0 -1]
        [ 1  1  1]
        sage: for m, phi in chi.decompose():
        ....:     print(m, phi.values())
        11 [1, -1, 1]
        11 [1, 1, 1]
        24 [2, 0, -1]

        sage: # long time
        sage: S211V = V.schur_functor([2, 1, 1])
        sage: S211V.dimension()
        105
        sage: SemistandardTableaux([2, 1, 1], max_entry=V.dimension()).cardinality()
        105
        sage: chi = G.character([S211V.representation_matrix(g).trace()
        ....:                    for g in G.conjugacy_classes_representatives()])
        sage: chi.values()
        [105, -3, 0]
        sage: for m, phi in chi.decompose():
        ....:     print(m, phi.values())
        19 [1, -1, 1]
        16 [1, 1, 1]
        35 [2, 0, -1]

    An example with the cyclic group::

        sage: C3 = groups.permutation.Cyclic(3)
        sage: V3 = C3.regular_representation()
        sage: S31V3 = V3.schur_functor([3, 1])
        sage: S31V3.dimension()
        15
        sage: chi = C3.character([S31V3.representation_matrix(g).trace()
        ....:                     for g in C3.conjugacy_classes_representatives()])
        sage: chi.values()
        [15, 0, 0]
        sage: C3.character_table()
        [         1          1          1]
        [         1      zeta3 -zeta3 - 1]
        [         1 -zeta3 - 1      zeta3]
        sage: for m, phi in chi.decompose():
        ....:     print(m, phi.values())
        5 [1, 1, 1]
        5 [1, -zeta3 - 1, zeta3]
        5 [1, zeta3, -zeta3 - 1]

    An example of `GL_3(\GF{2})`::

        sage: G = groups.matrix.GL(3, 2)
        sage: from sage.modules.with_basis.representation import Representation
        sage: N = G.natural_representation()
        sage: S21N = N.schur_functor([2, 1])
        sage: S21N.dimension()
        8

    An example with the Weyl/Coxeter group of type `C_3`::

        sage: G = WeylGroup(['C', 3], prefix='s')
        sage: R = G.reflection_representation()
        sage: S = R.schur_functor([3, 2, 1])
        sage: g = G.an_element(); g
        s1*s2*s3
        sage: v = S.an_element(); v
        2*S[0] + 2*S[1] + 3*S[2]
        sage: v * g
        -(2*a+1)*S[0] - (a+2)*S[1] - (2*a-2)*S[2] + (2*a-2)*S[3]
         - (-2*a+4)*S[4] + (-2*a+4)*S[5] + 2*S[6] + 2*a*S[7]
        sage: g * v
        3*S[0] + (-2*a+5)*S[2] + 3*a*S[4] - (5*a-2)*S[6] - 6*S[7]
    """
    @staticmethod
    def __classcall_private__(cls, V, shape):
        r"""
        Normalize input to ensure a unique representation.

        EXAMPLES::

            sage: from sage.modules.with_basis.representation import SchurFunctorRepresentation
            sage: G = groups.permutation.Cyclic(3)
            sage: V = G.regular_representation()
            sage: S1 = SchurFunctorRepresentation(V, Partition([2, 1, 1]))
            sage: S2 = SchurFunctorRepresentation(V, (2, 1, 1))
            sage: S1 is S2
            True
        """
        from sage.combinat.partition import _Partitions
        return super().__classcall__(cls, V, _Partitions(shape))

    def __init__(self, V, shape):
        r"""
        Initialize ``self``.

        INPUT:

        - ``V`` -- a representation
        - ``shape`` -- a partition

        EXAMPLES::

            sage: G = groups.permutation.Dihedral(3)
            sage: V = G.regular_representation()
            sage: S = V.schur_functor([2, 1])
            sage: TestSuite(S).run()

            sage: G = CoxeterGroup(['B', 3])
            sage: R = G.reflection_representation()
            sage: S = R.schur_functor([3, 2, 1])
            sage: TestSuite(S).run()

            sage: G = groups.matrix.GL(3, 2)
            sage: N = G.natural_representation()
            sage: S = N.schur_functor([2, 1])
            sage: TestSuite(S).run()
        """
        from sage.combinat.symmetric_group_algebra import SymmetricGroupAlgebra
        from sage.groups.perm_gps.permgroup_named import SymmetricGroup
        from sage.categories.tensor import tensor
        from sage.matrix.matrix_space import MatrixSpace

        R = V.base_ring()
        self._shape = shape
        self._module = V
        d = sum(shape)

        if not isinstance(V, CombinatorialFreeModule):
            keys = list(range(V.dimension()))
            V = CombinatorialFreeModule(R, keys)
        else:
            keys = list(V.basis().keys())

        ambient = tensor([V]*d)
        I = ambient.indices()
        cla = SymmetricGroupAlgebra(R, SymmetricGroup(d)).young_symmetrizer(shape)
        mc = cla.monomial_coefficients(copy=False)
        gens = [ambient.sum_of_terms((tuple([k[i-1] for i in p.tuple()]), coeff)
                                     for p, coeff in mc.items())
                for k in ambient.basis().keys()]
        support_order = ambient._compute_support_order(gens, None)
        from sage.sets.family import Family
        gens = Family(ambient.echelon_form(gens, order=support_order))
        cat = Modules(ambient.category().base_ring()).WithBasis().Subobjects()
        Subrepresentation.__init__(self, gens, support_order, ambient, unitriangular=False, category=cat, prefix='S')

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: G = groups.matrix.GL(4, 2)
            sage: N = G.natural_representation()
            sage: N.schur_functor([2, 1])
            Schur functor of shape [2, 1] image of Natural representation of
             General Linear Group of degree 4 over Finite Field of size 2
        """
        return "Schur functor of shape {} image of {}".format(self._shape, self._module)

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: G = groups.matrix.GL(4, 2)
            sage: N = G.natural_representation()
            sage: latex(N.schur_functor([2, 1]))
            \mathbb{S}_{{\def\lr#1{\multicolumn{1}{|@{\hspace{.6ex}}c@{\hspace{.6ex}}|}{\raisebox{-.3ex}{$#1$}}}
            \raisebox{-.6ex}{$\begin{array}[b]{*{2}c}\cline{1-2}
            \lr{\phantom{x}}&\lr{\phantom{x}}\\\cline{1-2}
            \lr{\phantom{x}}\\\cline{1-1}
            \end{array}$}
            }}(\Bold{F}_{2} ^{4})
        """
        from sage.misc.latex import latex
        return "\\mathbb{{S}}_{{{}}}({})".format(latex(self._shape), latex(self._module))

    Element = Subrepresentation.Element
