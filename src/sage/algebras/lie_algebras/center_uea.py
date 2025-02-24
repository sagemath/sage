"""
Center of a Universal Enveloping Algebra

AUTHORS:

- Travis Scrimshaw (2024-01-02): Initial version
"""

# ****************************************************************************
#       Copyright (C) 2024 Travis Scrimshaw <tcscrims at gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

# from sage.structure.unique_representation import UniqueRepresentation
# from sage.structure.parent import Parent
from sage.combinat.free_module import CombinatorialFreeModule
from sage.combinat.integer_lists.invlex import IntegerListsLex
from sage.matrix.constructor import matrix
from sage.monoids.indexed_free_monoid import IndexedFreeAbelianMonoid
from sage.monoids.indexed_free_monoid import IndexedMonoid
from sage.combinat.root_system.coxeter_group import CoxeterGroup
from sage.combinat.integer_vector_weighted import iterator_fast as intvecwt_iterator
from sage.sets.non_negative_integers import NonNegativeIntegers
from sage.sets.finite_enumerated_set import FiniteEnumeratedSet
from sage.sets.family import Family
from sage.rings.integer_ring import ZZ
from sage.categories.kac_moody_algebras import KacMoodyAlgebras
from sage.categories.finite_dimensional_lie_algebras_with_basis import FiniteDimensionalLieAlgebrasWithBasis
from sage.categories.graded_algebras_with_basis import GradedAlgebrasWithBasis
from sage.categories.fields import Fields
from sage.categories.monoids import Monoids
from sage.categories.enumerated_sets import EnumeratedSets
from sage.misc.cachefunc import cached_method
from sage.misc.lazy_attribute import lazy_attribute
from sage.data_structures.blas_dict import iaxpy
from collections import deque


class CenterIndices(IndexedFreeAbelianMonoid):
    r"""
    Set of basis indices for the center of a universal enveloping algebra.

    This also constructs the lift from the center to the universal enveloping
    algebra as part of computing the generators and basis elements. The
    basic algorithm is to construct the centralizer of each filtered
    component in increasing order (as each is a finite dimensional vector
    space). For more precise details, see [Motsak2006]_.
    """
    @staticmethod
    def __classcall__(cls, center):
        r"""
        Normalize input to ensure a unique representation.

        EXAMPLES::

            sage: from sage.algebras.lie_algebras.center_uea import CenterIndices
            sage: g = lie_algebras.pwitt(GF(3), 3)
            sage: U = g.pbw_basis()
            sage: Z = U.center()
            sage: CenterIndices(Z) is CenterIndices(Z)
            True
        """
        return super(IndexedMonoid, cls).__classcall__(cls, center)

    def __init__(self, center, indices=None):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: g = lie_algebras.pwitt(GF(5), 5)
            sage: U = g.pbw_basis()
            sage: Z = U.center()
            sage: I = Z.indices()
            sage: TestSuite(I).run(max_runs=7)
        """
        if indices is None:
            indices = NonNegativeIntegers()
        category = Monoids() & EnumeratedSets().Infinite()
        IndexedFreeAbelianMonoid.__init__(self, indices, prefix='Z', category=category)

        self._center = center
        self._envelop_alg = self._center._envelop_alg
        self._g = self._envelop_alg._g
        # The _lift_map will be a dict with the keys being the degree and the values
        #   given as dicts with the keys being the leading support. This will be
        #   used to do the corresponding reductions.
        self._lift_map = {0: {self._envelop_alg.one_basis(): self._envelop_alg.one()}}
        self._cur_deg = 0
        self._cur_vecs = deque()  # we do a lot of deletions in the middle
        # The _cur_basis is a mapping from the leading supports to a monoid element of self
        self._cur_basis = {self._envelop_alg.one_basis(): self.one()}
        self._cur_basis_inv = {self.one(): self._envelop_alg.one_basis()}
        self._gen_degrees = {}
        self._cur_num_gens = 0

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: g = lie_algebras.pwitt(GF(5), 5)
            sage: U = g.pbw_basis()
            sage: Z = U.center()
            sage: Z.indices()
            Basis indices of Center of Universal enveloping algebra of
             The 5-Witt Lie algebra over Finite Field of size 5 in the Poincare-Birkhoff-Witt basis
        """
        return "Basis indices of {}".format(self._center)

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: g = lie_algebras.pwitt(GF(5), 5)
            sage: U = g.pbw_basis()
            sage: Z = U.center()
            sage: I = Z.indices()
            sage: latex(I)
            B\left( Z\left( PBW\left( \mathcal{W}(5)_{\Bold{F}_{5}} \right) \right) \right)
        """
        from sage.misc.latex import latex
        return r"B\left( {} \right)".format(latex(self._center))

    def lift_on_basis(self, m):
        r"""
        Return the image of the basis element indexed by ``m`` in the
        universal enveloping algebra.

        EXAMPLES::

            sage: g = lie_algebras.Heisenberg(QQ, 3)
            sage: U = g.pbw_basis()
            sage: Z = U.center()
            sage: I = Z.indices()
            sage: z0 = I.monoid_generators()[0]
            sage: I._lift_map
            {0: {1: 1}}
            sage: I.lift_on_basis(z0)
            PBW['z']
            sage: I._lift_map
            {0: {1: 1}, 1: {PBW['z']: PBW['z']}}
            sage: I.lift_on_basis(z0^3)
            PBW['z']^3
            sage: I._lift_map
            {0: {1: 1}, 1: {PBW['z']: PBW['z']}}
            sage: I._construct_next_degree()
            sage: I._construct_next_degree()
            sage: I._lift_map
            {0: {1: 1},
             1: {PBW['z']: PBW['z']},
             2: {PBW['z']^2: PBW['z']^2},
             3: {PBW['z']^3: PBW['z']^3}}
            sage: I.lift_on_basis(z0^3)
            PBW['z']^3
        """
        while m not in self._cur_basis_inv:
            supp = m.support()
            # We might have not computed the correct degree, but we can lift the
            #   element if we have computed all of the corresponding generators.
            if all(i in self._gen_degrees and self._gen_degrees[i] in self._lift_map
                   for i in supp):
                ret = self._envelop_alg.one()
                divisors = [mp for mp in self._cur_basis_inv if mp.divides(m) and not mp.is_one()]
                while not m.is_one():
                    div = max(divisors, key=lambda elt: len(elt))
                    ls = self._cur_basis_inv[div]
                    deg = ls.length()
                    ret *= self._lift_map[deg][ls]
                    m = m // div
                    divisors = [mp for mp in divisors if mp.divides(m)]
                return ret
            self._construct_next_degree()
        ls = self._cur_basis_inv[m]
        deg = ls.length()
        return self._lift_map[deg][ls]

    def __iter__(self):
        r"""
        Iterate over ``self`` in degree increasing order.

        EXAMPLES::

            sage: g = lie_algebras.pwitt(GF(3), 6)
            sage: U = g.pbw_basis()
            sage: Z = U.center()
            sage: I = Z.indices()
            sage: it = iter(I)
            sage: [next(it) for _ in range(10)]
            [1, Z[0], Z[1], Z[2], Z[3], Z[4], Z[5], Z[6], Z[7], Z[0]^2]
        """
        yield self.one()  # start with the identity
        deg = 1
        while True:
            while deg not in self._lift_map:
                self._construct_next_degree()
            for ls in self._lift_map[deg]:
                yield self._cur_basis[ls]
            deg += 1

    def some_elements(self):
        r"""
        Return some elements of ``self``.

        EXAMPLES::

            sage: g = lie_algebras.pwitt(GF(3), 3)
            sage: U = g.pbw_basis()
            sage: Z = U.center()
            sage: I = Z.indices()
            sage: I.some_elements()
            [1, Z[0], Z[1], Z[2], Z[0]*Z[1]*Z[2], Z[0]*Z[2]^4, Z[0]^4*Z[1]^3]
        """
        it = iter(self)
        gens = [next(it) for _ in range(4)]
        # We construct it as a set in case we introduce duplicates.
        ret = set(gens)
        ret.update([self.prod(gens), gens[1] * gens[3]**4, gens[1]**4 * gens[2]**3])
        # Sort the output for uniqueness
        ret = sorted(ret, key=lambda m: (self.degree(m), m.to_word_list()))
        return ret

    def degree(self, m):
        r"""
        Return the degre of ``m`` in ``self``.

        EXAMPLES::

            sage: g = LieAlgebra(QQ, cartan_type=['E', 6])
            sage: U = g.pbw_basis()
            sage: Z = U.center()
            sage: I = Z.indices()
            sage: [I.degree(g) for g in I.monoid_generators()]
            [2, 5, 6, 8, 9, 12]
            sage: [(elt, I.degree(elt)) for elt in I.some_elements()]
            [(1, 0), (Z[0], 2), (Z[0]^2, 4), (Z[1], 5), (Z[0]^3*Z[1], 11),
             (Z[0]^10, 20), (Z[0]*Z[1]^4, 22)]
        """
        return ZZ.sum(e * self._gen_degrees[i] for i, e in m._monomial.items())

    def _construct_next_degree(self):
        r"""
        Construct the next elements of ``self`` for the next (uncomputed) degree.

        EXAMPLES::

            sage: g = lie_algebras.three_dimensional_by_rank(QQ, 2, 1)
            sage: U = g.pbw_basis()
            sage: Z = U.center()
            sage: I = Z.indices()
            sage: I._lift_map
            {0: {1: 1}}
            sage: I._construct_next_degree()
            sage: I._construct_next_degree()
            sage: I._construct_next_degree()
            sage: I._construct_next_degree()
            sage: I._lift_map
            {0: {1: 1}, 1: {}, 2: {}, 3: {}, 4: {}}
        """
        UEA = self._envelop_alg
        gens = UEA.algebra_generators()
        monoid = UEA.basis().keys()
        self._cur_deg += 1

        # We first update the lift map with all possible products
        # Note that we are using the fact that the elements are central
        #   so the product order doesn't matter.
        # Since we always update this, it is sufficient to compute it
        new_red = {}
        for i in range(1, self._cur_deg//2+1):
            for ls, lelt in self._lift_map[self._cur_deg-i].items():
                for rs, relt in self._lift_map[i].items():
                    supp = ls * rs
                    new_red[supp] = lelt * relt
                    mon = self._cur_basis[ls] * self._cur_basis[rs]
                    self._cur_basis[supp] = mon
                    self._cur_basis_inv[mon] = supp
        # TODO: Determine if we need to or benefit from another reduction of the new elements
        self._lift_map[self._cur_deg] = new_red

        # Determine the PBW elements of the current degree that are not reduced
        #   modulo the currently computed center.
        for exps in IntegerListsLex(n=self._cur_deg, length=len(gens)):
            elt = monoid.element_class(monoid, {k: p for k, p in zip(monoid._indices, exps) if p})
            if elt in new_red:  # already has a central element with this leading term
                continue
            # A new basis element to consider
            self._cur_vecs.append(UEA.monomial(elt))

        # Perform the centralization
        R = UEA.base_ring()
        vecs = list(self._cur_vecs)
        for g in gens:
            # TODO: We should hold onto previously computed values under the adjoint action
            ad = [g * v - v * g for v in vecs]
            # Compute the kernel
            supp = set()
            for v in ad:
                supp.update(v._monomial_coefficients)
            supp = sorted(supp, key=UEA._monomial_key, reverse=True)
            if not supp:  # no support for the image, so everything is in the kernel
                continue
            M = matrix(R, [[v[s] for v in ad] for s in supp])
            ker = M.right_kernel_matrix()
            vecs = [self._reduce(UEA.linear_combination((vecs[i], c) for i, c in kv.iteritems()))
                    for kv in ker.rows()]

        # Lastly, update the appropriate data
        if not vecs: # No new central elements, so nothing to do
            return
        new_gens = {}
        for v in vecs:
            v = self._reduce(v)  # possibly not needed to check this
            if not v:
                continue
            ls = v.trailing_support(key=UEA._monomial_key)
            self._cur_vecs.remove(UEA.monomial(ls))
            new_gens[ls] = self._reduce(v)
            assert (self._cur_num_gens not in self._gen_degrees
                    or self._gen_degrees[self._cur_num_gens] == self._cur_deg)
            self._gen_degrees[self._cur_num_gens] = self._cur_deg
            mon = self.gen(self._cur_num_gens)
            self._cur_basis[ls] = mon
            self._cur_basis_inv[mon] = ls
            self._cur_num_gens += 1
        self._lift_map[self._cur_deg].update(new_gens)

    def _reduce(self, vec):
        r"""
        Return the UEA vector ``vec`` by the currently computed center.

        EXAMPLES::

            sage: g = LieAlgebra(QQ, cartan_type=['A', 1])
            sage: U = g.pbw_basis()
            sage: Z = U.center()
            sage: I = Z.indices()
            sage: z0 = I.gen(0)
            sage: I._reduce(I.lift_on_basis(z0))
            0
            sage: max(I._lift_map)
            2
            sage: I._reduce(I.lift_on_basis(z0^2))
            4*PBW[alpha[1]]^2*PBW[-alpha[1]]^2
             + 2*PBW[alpha[1]]*PBW[alphacheck[1]]^2*PBW[-alpha[1]]
             + 1/4*PBW[alphacheck[1]]^4 - PBW[alphacheck[1]]^3
             + PBW[alphacheck[1]]^2
            sage: I._reduce(I.lift_on_basis(z0^2) - I.lift_on_basis(z0))
            4*PBW[alpha[1]]^2*PBW[-alpha[1]]^2
             + 2*PBW[alpha[1]]*PBW[alphacheck[1]]^2*PBW[-alpha[1]]
             + 1/4*PBW[alphacheck[1]]^4 - PBW[alphacheck[1]]^3
             + PBW[alphacheck[1]]^2
            sage: I._construct_next_degree()
            sage: I._construct_next_degree()
            sage: max(I._lift_map)
            4
            sage: I._reduce(I.lift_on_basis(z0^2) - I.lift_on_basis(z0))
            0
        """
        # This is replicating what SubmoduleWithBasis does
        ret = dict(vec._monomial_coefficients)
        for data in self._lift_map.values():
            for m, qv in data.items():
                if m not in ret:
                    continue
                iaxpy(-ret[m] / qv[m], qv._monomial_coefficients, ret)
        return self._envelop_alg._from_dict(ret, remove_zeros=False)


class SimpleLieCenterIndices(CenterIndices):
    r"""
    Set of basis indices for the center of a universal enveloping algebra of
    a simple Lie algebra.

    For more information, see
    :class:`~sage.algebras.lie_algebras.center_uea.CenterIndices`.
    """
    def __init__(self, center):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: g = LieAlgebra(QQ, cartan_type=['E', 6])
            sage: U = g.pbw_basis()
            sage: Z = U.center()
            sage: I = Z.indices()
            sage: TestSuite(I).run()
        """
        self._cartan_type = center._envelop_alg._g.cartan_type()
        r = self._cartan_type.rank()
        super().__init__(center, indices=FiniteEnumeratedSet(range(r)))
        W = CoxeterGroup(self._cartan_type)
        self._gen_degrees = dict(enumerate(W.degrees()))

    def __iter__(self):
        r"""
        Iterate over ``self`` in degree increasing order.

        EXAMPLES::

            sage: g = LieAlgebra(QQ, cartan_type=['E', 6])
            sage: U = g.pbw_basis()
            sage: Z = U.center()
            sage: I = Z.indices()
            sage: it = iter(I)
            sage: [next(it) for _ in range(10)]
            [1, Z[0], Z[0]^2, Z[1], Z[2], Z[0]^3, Z[0]*Z[1], Z[3], Z[0]*Z[2], Z[0]^4]
        """
        deg = 0
        n = len(self._gen_degrees)
        wts = sorted(self._gen_degrees.values(), reverse=True)
        while True:
            for exps in intvecwt_iterator(deg, wts):
                yield self.element_class(self, {n-1-i: e for i, e in enumerate(exps) if e})
            deg += 1


class CenterUEA(CombinatorialFreeModule):
    r"""
    The center of a universal enveloping algebra.

    .. TODO::

        Generalize this to be the centralizer of any set of the UEA.

    .. TODO::

        For characteristic `p > 0`, implement the `p`-center of a simple
        Lie algebra. See, e.g.,

        - Theorem 5.12 of [Motsak2006]_
        - http://www.math.kobe-u.ac.jp/icms2006/icms2006-video/slides/059.pdf

    EXAMPLES::

        sage: g = LieAlgebra(QQ, cartan_type=['A', 2])
        sage: U = g.pbw_basis()
        sage: Z = U.center()
        sage: B = Z.basis()
        sage: it = iter(B)
        sage: center_elts = [next(it) for _ in range(6)]; center_elts
        [1, Z[0], Z[1], Z[0]^2, Z[0]*Z[1], Z[1]^2]
        sage: elts = [U(v) for v in center_elts]  # long time
        sage: all(v * g == g * v for g in U.algebra_generators() for v in elts)  # long time
        True

    The Heisenberg Lie algebra `H_4` over a finite field; note the basis
    elements `b^p \in Z(U(H_4))` for the basis elements `b \in H_4`::

        sage: g = lie_algebras.Heisenberg(GF(3), 4)
        sage: U = g.pbw_basis()
        sage: Z = U.center()
        sage: B = Z.basis()
        sage: it = iter(B)
        sage: center_elts = [next(it) for _ in range(12)]; center_elts
        [1, Z[0], Z[0]^2, Z[0]^3, Z[1], Z[2], Z[3], Z[4], Z[5], Z[6], Z[7], Z[8]]
        sage: elts = [U(v) for v in center_elts]; elts
        [1, PBW['z'], PBW['z']^2, PBW['z']^3, PBW['p1']^3, PBW['p2']^3, PBW['p3']^3,
         PBW['p4']^3, PBW['q1']^3, PBW['q2']^3, PBW['q3']^3, PBW['q4']^3]
        sage: all(v * g == g * v for g in U.algebra_generators() for v in elts)
        True

    An example with a free 4-step nilpotent Lie algebras on 2 generators::

        sage: L = LieAlgebra(QQ, 2, step=4); L
        Free Nilpotent Lie algebra on 8 generators
         (X_1, X_2, X_12, X_112, X_122, X_1112, X_1122, X_1222) over Rational Field
        sage: U = L.pbw_basis()
        sage: Z = U.center()
        sage: it = iter(Z.basis())
        sage: center_elts = [next(it) for _ in range(10)]; center_elts
        [1, Z[0], Z[1], Z[2], Z[0]^2, Z[0]*Z[1], Z[0]*Z[2], Z[1]^2, Z[1]*Z[2], Z[2]^2]
        sage: elts = [U(v) for v in center_elts]; elts
        [1, PBW[(1, 1, 1, 2)], PBW[(1, 1, 2, 2)], PBW[(1, 2, 2, 2)], PBW[(1, 1, 1, 2)]^2,
         PBW[(1, 1, 1, 2)]*PBW[(1, 1, 2, 2)], PBW[(1, 1, 1, 2)]*PBW[(1, 2, 2, 2)],
         PBW[(1, 1, 2, 2)]^2, PBW[(1, 1, 2, 2)]*PBW[(1, 2, 2, 2)], PBW[(1, 2, 2, 2)]^2]
        sage: all(v * g == g * v for g in U.algebra_generators() for v in elts)
        True

    Using the Engel Lie algebra::

        sage: L.<X,Y,Z> = LieAlgebra(QQ, {('X','Y'): {'Z': 1}}, nilpotent=True)
        sage: U = L.pbw_basis()
        sage: Z = U.center()
        sage: it = iter(Z.basis())
        sage: center_elts = [next(it) for _ in range(6)]; center_elts
        [1, Z[0], Z[0]^2, Z[0]^3, Z[0]^4, Z[0]^5]
        sage: elts = [U(v) for v in center_elts]; elts
        [1, PBW['Z'], PBW['Z']^2, PBW['Z']^3, PBW['Z']^4, PBW['Z']^5]
        sage: all(v * g == g * v for g in U.algebra_generators() for v in elts)
        True
    """
    def __init__(self, g, UEA):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: g = LieAlgebra(ZZ['t'].fraction_field(), cartan_type=['D', 4])
            sage: U = g.pbw_basis()
            sage: Z = U.center()
            sage: TestSuite(Z).run()

            sage: g = lie_algebras.Heisenberg(GF(3), 4)
            sage: U = g.pbw_basis()
            sage: Z = U.center()
            sage: TestSuite(Z).run()
        """
        if g not in FiniteDimensionalLieAlgebrasWithBasis:
            raise NotImplementedError("only implemented for finite dimensional Lie algebras with a distinguished basis")

        R = UEA.base_ring()
        if R not in Fields():
            raise NotImplementedError("only implemented for the base ring a field")

        self._g = g
        self._envelop_alg = UEA
        if (self._g in KacMoodyAlgebras
            and self._g.cartan_type().is_finite()
            and R.characteristic() == 0):
            indices = SimpleLieCenterIndices(self)
        else:
            indices = CenterIndices(self)
        category = UEA.category()
        base = category.base()
        category = GradedAlgebrasWithBasis(base).Commutative() | category.Subobjects()
        CombinatorialFreeModule.__init__(self, R, indices, category=category,
                                         prefix='', bracket=False, latex_bracket=False,
                                         sorting_key=self._sorting_key)
        self.lift.register_as_coercion()

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: g = LieAlgebra(QQ, cartan_type=['A',2])
            sage: U = g.pbw_basis()
            sage: U.center()
            Center of Universal enveloping algebra of Lie algebra of ['A', 2]
             in the Chevalley basis in the Poincare-Birkhoff-Witt basis
        """
        return "Center of " + repr(self._envelop_alg)

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: g = lie_algebras.pwitt(GF(5), 5)
            sage: U = g.pbw_basis()
            sage: Z = U.center()
            sage: latex(Z)
            Z\left( PBW\left( \mathcal{W}(5)_{\Bold{F}_{5}} \right) \right)
        """
        from sage.misc.latex import latex
        return r"Z\left( {} \right)".format(latex(self._envelop_alg))

    def _sorting_key(self, m):
        r"""
        Return a key for ``m`` used in sorting elements of ``self``.

        EXAMPLES::

            sage: g = LieAlgebra(QQ, cartan_type=['A', 2])
            sage: U = g.pbw_basis()
            sage: Z = U.center()
            sage: z0, z1 = Z.algebra_generators()
            sage: z0 * z1  # indirect doctest
            Z[0]*Z[1]
            sage: z1^2 + z0*z1 + z0^2  # indirect doctest
            Z[0]^2 + Z[0]*Z[1] + Z[1]^2
            sage: z1^3 + z0*z1^2 + z0^2*z1 + z0^3  # indirect doctest
            Z[0]^3 + Z[0]^2*Z[1] + Z[0]*Z[1]^2 + Z[1]^3
        """
        return (-m.length(), m.to_word_list())

    @cached_method
    def algebra_generators(self):
        r"""
        Return the algebra generators of ``self``.

        .. WARNING::

            When the universal enveloping algebra is not known to have
            a finite generating set, the generating set will be the basis
            of ``self`` in a degree (weakly) increasing order indexed by
            `\ZZ_{\geq 0}`. In particular, the `0`-th generator will be
            the multiplicative identity `1`.

        EXAMPLES::

            sage: g = lie_algebras.Heisenberg(QQ, 3)
            sage: U = g.pbw_basis()
            sage: Z = U.center()
            sage: Z.algebra_generators()[0]
            1
            sage: Z.algebra_generators()[1]
            Z[0]

            sage: g = LieAlgebra(QQ, cartan_type=['G', 2])
            sage: U = g.pbw_basis()
            sage: Z = U.center()
            sage: Z.algebra_generators()
            Finite family {0: Z[0], 1: Z[1]}
        """
        mon_gens = self._indices.monoid_generators()
        if mon_gens.cardinality() == float("inf"):
            return Family(NonNegativeIntegers(), lambda m: self.monomial(self._indices.unrank(m)))
        return Family({i: self.monomial(mon_gens[i]) for i in mon_gens.keys()})

    @cached_method
    def one_basis(self):
        r"""
        Return the basis index of `1` in ``self``.

        EXAMPLES::

            sage: g = LieAlgebra(QQ['t'].fraction_field(), cartan_type=['B', 5])
            sage: U = g.pbw_basis()
            sage: Z = U.center()
            sage: ob = Z.one_basis(); ob
            1
            sage: ob.parent() is Z.indices()
            True
        """
        return self._indices.one()

    def ambient(self):
        r"""
        Return the ambient algebra of ``self``.

        EXAMPLES::

            sage: g = LieAlgebra(GF(5), cartan_type=['A', 2])
            sage: U = g.pbw_basis()
            sage: Z = U.center()
            sage: Z.ambient() is U
            True
        """
        return self._envelop_alg

    def product_on_basis(self, left, right):
        r"""
        Return the product of basis elements indexed by ``left`` and ``right``.

        EXAMPLES::

            sage: g = LieAlgebra(QQ, cartan_type=['E', 6])
            sage: U = g.pbw_basis()
            sage: Z = U.center()
            sage: mg = Z.indices().monoid_generators()
            sage: Z.product_on_basis(mg[1]*mg[2], mg[0]*mg[1]^3*mg[2]*mg[3]^3)
            Z[0]*Z[1]^4*Z[2]^2*Z[3]^3
        """
        return self.monomial(left * right)

    def degree_on_basis(self, m):
        r"""
        Return the degree of the basis element indexed by ``m`` in ``self``.

        EXAMPLES::

            sage: g = LieAlgebra(QQ, cartan_type=['E', 6])
            sage: U = g.pbw_basis()
            sage: Z = U.center()
            sage: I = Z.indices()
            sage: it = iter(I)
            sage: supports = [next(it) for _ in range(10)]; supports
            [1, Z[0], Z[0]^2, Z[1], Z[2], Z[0]^3, Z[0]*Z[1], Z[3], Z[0]*Z[2], Z[0]^4]
            sage: [Z.degree_on_basis(m) for m in supports]
            [0, 2, 4, 5, 6, 6, 7, 8, 8, 8]
        """
        return self._indices.degree(m)

    @lazy_attribute
    def lift(self):
        r"""
        The lift map from ``self`` to the universal enveloping algebra.

        EXAMPLES::

            sage: g = LieAlgebra(QQ, cartan_type=['A', 1])
            sage: U = g.pbw_basis()
            sage: Z = U.center()
            sage: gens = Z.algebra_generators()
            sage: U(gens[0]^2 + gens[0])
            4*PBW[alpha[1]]^2*PBW[-alpha[1]]^2
             + 2*PBW[alpha[1]]*PBW[alphacheck[1]]^2*PBW[-alpha[1]]
             + 1/4*PBW[alphacheck[1]]^4 - PBW[alphacheck[1]]^3
             - 2*PBW[alpha[1]]*PBW[-alpha[1]] + 1/2*PBW[alphacheck[1]]^2
             + PBW[alphacheck[1]]
            sage: U(-1/4*gens[0]) == U.casimir_element()
            True
        """
        # This is correct if we are using key=self._envelop_alg._monomial_key,
        #   but we are currently unable to pass such an option.
        return self.module_morphism(self._indices.lift_on_basis, codomain=self._envelop_alg, unitriangular='upper')

    def retract(self, elt):
        r"""
        The retraction map to ``self`` from the universal enveloping algebra.

        .. TODO::

            Implement a version of this that checks if the leading term of
            ``elt`` is divisible by a product of all of the currently known
            generators in order to avoid constructing the full centralizer
            of larger degrees than needed.

        EXAMPLES::

            sage: g = lie_algebras.Heisenberg(QQ, 3)
            sage: U = g.pbw_basis()
            sage: Z = U.center()
            sage: z0 = Z.algebra_generators()[1]; z0
            Z[0]
            sage: Z.retract(U(z0^2) - U(3*z0))
            Z[0]^2 - 3*Z[0]

            sage: g = LieAlgebra(QQ, cartan_type=['A', 2])
            sage: U = g.pbw_basis()
            sage: Z = U.center()
            sage: z0, z1 = Z.algebra_generators()
            sage: Z.retract(U(z0*z0) - U(z1))  # long time
            Z[0]^2 - Z[1]
            sage: zc = Z.retract(U.casimir_element()); zc
            -1/3*Z[0]
            sage: U(zc) == U.casimir_element()
            True
        """
        # This should work except it needs the monomials of the PBW basis to be
        # comparable. However, this does not work for, e.g., Lie algebras
        # in the Chevalley basis as ee are unable to pass a key for the
        # module morphism. Additionally, the implementation below does more
        # operations in-place than the module morphism.
        #return self.lift.section()
        UEA = self._envelop_alg
        elt = UEA(elt)
        # We manipulate the dictionary (in place) to avoid creating elements
        data = elt.monomial_coefficients(copy=True)
        indices = self._indices
        ret = {}
        while data:
            lm = min(data, key=UEA._monomial_key)
            while indices._cur_deg < UEA.degree_on_basis(lm):
                indices._construct_next_degree()
            ind = indices._cur_basis[lm]
            other = indices.lift_on_basis(ind).monomial_coefficients(copy=False)
            coeff = data[lm] / other[lm]
            ret[ind] = coeff
            iaxpy(-coeff, other, data)
        return self.element_class(self, ret)
