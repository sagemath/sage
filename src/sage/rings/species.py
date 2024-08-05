from itertools import chain
from sage.categories.graded_algebras_with_basis import GradedAlgebrasWithBasis
from sage.categories.infinite_enumerated_sets import InfiniteEnumeratedSets
from sage.categories.monoids import Monoids
from sage.categories.sets_with_grading import SetsWithGrading
from sage.combinat.free_module import CombinatorialFreeModule
from sage.combinat.integer_vector import IntegerVectors
from sage.groups.perm_gps.constructor import PermutationGroupElement
from sage.groups.perm_gps.permgroup import PermutationGroup, PermutationGroup_generic
from sage.groups.perm_gps.permgroup_named import SymmetricGroup
from sage.libs.gap.libgap import libgap
from sage.misc.cachefunc import cached_method
from sage.monoids.indexed_free_monoid import IndexedFreeAbelianMonoid, IndexedFreeAbelianMonoidElement, IndexedMonoid
from sage.rings.integer_ring import ZZ
from sage.structure.element import Element, parent
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.sets.finite_enumerated_set import FiniteEnumeratedSet

GAP_FAIL = libgap.eval('fail')


def _is_conjugate(G, H1, H2):
    r"""
    Test if ``H1`` and ``H2`` are conjugate subgroups in ``G``.

    EXAMPLES::

        sage: G = SymmetricGroup(3)
        sage: H1 = PermutationGroup([(1,2)])
        sage: H2 = PermutationGroup([(2,3)])
        sage: from sage.rings.species import _is_conjugate
        sage: _is_conjugate(G, H1, H2)
        True
    """
    return GAP_FAIL != libgap.RepresentativeAction(G, H1, H2)


class ElementCache():
    def __init__(self):
        r"""
        Class for caching elements.
        """
        self._cache = dict()

    def _cache_get(self, elm):
        r"""
        Return the cached element, or create it
        if it doesn't exist.

        ``elm`` must implement the following methods:
          ``_element_key`` - hashable type for dict lookup.
          ``__eq__`` - to compare two elements.
        """
        key = elm._element_key()
        if key in self._cache:
            lookup = self._cache[key]
            for other_elm in lookup:
                if elm == other_elm:
                    return other_elm
            lookup.append(elm)
        else:
            self._cache[key] = [elm]
        return elm


class ConjugacyClassOfDirectlyIndecomposableSubgroups(Element):
    def __init__(self, parent, C):
        r"""
        A conjugacy class of directly indecomposable subgroups.
        """
        Element.__init__(self, parent)
        self._C = C
        self._sorted_orbits = tuple(sorted(len(orbit) for orbit in C.orbits()))
        self._order = C.order()

    def __hash__(self):
        return hash(self._C)

    def _repr_(self):
        r"""
        Return a string representation of ``self``.
        """
        return f"{self._C.gens_small()}"

    def __eq__(self, other):
        r"""
        Return whether ``self`` is equal to ``other``.

        ``self`` is equal to ``other`` if they have the same degree (say `n`)
        and order and are conjugate within `S_n`.
        """
        return (isinstance(other, ConjugacyClassOfDirectlyIndecomposableSubgroups)
                and self._C.degree() == other._C.degree()
                and  _is_conjugate(SymmetricGroup(self._C.degree()),
                                   self._C, other._C))


class ConjugacyClassesOfDirectlyIndecomposableSubgroups(UniqueRepresentation, Parent):
    def __init__(self):
        r"""
        Conjugacy classes of directly indecomposable subgroups.
        """
        Parent.__init__(self, category=InfiniteEnumeratedSets())
        ElementCache.__init__(self)

    def _element_constructor_(self, x):
        r"""
        ``x`` is an element of ``self`` or a group `H` such that
        `H` is directly indecomposable.

        EXAMPLES::

            sage: G = PermutationGroup([[(1,3),(4,7)], [(2,5),(6,8)], [(1,4),(2,5),(3,7)]])
            sage: A = AtomicSpecies(1)
            sage: A(G)
            {[(5,6)(7,8), (1,2)(5,7)(6,8), (1,2)(3,4)]: (8,)}
        """
        if parent(x) == self:
            return x
        if isinstance(x, PermutationGroup_generic):
            if len(x.disjoint_direct_product_decomposition()) > 1:
                raise ValueError(f"{x} is not directly indecomposable")
            pi = self.canonical_label(x)
            return self.element_class(self, PermutationGroup(gap_group=libgap.ConjugateGroup(x, pi)))
        raise ValueError(f"unable to convert {x} to {self}")

    def _repr_(self):
        return "Infinite set of conjugacy classes of directly indecomposable subgroups"

    @cached_method
    def canonical_label(self, H):
        r"""
        Return the permutation group element `g` such that
        `g^{-1} H g` is a canonical representative.
        """
        return PermutationGroupElement([e for o in sorted([sorted(orbit) for orbit in H.orbits()], key=len) for e in o], check=False)

    Element = ConjugacyClassOfDirectlyIndecomposableSubgroups


class AtomicSpeciesElement(Element):
    def __init__(self, parent, dis, domain_partition):
        r"""
        Initialize an atomic species.

        INPUT:

            - ``dis``, a :class:`ConjugacyClassOfDirectlyIndecomposableSubgroups`

            - ``domain_partition``, a dict representing the
              assignment of each element of the domain of ``dis`` to
              a "variable".
        """
        # Notice that a molecular species (and therefore an atomic species)
        # must be centered on a multicardinality, otherwise it wouldn't be
        # molecular. So it is kind of redundant to say "atomic species
        # centered on a given multicardinality."
        Element.__init__(self, parent)
        self._dis = dis
        self._dompart = domain_partition
        L = [0 for _ in range(self.parent()._k)]
        for v in self._dompart.values():
            L[v - 1] += 1
        self._mc = tuple(L)
        self._tc = sum(self._mc)

    def _element_key(self):
        r"""
        Return a lookup key for ``self``.
        """
        return self._mc, self._dis._order, self._dis._sorted_orbits

    def __hash__(self):
        r"""
        Return the hash of the atomic species.
        """
        return hash((self._mc, self._dis))

    def __eq__(self, other):
        r"""
        Two atomic species are equal if the underlying groups are conjugate,
        and their multicardinalities are equal.
        """
        if parent(self) != parent(other):
            return False
        return self._mc == other._mc and self._dis == other._dis

    def _repr_(self):
        r"""
        Return a string representation of ``self``.
        """
        return "{" + f"{self._dis}: {self._mc}" + "}"


class AtomicSpecies(UniqueRepresentation, Parent, ElementCache):
    def __init__(self, k):
        r"""
        Infinite set of `k`-variate atomic species graded by
        integer vectors of length `k`.
        """
        category = SetsWithGrading().Infinite()
        Parent.__init__(self, category=category)
        ElementCache.__init__(self)
        self._k = k
        self._grading_set = IntegerVectors(length=k)
        self._dis_ctor = ConjugacyClassesOfDirectlyIndecomposableSubgroups()

    def __hash__(self):
        r"""
        Return a hash for ``self``.
        """
        return hash(self._k)

    @cached_method
    def an_element(self):
        """
        Return an element of ``self``.
        """
        G = SymmetricGroup(self._k).young_subgroup([1] * self._k)
        m = {e: i for i, e in enumerate(range(1, self._k + 1), 1)}
        return self._element_constructor_((G, m))

    def _normalize(self, H, M):
        r"""
        Normalize the domain of `H` and return `H` and the domain partition.
        """
        # TODO: Complete the documentation.
        if set(H.domain()) != set(M.keys()):
            raise ValueError(f"Keys of {M} do not match with domain of {H} (= {H.domain()})")
        if not set(M.values()).issubset(range(1, self._k + 1)):
            raise ValueError(f"Values of {M} must be in the range [1, {self._k}]")
        # normalize domain to {1..n}
        if sorted(M.keys()) == list(range(1, H.degree() + 1)):
            return H, M
        mapping = {v: i for i, v in enumerate(H.domain(), 1)}
        normalized_gens = [[tuple(mapping[x] for x in cyc)
                            for cyc in gen.cycle_tuples()]
                           for gen in H.gens()]
        P = PermutationGroup(gens=normalized_gens)
        # Fix for SymmetricGroup(0)
        if H.degree() == 0:
            P = SymmetricGroup(0)
        # create domain partition
        dompart = {mapping[k]: v for k, v in M.items()}
        return P, dompart

    def _element_constructor_(self, x):
        r"""
        Construct the `k`-variate atomic species with the given data.

        INPUT:

        - ``x`` can be any of the following:
            - an element of ``self``.
            - a tuple ``(H, M)`` where `H` is the permutation group
              representation for the atomic species and `M` is a
              ``dict`` mapping each element of the domain of `H` to
              integers in `\{ 1 \ldots k \}`, representing the set to
              which the element belongs.
            - if `k=1`, i.e. we are working with univariate atomic
              species, the mapping `M` may be omitted and just the
              group `H` may be passed.

        """
        if parent(x) == self:
            return x
        H, M = None, None
        if isinstance(x, tuple):
            H, M = x
        elif self._k == 1:
            if isinstance(x, PermutationGroup_generic):
                H = x
                M = {e: 1 for e in H.domain()}
            else:
                raise ValueError(f"{x} must be a permutation group")
        else:
            raise ValueError(f"{x} must be a tuple for multivariate species")
        H_norm, dompart = self._normalize(H, M)
        dis_elm = self._dis_ctor(H_norm)
        # Python sorts are stable so we can use sorted twice without worrying about a change
        # in the output.
        pi = self._dis_ctor.canonical_label(H_norm)
        dompart_norm = {pi(k): v for k, v in dompart.items()}
        elm = self.element_class(self, dis_elm, dompart_norm)
        return self._cache_get(elm)

    def __getitem__(self, x):
        r"""
        Call ``_element_constructor_`` on ``x``.
        """
        return self._element_constructor_(x)

    def __contains__(self, x):
        r"""
        Return if ``x`` is in ``self``.
        """
        if parent(x) == self:
            return True
        H, M = x
        return (set(H.domain()) == set(M.keys())
                and set(M.values()).issubset(range(1, self._k + 1))
                and len(H.disjoint_direct_product_decomposition()) == 1)

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        TESTS::

            sage: At1 = AtomicSpecies(1)
            sage: At1
            Infinite set of 1-variate atomic species
            sage: At2 = AtomicSpecies(2)
            sage: At2
            Infinite set of 2-variate atomic species
        """
        return f"Infinite set of {self._k}-variate atomic species"

    Element = AtomicSpeciesElement

class MolecularSpecies(IndexedFreeAbelianMonoid, ElementCache):
    @staticmethod
    def __classcall__(cls, indices, prefix, **kwds):
        return super(IndexedMonoid, cls).__classcall__(cls, indices, prefix, **kwds)

    def __init__(self, indices, prefix, **kwds):
        category = Monoids() & InfiniteEnumeratedSets()
        IndexedFreeAbelianMonoid.__init__(self, indices, prefix=prefix, category=category, **kwds)
        ElementCache.__init__(self)

    def _element_constructor_(self, x=None):
        elm = self._cache_get(super()._element_constructor_(x))
        if elm._group is None:
            elm._group_constructor()
        return elm

    @cached_method
    def one(self):
        elm = super().one()
        elm._group = SymmetricGroup(0)
        return elm

    def gen(self, x):
        elm = self._cache_get(super().gen(x))
        if elm._group is None:
            elm._group_constructor()
        return elm

    class Element(IndexedFreeAbelianMonoidElement):
        def __init__(self, F, x):
            super().__init__(F, x)
            self._group = None

        def _group_constructor(self):
            r"""
            Construct the group of ``self``.
            """
            if self._group is None:
                self._group = SymmetricGroup(0)
            for A, p in self._monomial.items():
                curgrp = self._groupexp(A._dis._C, p)
                self._group = self._groupmul(self._group, curgrp)

        def _groupmul(self, G1, G2):
            r"""
            Multiply two groups.
            """
            # Would be great to be able to cache the intermediate results
            if G1.degree() + G2.degree() == 0:
                return SymmetricGroup(0)
            if G1.degree() == 0:
                return G2
            if G2.degree() == 0:
                return G1
            gens1 = G1.gens()
            if len(gens1) > 50:
                gens1 = G1.gens_small()
            gens2 = G2.gens()
            if len(gens2) > 50:
                gens1 = G2.gens_small()
            # always loop over the smaller gens list
            if len(gens1) < len(gens2):
                gens1, gens2 = gens2, gens1
                G1, G2 = G2, G1
            gens = list(gens1)
            for gen in gens2:
                gens.append([tuple(G1.degree() + k for k in cyc) for cyc in gen.cycle_tuples()])
            return PermutationGroup(gens, domain=range(1, G1.degree() + G2.degree() + 1))

        def _groupexp(self, G, n):
            r"""
            Exponentiate a group.
            """
            grp = SymmetricGroup(0)
            while n > 0:
                if n % 2 == 1:
                    grp = self._groupmul(grp, G)
                G = self._groupmul(G, G)
                n //= 2
            return grp

        # TODO: didn't override __floordiv__
        def _mul_(self, other):
            res = super()._mul_(other)
            elm = self.parent()._cache_get(res)
            if self._group is None:
                self._group = SymmetricGroup(0)
            if elm._group is None:
                # Multiply two groups
                elm._group = elm._groupmul(self._group, other._group)
            return elm

        def __pow__(self, n):
            res = super().__pow__(n)
            elm = self.parent()._cache_get(res)
            if self._group is None:
                self._group = SymmetricGroup(0)
            if elm._group is None:
                # Exponentiate the group
                elm._group = elm._groupexp(self._group, n)
            return elm

        def _element_key(self):
            return self

        def grade(self):
            r"""
            Return the grade of ``self``.
            """
            return sum(A._tc * p for A, p in self._monomial.items())

        def domain(self):
            r"""
            Return the domain of ``self``.
            """
            return FiniteEnumeratedSet(range(1, self.grade() + 1))


class PolynomialSpecies(CombinatorialFreeModule):
    def __init__(self, k, base_ring=ZZ):
        r"""
        Ring of `k`-variate polynomial (virtual) species.

        TESTS::

            sage: P = PolynomialSpecies(1)
            sage: TestSuite(P).run()
            sage: P2 = PolynomialSpecies(2)
            sage: TestSuite(P2).run()
        """
        # should we pass a category to basis_keys?
        basis_keys = MolecularSpecies(AtomicSpecies(k),
                                              prefix='', bracket=False)
        category = GradedAlgebrasWithBasis(base_ring).Commutative()
        CombinatorialFreeModule.__init__(self, base_ring,
                                        basis_keys=basis_keys,
                                        category=category,
                                        element_class=self.Element,
                                        prefix='', bracket=False)
        self._k = k
        self._atomic_basis = basis_keys.indices()

    def degree_on_basis(self, m):
        r"""
        Return the degree of the molecular species indexed by ``m``.
        """
        return m.grade()

    def _project(self, H, f, part):
        r"""
        Project `H` onto a subset ``part`` of its domain.

        ``part`` must be a union of cycles, but this is not checked.
        """
        restricted_gens = [[cyc for cyc in gen.cycle_tuples() if cyc[0] in part] for gen in H.gens()]
        mapping = {p: f[p] for p in part}
        return PermutationGroup(gens=restricted_gens, domain=part), mapping

    def _element_constructor_(self, x):
        r"""
        Construct the `k`-variate molecular species with the given data.

        INPUT:

        - ``x`` can be any of the following:
            - an element of ``self``.
            - a tuple ``(H, M)`` where `H` is the permutation group
              representation for the species and `M` is a ``dict``
              mapping each element of the domain of `H` to integers
              in `\{ 1 \ldots k \}`, representing the set to which
              the element belongs.
            - if `k=1`, i.e. we are working with univariate species,
              the mapping `M` may be omitted and just the group `H`
              may be passed.

        """
        if parent(x) == self:
            return x
        H, M = None, None
        if isinstance(x, tuple):
            H, M = x
        elif self._k == 1:
            if isinstance(x, PermutationGroup_generic):
                H = x
                M = {e: 1 for e in H.domain()}
            else:
                raise ValueError(f"{x} must be a permutation group")
        else:
            raise ValueError(f"{x} must be a tuple for multivariate species")
        domain_partition = H.disjoint_direct_product_decomposition()
        term = self._indices.one()
        for part in domain_partition:
            term *= self._indices.gen(self._project(H, M, part))
        return self._from_dict({term: 1})

    def __getitem__(self, x):
        r"""
        Calls _element_constructor_ on x.

        TESTS::

            sage: P = PolynomialSpecies(1)
            sage: At1 = AtomicSpecies(1)
            sage: At1(SymmetricGroup(1)).rename("X")
            sage: X2 = SymmetricGroup(2).young_subgroup([1, 1])
            sage: P[X2]
            X^2
            sage: P2 = PolynomialSpecies(2)
            sage: At2 = AtomicSpecies(2)
            sage: At2((SymmetricGroup(1), {1: 1})).rename("X")
            sage: At2((SymmetricGroup(1), {1: 2})).rename("Y")
            sage: XY = (SymmetricGroup(2).young_subgroup([1, 1]), {1: 1, 2: 2})
            sage: P2[XY]
            X*Y
        """
        return self._element_constructor_(x)

    @cached_method
    def one_basis(self):
        r"""
        Returns SymmetricGroup(0), which indexes the one of this algebra,
        as per :meth:`AlgebrasWithBasis.ParentMethods.one_basis`.

        EXAMPLES::

            sage: P = PolynomialSpecies(1)
            sage: P.one_basis()
            1
            sage: P2 = PolynomialSpecies(2)
            sage: P2.one_basis()
            1
        """
        return self._indices.one()

    @cached_method
    def an_element(self):
        """
        Return an element of ``self``.

            sage: P = PolynomialSpecies(1)
            sage: P.an_element()
            1
            sage: P2 = PolynomialSpecies(2)
            sage: P2.an_element()
            1
        """
        return self.one()

    def product_on_basis(self, H, K):
        r"""
        Return the product of the basis elements indexed by `H` and `K`.

        EXAMPLES::

            sage: A = AtomicSpecies(1)
            sage: A(SymmetricGroup(1)).rename("X")
            sage: [A(SymmetricGroup(n)).rename(f"E_{n}") for n in range(2, 5)]
            [None, None, None]
            sage: [A(CyclicPermutationGroup(n)).rename(f"C_{n}") for n in range(3, 5)]
            [None, None]
            sage: P = PolynomialSpecies(1)
            sage: L1 = [P(H) for H in SymmetricGroup(3).conjugacy_classes_subgroups()]
            sage: L2 = [P(H) for H in SymmetricGroup(2).conjugacy_classes_subgroups()]
            sage: matrix([[F * G for F in L1] for G in L2])
            [    X^5 X^3*E_2 C_3*X^2 E_3*X^2]
            [X^3*E_2 X*E_2^2 C_3*E_2 E_3*E_2]
        """
        return self._from_dict({H * K: 1})

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: P = PolynomialSpecies(1)
            sage: P
            Ring of 1-variate virtual species
            sage: P2 = PolynomialSpecies(2)
            sage: P2
            Ring of 2-variate virtual species
        """
        return f"Ring of {self._k}-variate virtual species"

    class Element(CombinatorialFreeModule.Element):
        def is_virtual(self):
            return any(x < 0 for x in self.coefficients(sort=False))

        def is_molecular(self):
            return len(self.coefficients(sort=False)) == 1 and self.coefficients(sort=False)[0] == 1

        def is_atomic(self):
            return self.is_molecular() and len(self.support()[0]) == 1

        def __call__(self, *args):
            # should this also have a base ring check or is it coerced?
            # the usage of type(self)(...), is it correct?
            if len(args) != self.parent()._k:
                raise ValueError(f"Number of args (= {len(args)}) must equal arity of self (= {len(args)})")
            if any(not isinstance(arg, PolynomialSpecies.Element) for arg in args):
                raise ValueError("All args must be elements of PolynomialSpecies")
            if self.parent()._k > 1 or max(arg.parent()._k for arg in args) > 1:
                raise NotImplementedError("Only univariate species are supported")
            if self.is_virtual() or any(arg.is_virtual() for arg in args):
                raise NotImplementedError(f"Only non-virtual species are supported")
            # Now we only have non-virtual univariate species
            res = 0
            Gm = args[0].monomial_coefficients()
            for M, fM in self.monomial_coefficients().items():
                term = 0
                powvecs = IntegerVectors(M.grade(), len(Gm))
                for vec in powvecs:
                    coeff = 1
                    for fN, exponent in zip(Gm.values(), vec):
                        coeff *= fN ** exponent
                    N_list = list(chain.from_iterable([[N for _ in range(c)] for N, c in zip(Gm.keys(), vec)]))
                    R = self.parent()(wreath_imprimitive_general(N_list, M))
                    term += coeff * R
                res += fM * term
            return res


def wreath_imprimitive_general(G_list, H):
    r"""
    H([G[0]] - [G[1]] - ... - [G[m]])
    All members of MolecularSpecies.
    """
    if len(G_list) != H.grade():
        raise ValueError(f"length of G_list (= {len(G_list)}) must be equal to degree of H (= {H.grade()})")

    gens, dlist = [], []
    dsum = 0
    for G in G_list:
        dlist.append(list(range(dsum + 1, dsum + G.grade() + 1)))
        dsum += G.grade()

    Hgens = H._group.gens()
    # First find gens of H
    for gen in Hgens:
        # each cycle must be homogenous in the degrees of groups shifted
        # bad things happen if this isn't checked
        # example: partitional_composition(E2, X+X^2)
        # look at the structure with 1 X and 1 X^2, like {(1), (2,3)}
        all_homogenous = True
        for cyc in gen.cycle_tuples():
            cyclens = [len(dlist[x - 1]) for x in cyc]
            all_homogenous &= min(cyclens) == max(cyclens)
        if not all_homogenous:
            continue
        perm = libgap.Flat(libgap.Permuted(dlist, gen))
        gens.append(perm)

    # Then find gens of G_i
    for i in range(len(G_list)):
        Ggens = G_list[i]._group.gens()
        for gen in Ggens:
            images = libgap.MappingPermListList(G_list[i].domain().list(), dlist[i])
            gens.append(gen ** images)

    # Finally create group
    return PermutationGroup(gens, domain = range(1, dsum + 1))
