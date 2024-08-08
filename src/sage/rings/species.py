from itertools import accumulate, chain, combinations
from sage.categories.cartesian_product import cartesian_product
from sage.categories.graded_algebras_with_basis import GradedAlgebrasWithBasis
from sage.categories.infinite_enumerated_sets import InfiniteEnumeratedSets
from sage.categories.monoids import Monoids
from sage.categories.sets_with_grading import SetsWithGrading
from sage.combinat.free_module import CombinatorialFreeModule
from sage.combinat.integer_vector import IntegerVectors
from sage.combinat.partition import Partitions
from sage.combinat.permutation import Permutations
from sage.groups.perm_gps.constructor import PermutationGroupElement
from sage.groups.perm_gps.permgroup import PermutationGroup, PermutationGroup_generic
from sage.groups.perm_gps.permgroup_named import SymmetricGroup
from sage.libs.gap.libgap import libgap
from sage.misc.cachefunc import cached_method
from sage.misc.lazy_attribute import lazy_attribute
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
        Return the cached element, or add it
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
            elm = self._canonical_label(elm)
            lookup.append(elm)
        else:
            elm = self._canonical_label(elm)
            self._cache[key] = [elm]
        return elm

    @cached_method
    def _canonical_label(self, elm):
        return elm

class ConjugacyClassOfDirectlyIndecomposableSubgroups(Element):
    def __init__(self, parent, C):
        r"""
        A conjugacy class of directly indecomposable subgroups.
        """
        Element.__init__(self, parent)
        self._C = C
        self._sorted_orbits = sorted([sorted(orbit) for orbit in C.orbits()], key=len)
        self._orbit_lens = tuple(len(orbit) for orbit in self._sorted_orbits)
        self._order = C.order()

    @lazy_attribute
    def _canonicalizing_perm(self):
        return PermutationGroupElement([e for o in self._sorted_orbits for e in o], check=False)

    def __hash__(self):
        return hash(self._C)

    def _repr_(self):
        r"""
        Return a string representation of ``self``.
        """
        return f"{self._C.gens_small()}"

    def _element_key(self):
        return self._C.degree(), self._order, self._orbit_lens

    def __eq__(self, other):
        r"""
        Return whether ``self`` is equal to ``other``.

        ``self`` is equal to ``other`` if they have the same degree (say `n`)
        and order and are conjugate within `S_n`.
        """
        return (isinstance(other, ConjugacyClassOfDirectlyIndecomposableSubgroups)
                and self._C.degree() == other._C.degree()
                and (self._C == other._C
                     or _is_conjugate(SymmetricGroup(self._C.degree()),
                                      self._C, other._C)))


class ConjugacyClassesOfDirectlyIndecomposableSubgroups(UniqueRepresentation, Parent, ElementCache):
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
            elm = self.element_class(self, x)
            return self._cache_get(elm)
        raise ValueError(f"unable to convert {x} to {self}")

    def _repr_(self):
        return "Infinite set of conjugacy classes of directly indecomposable subgroups"

    @cached_method
    def canonical_label(self, elm):
        r"""

        EXAMPLES::

            sage: G = PermutationGroup([[(1,3),(4,7)], [(2,5),(6,8)], [(1,4),(2,5),(3,7)]])
            sage: A = AtomicSpecies(1)
            sage: A(G)
            {[(5,6)(7,8), (1,2)(5,7)(6,8), (1,2)(3,4)]: (8,)}
        """
        return self.element_class(self, PermutationGroup(gap_group=libgap.ConjugateGroup(elm._C, elm._canonicalizing_perm)))

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
        return self._mc, self._dis

    def __hash__(self):
        r"""
        Return the hash of the atomic species.
        """
        return hash(self._element_key())

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

# How to remember the names without ElementCache?
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
            raise ValueError(f"Keys of mapping do not match with domain of {H}")
        if not set(M.values()).issubset(range(1, self._k + 1)):
            raise ValueError(f"Values of mapping must be in the range [1, {self._k}]")
        # each orbit of H must only contain elements from one sort
        for orbit in H.orbits():
            if len(set(M[k] for k in orbit)) > 1:
                raise ValueError(f"For each orbit of {H}, all elements must belong to the same set")
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
        # Trying to avoid libgap.RepresentativeAction; it is slow
        pi = dis_elm._canonicalizing_perm
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
        self._k = indices._k

    # I don't think _element_constructor is ever called, so
    def _element_constructor_(self, x=None):
        print("molecular species element constructor is called!")
        raise NotImplementedError

    @cached_method
    def one(self):
        elm = super().one()
        elm._group = SymmetricGroup(0)
        elm._dompart = dict()
        elm._mc = [0 for _ in range(self._k)]
        elm._tc = 0
        return elm

    def gen(self, x):
        r"""
        Create the molecular species from an atomic species.
        """
        if x not in self._indices:
            raise IndexError(f"{x} is not in the index set")
        at = self._indices(x)
        elm = self._cache_get(self.element_class(self, {at: ZZ.one()}))
        if elm._group is None:
            elm._group = at._dis._C
            elm._dompart = at._dompart
            elm._mc = at._mc
            elm._tc = at._tc
        return elm

    class Element(IndexedFreeAbelianMonoidElement):
        def __init__(self, F, x):
            super().__init__(F, x)
            self._group = None
            self._dompart = None
            self._mc = None
            self._tc = None

        def _assign_group_info(self, other):
            self._group = other._group
            self._dompart = other._dompart
            self._mc = other._mc
            self._tc = other._tc

        def _group_constructor(self):
            r"""
            Construct the group of ``self``.
            """
            temp = self.parent().one()
            for A, p in self._monomial.items():
                at = self.parent().gen(A)
                at_p = at ** p
                temp = temp * at_p
            self._assign_group_info(temp)

        def _elmmul(self, elm1, elm2):
            r"""
            Populate the group info of ``self`` by multiplying
            the groups of ``elm1`` and ``elm2``.
            """
            if elm1._group is None:
                elm1._group_constructor()
            if elm2._group is None:
                elm2._group_constructor()
            if elm1.grade() == 0:
                self._assign_group_info(elm2)
                return
            if elm2.grade() == 0:
                self._assign_group_info(elm1)
                return
            self._mc = tuple(elm1._mc[i] + elm2._mc[i] for i in range(self.parent()._k))
            self._tc = elm1._tc + elm2._tc
            gens1 = elm1._group.gens()
            # Try to avoid gens_small unless necessary
            if len(gens1) > 50:
                gens1 = elm1._group.gens_small()
            gens2 = elm2._group.gens()
            if len(gens2) > 50:
                gens2 = elm2._group.gens_small()
            # always loop over the smaller gens list
            if len(gens1) < len(gens2):
                gens1, gens2 = gens2, gens1
                elm1, elm2 = elm2, elm1
            gens = list(gens1)
            self._dompart = elm1._dompart | {(elm1._tc + k): v for k, v in elm2._dompart.items()}
            for gen in gens2:
                gens.append([tuple(elm1._tc + k for k in cyc) for cyc in gen.cycle_tuples()])
            self._group = PermutationGroup(gens, domain=range(1, elm1._tc + elm2._tc + 1))

        def __floordiv__(self, elt):
            raise NotImplementedError("Cannot cancel in this monoid")

        def _mul_(self, other):
            r"""
            Multiply ``self`` by ``other``.
            """
            res = super()._mul_(other)
            elm = self.parent()._cache_get(res)
            if elm._group is None:
                elm._elmmul(self, other)
            return elm

        def _elmexp(self, other, n):
            r"""
            Populate the group info of ``self`` by exponentiating
            the group of ``other``.
            """
            if other._group is None:
                other._group_constructor()
            temp = self.parent().one()
            while n > 0:
                if n % 2 == 1:
                    temp = temp * other
                other = other * other
                n //= 2
            self._assign_group_info(temp)

        def __pow__(self, n):
            r"""
            Raise ``self`` to the power of ``n``.
            """
            res = super().__pow__(n)
            elm = self.parent()._cache_get(res)
            if elm._group is None:
                elm._elmexp(self, n)
            return elm

        def _element_key(self):
            return self

        def grade(self):
            r"""
            Return the grade of ``self``.
            """
            return self._tc

        def domain(self):
            r"""
            Return the domain of ``self``.
            """
            return FiniteEnumeratedSet(range(1, self._tc + 1))


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
                M = {e: ZZ.one() for e in H.domain()}
            else:
                raise ValueError(f"{x} must be a permutation group")
        else:
            raise ValueError(f"{x} must be a tuple for multivariate species")
        domain_partition = H.disjoint_direct_product_decomposition()
        term = self._indices.one()
        for part in domain_partition:
            term *= self._indices.gen(self._project(H, M, part))
        return self._from_dict({term: ZZ.one()})

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
            r"""
            Return the (partitional) composition of ``self`` with ``args``.
            Uses the word class expansion (Theorem 2.3, Auger paper).
            """
            # should this also have a base ring check or is it coerced?
            if len(args) != self.parent()._k:
                raise ValueError(f"Number of args (= {len(args)}) must equal arity of self (= {len(args)})")
            if any(not isinstance(arg, PolynomialSpecies.Element) for arg in args):
                raise ValueError("All args must be elements of PolynomialSpecies")
            if min(arg.parent()._k for arg in args) != max(arg.parent()._k for arg in args):
                raise ValueError("All args must have same arity")
            if self.is_virtual() or any(arg.is_virtual() for arg in args):
                raise NotImplementedError(f"Only non-virtual species are supported")

            # Now we only have non-virtual species
            res = 0
            for F, coeff in self.monomial_coefficients().items():
                term = 0
                # F(G1;G2;...) each Gi has some molecular decomposition
                # Find the automorphisms
                combos = [Partitions(p, max_length=len(arg)) for arg, p in zip(args, F._mc)]
                # make F such that sort 1 is [123...], and so on
                # cache? lol, cache everything
                dominverted = sorted(F._dompart.keys(), key=lambda x: F._dompart[x])
                pi = libgap.MappingPermListList(dominverted, libgap.eval(f'[1..{F.grade()}]'))
                Fconj = libgap.ConjugateGroup(F._group, pi)
                print('Fconj',Fconj)
                for IVcombo in cartesian_product(combos):
                    # TODO: calculate coefficient (because each group has some multiplicity)
                    # Use F(...) x E(..)^m1 ...
                    print('IVcombo',IVcombo)
                    dsumpartial = list(accumulate(chain.from_iterable(IVcombo),initial=0))
                    print('dsumpartial',dsumpartial)
                    Gtop = libgap.eval(f'SymmetricGroup({F._tc})')
                    Autlist = self._word_automorphism_groups(Fconj, F._mc, IVcombo, Gtop, dsumpartial)
                    print('Autlist',Autlist)
                    Gcombos = cartesian_product([[zip(Gs, permedpart)
                                                  # for each permutation of the partition
                                                  for permedpart in Permutations(IV)
                                                  # we want an embedding of permedpart into the Gs
                                                  for Gs in combinations(arg.support(), r=len(IV))]
                                                  # for each sort
                                                  for arg, IV in zip(args, IVcombo)])
                    for Gcombo in Gcombos:
                        print('Gcombo',Gcombo)
                        gens = []
                        mapping = dict()
                        Fdom = []
                        # Find the gens for each group
                        # Also handles mapping
                        dsum = 0
                        # For the i-th sort
                        for Gs_sort in Gcombo:
                            # For each (group, power) pair
                            for G, cnt in Gs_sort:
                                dlist = G.domain().list()
                                Gdeg = G.grade()
                                curgens = G._group.gens()
                                cur_mapping = G._dompart
                                # calculate the gens
                                for i in range(cnt):
                                    Fdom.append([k + dsum + Gdeg * i for k in sorted(cur_mapping.keys(), key=lambda x: cur_mapping[x])])
                                    images = libgap.MappingPermListList(dlist, [k + dsum + Gdeg * i for k in dlist])
                                    mapping |= {(k + dsum + Gdeg * i): v for k, v in cur_mapping.items()}
                                    # Find the gens for this group
                                    for gen in curgens:
                                        gens.append(gen ** images)
                                dsum += Gdeg * cnt

                        Fdom = libgap(Fdom)
                        Fdomf = libgap.Flat(Fdom)
                        # for each automorphism group
                        for Aut in Autlist:
                            gensrem = []
                            # Find the gens for F
                            for gen in libgap.GeneratorsOfGroup(Aut):
                                # Since we picked F as the stabilizer subgroup of our current class,
                                # we don't have to worry about any extra conditions.
                                perm = libgap.MappingPermListList(Fdomf, libgap.Flat(libgap.Permuted(Fdom, gen)))
                                gensrem.append(perm)
                            totgens = gens + gensrem
                            term += args[0].parent()((PermutationGroup(totgens, domain=range(1, dsum + 1)), mapping))
                res += coeff * term
            return res
        
        def _word_automorphism_groups(self, F, Fmc, IVcombo, Gtop, dsumpartial):
            r"""
            Given a group F and a tuple of partitions on sorts,
            find the word classes and corresponding stabilizer subgroups.
            See Theorem 2.3 of the Auger paper.
            IVcombo is a tuple of integer partitions. 
            (I1[x1, x2,...], ..., In[z1,z2,...])
            """
            # It feels like the output of this function can be HEAVILY optimised.
            # For example, we only care about unordered integer vectors here.
            # Wait, we don't want 0 either. So we actually just want normal integer
            # partitions.

            # create domain
            # can this be taken out? you could cache this too
            domain = list(list(chain.from_iterable(x)) for x in
                          cartesian_product([Permutations(list(chain.from_iterable([i + precard] * v for i, v in enumerate(part, 1))))
                                    for precard, part in zip(accumulate(Fmc, initial=0), IVcombo)]))
            print('domain', domain)

            orig = domain[0]
            orbits = libgap.OrbitsDomain(F, domain, libgap.Permuted)
            grps = []
            for orbit in orbits:
                # Ok, now I have to find a mapping from orbrep to the original
                pi = libgap.RepresentativeAction(Gtop, domain, orbit[0], orig, libgap.Permuted)
                print(pi)
                stab = libgap.Stabilizer(F, domain, orbit[0], libgap.Permuted)
                print('stab',stab)
                print('cj',libgap.ConjugateGroup(stab, pi))
                grps.append(libgap.ConjugateGroup(stab, pi))
            return grps

        def cartesian_product(self, other):
            r"""
            Return the cartesian product of ``self`` with ``other``.
            """
            if not isinstance(other, PolynomialSpecies.Element):
                raise ValueError(f"{other} must be a polynomial species")

            # Do we want to allow cartesian products between different k-variates?
            if self.parent()._k != other.parent()._k:
                return self.parent().zero()

            terms = cartesian_product([self.terms(), other.terms()])
            res = 0
            for t1, t2 in terms:
                H, coeffH = t1.support()[0], t1.coefficients()[0]
                K, coeffK = t2.support()[0], t2.coefficients()[0]
                if H._mc != K._mc:
                    continue
                coeff = coeffH * coeffK
                Sn = SymmetricGroup(H._tc).young_subgroup(H._mc)
                Hgap, Kgap = libgap(H._group), libgap(K._group)
                # We need to normalize H and K to have the same domparts
                Hd = sorted(H._dompart.keys(), key=lambda x: H._dompart[x])
                Kd = sorted(K._dompart.keys(), key=lambda x: K._dompart[x])
                piH = libgap.MappingPermListList(Hd, Kd)
                taus = libgap.DoubleCosetRepsAndSizes(Sn, Hgap, Kgap)
                Hgap = libgap.ConjugateGroup(Hgap, piH)
                for tau, _ in taus:
                    tHt = libgap.ConjugateSubgroup(Hgap, tau)
                    G = PermutationGroup(gap_group=libgap.Intersection(tHt, Kgap), domain=K._dompart.keys())
                    res += coeff * self.parent()((G, K._dompart))
            return res

        def addition_formula(self, arg):
            r"""
            args is a list of the number of terms in each sort.
            Returns the addition formula decomposition of
            H(X1+X2..+Xk, Y1+Y2+..., ...).
            """
            if len(arg) != self.parent()._k:
                raise ValueError(f"Number of args (= {len(arg)}) must equal arity of self (= {len(arg)})")
            if self.is_virtual() or any(x < 0 for x in arg):
                raise NotImplementedError(f"Only non-virtual species are supported")

            # Now we only have non-virtual species
            res = 0
            Parg = PolynomialSpecies(sum(arg))
            for F, coeff in self.monomial_coefficients().items():
                term = 0
                S_top = SymmetricGroup(F._tc).young_subgroup(F._mc)
                dominv = sorted(F._dompart.keys(),key=lambda x: F._dompart[x])
                pi = libgap.MappingPermListList(dominv, libgap.eval(f'[1..{F.grade()}]'))
                Fconj = libgap.ConjugateGroup(F._group, pi)
                for parts in cartesian_product([IntegerVectors(k, l) for k, l in zip(F._mc, arg)]):
                    term2 = 0
                    # parts is a tuple of partitions
                    part = list(chain.from_iterable(parts))
                    S_bottom = SymmetricGroup(F._tc).young_subgroup(part)
                    taus = libgap.DoubleCosetRepsAndSizes(S_top, Fconj, S_bottom)
                    domvals = chain.from_iterable([[i] * v for i, v in enumerate(part, 1)])
                    newdom = {k: v for k, v in zip(range(1, F.grade() + 1), domvals)}
                    for tau, _ in taus:
                        tHt = libgap.ConjugateGroup(Fconj, tau)
                        G = PermutationGroup(gap_group=libgap.Intersection(tHt, S_bottom), domain=F.domain())
                        term2 += Parg((G, newdom))
                    term += term2
                res += coeff * term
            return res