from sage.combinat.integer_vector import IntegerVectors
from sage.misc.cachefunc import cached_method
from sage.monoids.indexed_free_monoid import IndexedFreeAbelianMonoid
from sage.structure.parent import Parent
from sage.structure.element import Element
from sage.structure.element import parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.rings.integer_ring import ZZ
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.categories.sets_with_grading import SetsWithGrading
from sage.categories.graded_algebras_with_basis import GradedAlgebrasWithBasis
from sage.categories.algebras import Algebras
from sage.groups.perm_gps.permgroup import PermutationGroup, PermutationGroup_generic
from sage.groups.perm_gps.permgroup_named import SymmetricGroup
from sage.libs.gap.libgap import libgap
from sage.combinat.free_module import CombinatorialFreeModule
from sage.sets.set import Set

GAP_FAIL = libgap.eval('fail')

def _is_conjugate(G, H1, H2):
    r"""
    Test if ``H1`` and ``H2`` are conjugate subgroups in ``G``.

    EXAMPLES::

        sage: G = SymmetricGroup(3)
        sage: H1 = PermutationGroup([(1,2)])
        sage: H2 = PermutationGroup([(2,3)])
        sage: from sage.rings.burnside import _is_conjugate
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
            else:
                lookup.append(elm)
        else:
            self._cache[key] = [elm]
        return elm

class ConjugacyClassOfSubgroups(Element):
    def __init__(self, parent, C):
        r"""
        Initialize the conjugacy class of ``C``.

        TESTS::

            sage: G = SymmetricGroup(4)
            sage: B = BurnsideRing(G)
            sage: Z4 = CyclicPermutationGroup(4)
            sage: TestSuite(B(Z4)).run()
        """
        Element.__init__(self, parent)
        self._C = C

    def subgroup_of(self):
        r"""
        Return the group which this conjugacy class of subgroups
        belongs to.

        TESTS::

            sage: G = SymmetricGroup(3)
            sage: from sage.rings.burnside import ConjugacyClassesOfSubgroups
            sage: C = ConjugacyClassesOfSubgroups(G)
            sage: Z3 = CyclicPermutationGroup(3)
            sage: CZ3 = C(Z3)
            sage: CZ3.subgroup_of()
            Symmetric group of order 3! as a permutation group
        """
        return self.parent()._G

    def __hash__(self):
        r"""
        Return the hash of the representative of the conjugacy class.

        TESTS::

            sage: G = SymmetricGroup(3)
            sage: B = BurnsideRing(G)
            sage: H1 = B(PermutationGroup([(1,2)]))
            sage: H2 = B(PermutationGroup([(2,3)]))
            sage: hash(H1) == hash(H2)
            True
        """
        return hash(self._C)

    @cached_method
    def _element_key(self):
        r"""
        Return tuple of group invariants associated with H.
        """
        # should this be a lazy attribute instead?
        return tuple([self._C.order(), self._C.degree()])
    
    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        TESTS::

            sage: G = SymmetricGroup(3)
            sage: from sage.rings.burnside import ConjugacyClassesOfSubgroups
            sage: C = ConjugacyClassesOfSubgroups(G)
            sage: Z3 = CyclicPermutationGroup(3)
            sage: CZ3 = C(Z3)
            sage: CZ3
            [(1,2,3)]
        """
        return repr(self._C.gens_small())

    def __le__(self, other):
        r"""
        Return if this element is less than or equal to ``other``.

        ``self`` is less or equal to ``other`` if it is conjugate to
        a subgroup of ``other`` in the parent group.

        EXAMPLES:

        We recreate the poset of conjugacy classes of subgroups of
        `S_4`::

            sage: G = SymmetricGroup(4)
            sage: B = BurnsideRing(G)
            sage: P = Poset([B.basis(), lambda b, c: b <= c])
            sage: len(P.cover_relations())
            17
        """
        return (isinstance(other, ConjugacyClassOfSubgroups)
                and (GAP_FAIL != libgap.ContainedConjugates(self.subgroup_of(),
                                                            other._C,
                                                            self._C,
                                                            True)))

    def __eq__(self, other):
        r"""
        Return if this element is equal to ``other``.

        Two elements compare equal if they are conjugate subgroups in the parent group.

        TESTS::

            sage: G = SymmetricGroup(4)
            sage: B = BurnsideRing(G)
            sage: H1 = PermutationGroup([(1,2)])
            sage: H2 = PermutationGroup([(2,3)])
            sage: B[H1] == B[H2]
            True
        """
        return (isinstance(other, ConjugacyClassOfSubgroups)
                and _is_conjugate(self.subgroup_of(), self._C, other._C))

    def __lt__(self, other):
        r"""
        Return if this element is less than ``other``.

        TESTS::

            sage: G = SymmetricGroup(4)
            sage: B = BurnsideRing(G)
            sage: H1 = PermutationGroup([(1,2)])
            sage: H2 = PermutationGroup([(2,3)])
            sage: H3 = PermutationGroup([(1,2),(3,4)])
            sage: B[H1] < B[H2]
            False
            sage: B[H1] == B[H2]
            True
            sage: B[H1] < B[H3]
            True
            sage: B[H2] < B[H3]
            True
        """
        return self <= other and self != other

class ConjugacyClassesOfSubgroups(Parent, ElementCache):
    def __init__(self, G):
        r"""
        Initialize the set of conjugacy classes of ``G``.

        INPUT:

        ``G`` -- a group.

        TESTS::

            sage: from sage.rings.burnside import ConjugacyClassesOfSubgroups
            sage: G = CyclicPermutationGroup(4)
            sage: C = ConjugacyClassesOfSubgroups(G)
            sage: TestSuite(C).run()
        """
        self._G = G
        Parent.__init__(self, category=FiniteEnumeratedSets())
        ElementCache.__init__(self)

    def __eq__(self, other):
        r"""
        Return if ``self`` is equal to ``other``.

        TESTS::

            sage: from sage.rings.burnside import ConjugacyClassesOfSubgroups
            sage: Z4 = CyclicPermutationGroup(4)
            sage: D4 = DiCyclicGroup(4)
            sage: C1 = ConjugacyClassesOfSubgroups(Z4)
            sage: C2 = ConjugacyClassesOfSubgroups(D4)
            sage: C1 == C2
            False
        """
        if not isinstance(other, ConjugacyClassesOfSubgroups):
            return False
        return self._G == other._G

    def __hash__(self):
        r"""
        Return the hash of the underlying group.

        TESTS::

            sage: from sage.rings.burnside import ConjugacyClassesOfSubgroups
            sage: Z4 = CyclicPermutationGroup(4)
            sage: D4 = DiCyclicGroup(4)
            sage: C1 = ConjugacyClassesOfSubgroups(Z4)
            sage: C2 = ConjugacyClassesOfSubgroups(D4)
            sage: hash(C1) == hash(C2)
            False
        """
        return hash(self._G)

    def _element_constructor_(self, x):
        r"""
        Construct the conjugacy class of subgroups containing ``x``.

        TESTS::

            sage: from sage.rings.burnside import ConjugacyClassesOfSubgroups
            sage: G = SymmetricGroup(4)
            sage: C = ConjugacyClassesOfSubgroups(G)
            sage: Z4 = CyclicPermutationGroup(4)
            sage: C(Z4)
            [(1,2,3,4)]
        """
        if parent(x) == self:
            return x
        if isinstance(x, PermutationGroup_generic):
            if x.is_subgroup(self._G):
                elm = self.element_class(self, self._G.subgroup(x.gens()))
                return self._cache_get(elm)
            raise ValueError(f"unable to convert {x} into {self}: not a subgroup of {self._G}")
        raise ValueError("unable to convert {x} into {self}")


    def __iter__(self):
        r"""
        Return iterator over conjugacy classes of subgroups of the group.

        TESTS::

            sage: G = SymmetricGroup(3)
            sage: B = BurnsideRing(G)
            sage: [g for g in B._indices]
            [[()], [(1,2)], [(1,2,3)], 1]
        """
        return iter(self(H) for H in self._G.conjugacy_classes_subgroups())

    def __contains__(self, H):
        r"""
        Return if ``H`` is a subgroup of the group.

        TESTS::

            sage: G = SymmetricGroup(4)
            sage: B = BurnsideRing(G)
            sage: Z4 = CyclicPermutationGroup(4)
            sage: Z4 in B._indices
            True
        """
        if parent(H) == self:
            return True
        return (isinstance(H, PermutationGroup_generic)
                and H.is_subgroup(self._G))

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        TESTS::

            sage: G = SymmetricGroup(4)
            sage: B = BurnsideRing(G)
            sage: B._indices
            Conjugacy classes of subgroups of Symmetric group of order 4! as a permutation group
        """
        return "Conjugacy classes of subgroups of " + repr(self._G)

    Element = ConjugacyClassOfSubgroups

class BurnsideRing(CombinatorialFreeModule):
    def __init__(self, G, base_ring=ZZ):
        r"""
        The Burnside ring of the group ``G``.

        INPUT:

        ``G`` -- a group.
        ``base_ring`` -- the ring of coefficients. Default value is ``ZZ``.

        EXAMPLES::

            sage: G = SymmetricGroup(6)
            sage: B = BurnsideRing(G)
            sage: X = Subsets(6, 2)
            sage: b = B.construct_from_action(lambda g, x: X([g(e) for e in x]), X)
            sage: B.basis().keys()._cache
            {(48,
            6): [(Subgroup generated by [(3,5,4,6), (1,2)(4,5,6)] of (Symmetric group of order 6! as a permutation group),
            [(1,2)(4,5,6), (3,5,4,6)])],
            (720,
            6): [(Subgroup generated by [(1,2,3,4,5,6), (1,2)] of (Symmetric group of order 6! as a permutation group),
            1)]}
            sage: b^2
            B[[(4,6,5), (4,6)]] + B[[(5,6), (3,4)(5,6), (1,2)(5,6)]] + B[[(1,2)(4,5,6), (3,5,4,6)]]
            sage: B.basis().keys()._cache
            {(6,
            6): [(Subgroup generated by [(4,6), (4,6,5)] of (Symmetric group of order 6! as a permutation group),
            [(4,6,5), (4,6)])],
            (8,
            6): [(Subgroup generated by [(5,6), (1,2)(5,6), (3,4)(5,6)] of (Symmetric group of order 6! as a permutation group),
            [(5,6), (3,4)(5,6), (1,2)(5,6)])],
            (48,
            6): [(Subgroup generated by [(3,5,4,6), (1,2)(4,5,6)] of (Symmetric group of order 6! as a permutation group),
            [(1,2)(4,5,6), (3,5,4,6)])],
            (720,
            6): [(Subgroup generated by [(1,2,3,4,5,6), (1,2)] of (Symmetric group of order 6! as a permutation group),
            1)]}

            sage: G = SymmetricGroup(4)
            sage: B = BurnsideRing(G)
            sage: X = Subsets(4, 2)
            sage: b = B.construct_from_action(lambda g, x: X([g(e) for e in x]), X)
            sage: b.tensor(b)
            B[[(3,4), (1,2)(3,4)]] # B[[(3,4), (1,2)(3,4)]]
            sage: (b.tensor(b))^2
            B[[()]] # B[[()]] + 2*B[[()]] # B[[(3,4), (1,2)(3,4)]] + 2*B[[(3,4), (1,2)(3,4)]] # B[[()]] + 4*B[[(3,4), (1,2)(3,4)]] # B[[(3,4), (1,2)(3,4)]]

        TESTS::

            sage: G = SymmetricGroup(4)
            sage: B = BurnsideRing(G)
            sage: TestSuite(B).run()
        """
        self._G = G
        basis_keys = ConjugacyClassesOfSubgroups(G)
        category = Algebras(base_ring).Commutative().WithBasis()
        CombinatorialFreeModule.__init__(self, base_ring, basis_keys,
                                        category=category, prefix="B")
        self._indices(G).rename("1")

    def __getitem__(self, H):
        r"""
        Return the basis element indexed by ``H``.

        ``H`` must be a subgroup of the group.

        EXAMPLES::

            sage: G = SymmetricGroup(4)
            sage: B = BurnsideRing(G)
            sage: Z4 = CyclicPermutationGroup(4)
            sage: B[Z4]
            B[[(1,2,3,4)]]
        """
        return self._from_dict({self._indices(H): 1})

    def construct_from_action(self, action, domain):
        r"""
        Construct an element of the Burnside ring from a group action.

        INPUT:

        - ``action`` - an action on ``domain``
        - ``domain`` - a finite set

        EXAMPLES::

            sage: G = SymmetricGroup(4)
            sage: B = BurnsideRing(G)

        We create a group action of `S_4` on two-element subsets::

            sage: X = Subsets(4, 2)
            sage: a = lambda g, x: X([g(e) for e in x])
            sage: B.construct_from_action(a, X)
            B[[(3,4), (1,2)(3,4)]]

        Next, we create a group action of `S_4` on itself via conjugation::

            sage: X = G
            sage: a = lambda g, x: g*x*g.inverse()
            sage: B.construct_from_action(a, X)
            B[[(3,4), (1,3)(2,4), (1,4)(2,3)]] + B[[(2,4,3)]] + B[[(3,4), (1,2)(3,4)]] + B[[(1,2,3,4)]] + B[1]

        TESTS::

            sage: G = SymmetricGroup(4)
            sage: B = BurnsideRing(G)
            sage: [b.support()[0]._C.order() for b in B.gens()]
            [1, 2, 2, 3, 4, 4, 4, 6, 8, 12, 24]
            sage: sorted((o, len(l)) for o, l in B._indices._cache.items())
            [((1, 4), 1),
            ((2, 4), 2),
            ((3, 4), 1),
            ((4, 4), 3),
            ((6, 4), 1),
            ((8, 4), 1),
            ((12, 4), 1),
            ((24, 4), 1)]

            sage: G = SymmetricGroup(4)
            sage: B = BurnsideRing(G)
            sage: B(-3)
            -3*B[1]
        """
        def find_stabilizer(action, pnt):
            stabilizer = []
            for g in self._G:
                if action(g, pnt) == pnt:
                    stabilizer.append(g)
            H = self._G.subgroup(stabilizer)
            gens = H.gens_small()
            return self._G.subgroup(gens)

        H = PermutationGroup(self._G.gens(), action=action, domain=domain)
        # decompose H into orbits
        orbit_list = H.orbits()
        # find the stabilizer subgroups
        stabilizer_list = [find_stabilizer(action, orbit[0]) for orbit in orbit_list]
        # normalize each summand and collect terms
        from collections import Counter
        C = Counter([self._indices(stabilizer) for stabilizer in stabilizer_list])
        return self._from_dict(dict(C))

    @cached_method
    def one_basis(self):
        r"""
        Returns the underlying group, which indexes the one of this algebra,
        as per :meth:`AlgebrasWithBasis.ParentMethods.one_basis`.

        EXAMPLES::

            sage: G = DiCyclicGroup(4)
            sage: B = BurnsideRing(G)
            sage: B.one_basis()
            1
        """
        return self._indices(self._G)

    def product_on_basis(self, H, K):
        r"""
        Return the product of the basis elements indexed by ``H`` and ``K``.

        For the symmetric group, this is also known as the Hadamard
        or tensor product of group actions.

        EXAMPLES::

            sage: G = SymmetricGroup(3)
            sage: B = BurnsideRing(G)
            sage: matrix([[b * c for b in B.gens()] for c in B.gens()])
            [           6*B[[()]]            3*B[[()]]            2*B[[()]]              B[[()]]]
            [           3*B[[()]] B[[()]] + B[[(1,2)]]              B[[()]]           B[[(1,2)]]]
            [           2*B[[()]]              B[[()]]       2*B[[(1,2,3)]]         B[[(1,2,3)]]]
            [             B[[()]]           B[[(1,2)]]         B[[(1,2,3)]]                 B[1]]

        TESTS::

            sage: G = SymmetricGroup(3)
            sage: B = BurnsideRing(G)
            sage: Z3 = CyclicPermutationGroup(3)
            sage: Z2 = CyclicPermutationGroup(2)
            sage: from sage.rings.burnside import ConjugacyClassesOfSubgroups
            sage: C = ConjugacyClassesOfSubgroups(G)
            sage: B.product_on_basis(C(Z2), C(Z3))
            B[[()]]
        """
        g_reps = [rep for rep, size in libgap.DoubleCosetRepsAndSizes(self._G, H._C, K._C)]
        from collections import Counter
        C = Counter()
        for g in g_reps:
            g_sup_K = libgap.ConjugateSubgroup(K._C, g)
            P = self._G.subgroup(gap_group=libgap.Intersection(H._C, g_sup_K))
            C[self._indices(P)] += 1
        return self._from_dict(dict(C))

    def group(self):
        r"""
        Return the underlying group.

        EXAMPLES::

            sage: G = DiCyclicGroup(4)
            sage: B = BurnsideRing(G)
            sage: B.group()
            Dicyclic group of order 16 as a permutation group
        """
        return self._G

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: G = SymmetricGroup(4)
            sage: B = BurnsideRing(G)
            sage: B
            Burnside ring of Symmetric group of order 4! as a permutation group
        """
        return "Burnside ring of " + repr(self._G)

class ConjugacyClassOfDirectlyIndecomposableSubgroups(ConjugacyClassOfSubgroups):
    def __init__(self, parent, C):
        r"""
        A conjugacy class of directly indecomposable subgroups.
        """
        if len(C.disjoint_direct_product_decomposition()) != 1:
            raise ValueError(f"{C} is not directly indecomposable")
        ConjugacyClassOfSubgroups.__init__(self, parent, C)

    @cached_method
    def subgroup_of(self):
        r"""
        Return the group which this conjugacy class of
        directly indecomposable subgroups belongs to.
        """
        return SymmetricGroup(self._C.degree())
    
    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: At = UnivariateAtomicConjugacyClasses()
            sage: g = At(CyclicPermutationGroup(3)); g
            {3, [(1,2,3)]}
        """
        return "{" + f"{self._C.degree()}, {self._C.gens_small()}" + "}"

class AtomicSpecies(ConjugacyClassOfDirectlyIndecomposableSubgroups):
    def __init__(self, parent, C, multicardinality):
        r"""
        An atomic species centered on a given multicardinality: a directly
        indecomposable subgroup paired with an integer vector.
        """
        total_cardinality = sum(multicardinality)
        if C.degree() != sum(multicardinality):
            raise ValueError(f"Degree of {C} (= {C.degree()}) must equal the total cardinality {total_cardinality}")
        G = SymmetricGroup(total_cardinality).young_subgroup(multicardinality)
        if not C.is_subgroup(G):
            raise ValueError(f"{C} is not a subgroup of {G}")
        ConjugacyClassOfDirectlyIndecomposableSubgroups.__init__(self, parent, C)
        self._mc = multicardinality

    @cached_method
    def subgroup_of(self):
        r"""
        Return the group which the underlying permutation group
        of this atomic species belongs to.
        """
        return SymmetricGroup(sum(self._mc)).young_subgroup(self._mc)

    def __hash__(self):
        r"""
        Return the hash of the atomic species.
        """
        return hash(tuple([self._C, self._mc]))

    @cached_method
    def _element_key(self):
        r"""
        Return a lookup key for ``self``.
        """
        return tuple([*super()._element_key(), self._mc])

    def __eq__(self, other):
        r"""
        Two atomic species are equal if the underlying groups are conjugate,
        and their multicardinalities are equal.
        """
        if parent(self) != parent(other):
            return False
        return self._mc == other._mc and super().__eq__(other)

    def _repr_(self):
        r"""
        Return a string representation of ``self``.
        """
        return "{" + f"{self._C.gens_small()}: {self._mc}" + "}"

class MultivariateAtomicSpecies(UniqueRepresentation, Parent, ElementCache):
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

    def __hash__(self):
        r"""
        Return a hash for ``self``.
        """
        return hash(self._k)

    def __eq__(self, other):
        r"""
        Needed for unique representation behaviour.

        TESTS::

            sage: At1 = MultivariateAtomicSpecies(1); At1
            sage: At2 = MultivariateAtomicSpecies(2); At2
            sage: At1_2 = MultivariateAtomicSpecies(1); At1_2
            sage: At1 is At1_2
            True
            sage: At1 is At2
            False
            sage: At1_2 is At2
            False
        """
        if type(self) != type(other):
            return False
        return self._k == other._k
    
    @cached_method
    def an_element(self):
        """
        Return an element of ``self``.
        """
        return self._element_constructor_((SymmetricGroup(self._k).young_subgroup([1] * self._k),
                                           [1] * self._k))

    def _normalize(self, H, f):
        r"""
        Normalize `H` and return `H` and the multicardinality.
        """
        # TODO: Complete the documentation.
        L = [[] for _ in range(self._k)]
        for k, v in f.items():
            L[v - 1].append(k)
        Lc = sum(L, [])
        Ls = [len(l) for l in L]
        mapping = {v: i for i, v in enumerate(Lc, 1)}
        normalized_gens = [[tuple(mapping[x] for x in cyc) for cyc in gen.cycle_tuples()] for gen in H.gens_small()]
        return PermutationGroup(gens=normalized_gens), self._grading_set(Ls)

    def _element_constructor_(self, x):
        r"""
        Construct the `k`-variate molecular species with the given data.

        INPUT:

        - ``x`` - an element of ``self`` or a tuple ``(H, f)`` where `H` is
        the directly indecomposable permutation group representation for the
        `k`-variate atomic species and `f` is a ``dict`` mapping each element
        of the domain of `H` to integers in `\{ 1 \ldots k \}`, representing
        the set to which the element belongs.
        """
        # BUG: SymmetricGroup(0) fails to be constructed.
        if parent(x) == self:
            return x
        H, f = x
        if Set(H.domain()) != Set(f.keys()):
            raise ValueError(f"Keys of {f} do not match with domain of {H}")
        if not Set(f.values()).issubset(Set(range(1, self._k + 1))):
            raise ValueError(f"Values of {f} must be in the range [1, {self._k}]")
        if isinstance(H, PermutationGroup_generic) and isinstance(f, dict):
            elm = self.element_class(self, *self._normalize(H, f))
            return self._cache_get(elm)
        raise ValueError("unable to convert {x} into {self}")
    
    def __getitem__(self, x):
        r"""
        Call ``_element_constructor_`` on ``x``.
        """
        return self._element_constructor_(x)
    
    def _repr_(self):
        r"""
        Return a string representation of ``self``.
        """
        return f"Infinite set of {self._k}-variate atomic species"
    
    Element = AtomicSpecies

class PolynomialMolecularDecomposition(CombinatorialFreeModule):
    def __init__(self, k, base_ring=ZZ):
        r"""
        Ring of `k`-variate virtual species.
        """
        # should we pass a category to basis_keys?
        basis_keys = IndexedFreeAbelianMonoid(MultivariateAtomicSpecies(k),
                                              prefix='', bracket=False)
        category = GradedAlgebrasWithBasis(base_ring).Commutative()
        CombinatorialFreeModule.__init__(self, base_ring,
                                        basis_keys=basis_keys,
                                        category=category,
                                        prefix='', bracket=False)
        self._k = k
        self._atomic_basis = basis_keys.indices()

    def _project(self, H, f, part):
        r"""
        Project `H` onto a subset ``part`` of its domain.
        ``part`` must be a union of cycles, but this is not checked.
        """
        restricted_gens = [[cyc for cyc in gen.cycle_tuples() if cyc[0] in part] for gen in H.gens_small()]
        mapping = {p: f[p] for p in part}
        return tuple([PermutationGroup(gens=restricted_gens), mapping])
    
    def _element_constructor_(self, x):
        r"""
        Construct the `k`-variate molecular species with the given data.

        INPUT:

        - ``x`` - an element of ``self`` or a tuple ``(H, f)`` where `H` is
        the permutation group representation for the species and `f` is a
        ``dict`` mapping each element of the domain of `H` to integers in
        `\{ 1 \ldots k \}`, representing the set to which the element belongs.
        """
        if parent(x) == self:
            return x
        H, f = x
        if Set(H.domain()) != Set(f.keys()):
            raise ValueError(f"Keys of {f} do not match with domain of {H}")
        if not Set(f.values()).issubset(Set(range(1, self._k + 1))):
            raise ValueError(f"Values of {f} must be in the range [1, {self._k}]")
        if isinstance(H, PermutationGroup_generic) and isinstance(f, dict):
            domain_partition = H.disjoint_direct_product_decomposition()
            from collections import Counter
            C = Counter([self._atomic_basis(self._project(H, f, part)) for part in domain_partition])
            return self._from_dict({self._indices(dict(C)): 1})
        raise ValueError("unable to convert {x} into {self}")

    def __getitem__(self, x):
        r"""
        Calls _element_constructor_ on x.
        """
        return self._element_constructor_(x)

    @cached_method
    def one_basis(self):
        r"""
        Returns SymmetricGroup(0), which indexes the one of this algebra,
        as per :meth:`AlgebrasWithBasis.ParentMethods.one_basis`.

        EXAMPLES::

            sage: P = PolynomialMolecularDecomposition()
            sage: P.one_basis()
            1
        """
        return self._element_constructor_(tuple([SymmetricGroup(0), dict()]))

    def product_on_basis(self, H, K):
        r"""
        Return the product of the basis elements indexed by `H` and `K`.

        EXAMPLES::

            sage: Mol = UnivariateMolecularConjugacyClasses()
            sage: P = PolynomialMolecularDecomposition()
            sage: L1 = SymmetricGroup(3).conjugacy_classes_subgroups()
            sage: L2 = SymmetricGroup(2).conjugacy_classes_subgroups()
            sage: matrix([[P(x) * P(y) for x in L1] for y in L2])
            [                       {1, [()]}^5           {1, [()]}^3*{2, [(1,2)]}         {3, [(1,2,3)]}*{1, [()]}^2  {3, [(1,2,3), (2,3)]}*{1, [()]}^2]
            [          {1, [()]}^3*{2, [(1,2)]}           {1, [()]}*{2, [(1,2)]}^2        {3, [(1,2,3)]}*{2, [(1,2)]} {3, [(1,2,3), (2,3)]}*{2, [(1,2)]}]
        """
        return self._from_dict({H * K: 1})

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: P = PolynomialMolecularDecomposition(2)
            sage: P
            Polynomial Molecular Decomposition
        """
        return f"Ring of {self._k}-variate virtual species"
