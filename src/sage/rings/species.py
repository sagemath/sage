from sage.combinat.integer_vector import IntegerVectors
from sage.misc.cachefunc import cached_method
from sage.monoids.indexed_free_monoid import IndexedFreeAbelianMonoid
from sage.structure.parent import Parent
from sage.structure.element import parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.rings.integer_ring import ZZ
from sage.categories.sets_with_grading import SetsWithGrading
from sage.categories.graded_algebras_with_basis import GradedAlgebrasWithBasis
from sage.groups.perm_gps.permgroup import PermutationGroup, PermutationGroup_generic
from sage.groups.perm_gps.permgroup_named import SymmetricGroup
from sage.libs.gap.libgap import libgap
from sage.combinat.free_module import CombinatorialFreeModule
from sage.sets.set import Set
from sage.rings.burnside import ConjugacyClassOfSubgroups, ElementCache, _is_conjugate

GAP_FAIL = libgap.eval('fail')

class ConjugacyClassOfDirectlyIndecomposableSubgroups(ConjugacyClassOfSubgroups):
    def __init__(self, parent, C):
        r"""
        A conjugacy class of directly indecomposable subgroups.
        """
        if len(C.disjoint_direct_product_decomposition()) > 1:
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
        """
        return "{" + f"{self._C.degree()}, {self._C.gens_small()}" + "}"

    def __le__(self, other):
        r"""
        Return if this element is less than or equal to ``other``.

        ``self`` is less or equal to ``other`` if it is conjugate to
        a subgroup of ``other`` in the parent group of larger degree.
        """
        # If selfdeg > otherdeg, return False
        # then selfdeg <= otherdeg
        return (isinstance(other, ConjugacyClassOfSubgroups)
                and self._C.degree() <= other._C.degree()
                and (self._C.degree() < other._C.degree() or
                     (GAP_FAIL != libgap.ContainedConjugates(self.subgroup_of(),
                                                             other._C, self._C, True))))

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
                and self._C.degree() == other._C.degree() and
                _is_conjugate(self.subgroup_of(), self._C, other._C))

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
        self._tc = total_cardinality

    @cached_method
    def subgroup_of(self):
        r"""
        Return the group which the underlying permutation group
        of this atomic species belongs to.
        """
        return SymmetricGroup(self._tc).young_subgroup(self._mc)

    def __hash__(self):
        r"""
        Return the hash of the atomic species.
        """
        return hash(tuple([self._C, self._mc, self._tc]))

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
                                           {e: i for i, e in enumerate(range(1, self._k + 1), 1)}))

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
        P = PermutationGroup(gens=normalized_gens)
        # Fix for SymmetricGroup(0)
        if H.degree() == 0:
            P = SymmetricGroup(0)
        return P, self._grading_set(Ls)

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
        if parent(x) == self:
            return x
        H, f = x
        if Set(H.domain()) != Set(f.keys()):
            raise ValueError(f"Keys of {f} do not match with domain of {H} (= {H.domain()})")
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

    def __contains__(self, x):
        r"""
        Return if ``x`` is in ``self``.
        """
        if parent(x) == self:
            return True
        try:
            self._element_constructor_(x)
        except:
            return False
        return True
    
    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        TESTS::

            sage: At1 = MultivariateAtomicSpecies(1)
            sage: At1
            Infinite set of 1-variate atomic species
            sage: At2 = MultivariateAtomicSpecies(2)
            sage: At2
            Infinite set of 2-variate atomic species
        """
        return f"Infinite set of {self._k}-variate atomic species"
    
    Element = AtomicSpecies

class PolynomialMolecularDecomposition(CombinatorialFreeModule):
    def __init__(self, k, base_ring=ZZ):
        r"""
        Ring of `k`-variate virtual species.

        TESTS::

            sage: P = PolynomialMolecularDecomposition(1)
            sage: TestSuite(P).run()
            sage: P2 = PolynomialMolecularDecomposition(2)
            sage: TestSuite(P2).run()
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
        return tuple([PermutationGroup(gens=restricted_gens, domain=part), mapping])
    
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

        TESTS::

            sage: P = PolynomialMolecularDecomposition(1)
            sage: At1 = MultivariateAtomicSpecies(1)
            sage: At1((SymmetricGroup(1), {1: 1})).rename("X")
            sage: X = (SymmetricGroup(2).young_subgroup([1, 1]), {1: 1, 2: 1})
            sage: P[X]
            X^2
            sage: P2=PolynomialMolecularDecomposition(2)
            sage: At2=MultivariateAtomicSpecies(2)
            sage: At2((SymmetricGroup(1), {1: 1})).rename("X")
            sage: At2((SymmetricGroup(1), {1: 2})).rename("Y")
            sage: XY = (SymmetricGroup(2).young_subgroup([1, 1]), {1: 1, 2: 2})
            sage: P2[XY]
            Y*X
        """
        return self._element_constructor_(x)

    @cached_method
    def one_basis(self):
        r"""
        Returns SymmetricGroup(0), which indexes the one of this algebra,
        as per :meth:`AlgebrasWithBasis.ParentMethods.one_basis`.

        EXAMPLES::

            sage: P = PolynomialMolecularDecomposition(1)
            sage: P.one_basis()
            {[]: [0]}
            sage: P2 = PolynomialMolecularDecomposition(2)
            sage: P2.one_basis()
            {[]: [0, 0]}
        """
        return self._indices({self._atomic_basis(tuple([SymmetricGroup(0), dict()])): 1})

    @cached_method
    def an_element(self):
        """
        Return an element of ``self``.

            sage: P=PolynomialMolecularDecomposition(1)
            sage: P.an_element()
            {[]: [0]}
            sage: P2=PolynomialMolecularDecomposition(2)
            sage: P2.an_element()
            {[]: [0, 0]}
        """
        return self.one()

    def product_on_basis(self, H, K):
        r"""
        Return the product of the basis elements indexed by `H` and `K`.

        EXAMPLES::

            sage: P=PolynomialMolecularDecomposition(1)
            sage: d2 = {e: 1 for e in range(1, 3)}
            sage: d3 = {e: 1 for e in range(1, 4)}
            sage: L1 = [(H, d3) for H in SymmetricGroup(3).conjugacy_classes_subgroups()]
            sage: L2 = [(H, d2) for H in SymmetricGroup(2).conjugacy_classes_subgroups()]
            sage: matrix([[P(x) * P(y) for x in L1] for y in L2])
            [                         {[()]: [1]}^5           {[()]: [1]}^3*{[(1,2)]: [2]}         {[(1,2,3)]: [3]}*{[()]: [1]}^2  {[(1,2,3), (2,3)]: [3]}*{[()]: [1]}^2]
            [          {[()]: [1]}^3*{[(1,2)]: [2]}           {[()]: [1]}*{[(1,2)]: [2]}^2        {[(1,2,3)]: [3]}*{[(1,2)]: [2]} {[(1,2,3), (2,3)]: [3]}*{[(1,2)]: [2]}]
        """
        # Hacky workaround for handling the one of the algebra :(
        if H == self.one_basis() and K == self.one_basis():
            return self._from_dict({self.one_basis(): 1})
        elif H == self.one_basis():
            return self._from_dict({K: 1})
        elif K == self.one_basis():
            return self._from_dict({H: 1})
        return self._from_dict({H * K: 1})

    def degree_on_basis(self, m):
        r"""
        Return the degree of the basis element indexed by ``m``
        in ``self``.
        """
        d = m.dict()
        return sum(at._tc * p for at, p in d.items())

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: P = PolynomialMolecularDecomposition(1)
            sage: P
            Ring of 1-variate virtual species
            sage: P2 = PolynomialMolecularDecomposition(2)
            sage: P2
            Ring of 2-variate virtual species
        """
        return f"Ring of {self._k}-variate virtual species"
