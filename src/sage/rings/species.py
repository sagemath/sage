from itertools import chain
from sage.categories.graded_algebras_with_basis import GradedAlgebrasWithBasis
from sage.categories.infinite_enumerated_sets import InfiniteEnumeratedSets
from sage.categories.sets_with_grading import SetsWithGrading
from sage.combinat.free_module import CombinatorialFreeModule
from sage.combinat.integer_vector import IntegerVectors
from sage.groups.perm_gps.constructor import PermutationGroupElement
from sage.groups.perm_gps.permgroup import PermutationGroup, PermutationGroup_generic
from sage.groups.perm_gps.permgroup_named import SymmetricGroup
from sage.libs.gap.libgap import libgap
from sage.misc.cachefunc import cached_method
from sage.misc.lazy_attribute import lazy_attribute
from sage.monoids.indexed_free_monoid import IndexedFreeAbelianMonoid, IndexedFreeAbelianMonoidElement
from sage.rings.integer_ring import ZZ
from sage.structure.element import Element, parent
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation

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
            elm = self._canonical_label(elm)
            lookup.append(elm)
        else:
            elm = self._canonical_label(elm)
            self._cache[key] = [elm]
        return elm

    def _canonical_label(self, elm):
        return elm


class ConjugacyClassOfDirectlyIndecomposableSubgroups(Element):
    def __init__(self, parent, C):
        r"""
        A conjugacy class of directly indecomposable subgroups.
        """
        Element.__init__(self, parent)
        self._C = C
        self._smallgen = C.gens_small()

    @cached_method
    def _element_key(self):
        r"""
        Return the key for this element.
        """
        return self._C.degree(), self._C.order()

    def __hash__(self):
        return hash(self._element_key)

    def _repr_(self):
        r"""
        Return a string representation of ``self``.
        """
        return f"{self._smallgen}"

    def __eq__(self, other):
        r"""
        Return whether ``self`` is equal to ``other``.

        ``self`` is equal to ``other`` if they have the same degree (say `n`)
        and order and are conjugate within `S_n`.
        """
        return (isinstance(other, ConjugacyClassOfDirectlyIndecomposableSubgroups)
                and self._element_key() == other._element_key()
                and _is_conjugate(SymmetricGroup(self._C.degree()),
                                  self._C, other._C))


class ConjugacyClassesOfDirectlyIndecomposableSubgroups(UniqueRepresentation,
                                                        Parent, ElementCache):
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

    def _canonical_label(self, elm):
        """

        EXAMPLES::

            sage: G = PermutationGroup([[(1,3),(4,7)], [(2,5),(6,8)], [(1,4),(2,5),(3,7)]])
            sage: A = AtomicSpecies(1)
            sage: A(G)
            {[(5,6)(7,8), (1,2)(5,7)(6,8), (1,2)(3,4)]: (8,)}
        """
        orbits = sorted(elm._C.orbits(), key=len)
        pi = PermutationGroupElement([e for o in orbits for e in o])
        pi_inv = pi.inverse()
        gens = [pi * g * pi_inv for g in elm._smallgen]
        return self.element_class(self, PermutationGroup(gens))

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

    @cached_method
    def _element_key(self):
        r"""
        Return a lookup key for ``self``.
        """
        return self._dis._C, self._mc

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
        # Here it is enough to compare mc
        return self._mc == other._mc and self._dis == self._dis

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
        mapping = {v: i for i, v in enumerate(H.domain(), 1)}
        normalized_gens = [[tuple(mapping[x] for x in cyc)
                            for cyc in gen.cycle_tuples()]
                           for gen in H.gens_small()]
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
        perm = libgap.RepresentativeAction(SymmetricGroup(H_norm.degree()), H_norm, dis_elm._C)
        dompart_norm = {ZZ(k ** perm): v for k, v in dompart.items()}
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
        # This needs to be improved.
        if parent(x) == self:
            return True
        H, M = x
        try:
            H_norm, _ = self._normalize(H, M)
            self._dis_ctor(H_norm)
        except ValueError:
            return False
        return True

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

class MolecularSpecies(IndexedFreeAbelianMonoid):
    class Element(IndexedFreeAbelianMonoidElement):
        @lazy_attribute
        def _group(self):
            def gmul(H, K):
                n, m = H.degree(), K.degree()
                gens = H.gens_small()
                for gen in K.gens_small():
                    shift = libgap.MappingPermListList(K.domain().list(), [n + k for k in K.domain()])
                    gens.append(shift ** -1 * gen * shift)
                return PermutationGroup(gens, domain=range(1, n + m + 1))
            G = SymmetricGroup(0)
            for A, p in self._monomial.items():
                for _ in range(p):
                    G = gmul(G, A._dis._C)
            return G

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

    def _project(self, H, f, part):
        r"""
        Project `H` onto a subset ``part`` of its domain.

        ``part`` must be a union of cycles, but this is not checked.
        """
        restricted_gens = [[cyc for cyc in gen.cycle_tuples() if cyc[0] in part] for gen in H.gens_small()]
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

            sage: P = PolynomialSpecies(1)
            sage: P
            Ring of 1-variate virtual species
            sage: P2 = PolynomialSpecies(2)
            sage: P2
            Ring of 2-variate virtual species
        """
        return f"Ring of {self._k}-variate virtual species"

    class Element(CombinatorialFreeModule.Element):
        @cached_method
        def is_virtual(self):
            return any(x < 0 for x in self.coefficients(sort=False))

        @cached_method
        def is_molecular(self):
            return len(self.coefficients(sort=False)) == 1 and self.coefficients(sort=False)[0] == 1

        @cached_method
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
                # maybe move degree_on_basis to MolecularSpecies (and it will become "grade")
                n = self.parent().degree_on_basis(M)
                powvecs = IntegerVectors(n, len(Gm))
                for vec in powvecs:
                    coeff = 1
                    for fN, exponent in zip(Gm.values(), vec):
                        coeff *= fN ** exponent
                    N_list = list(chain.from_iterable([[N._group for _ in range(c)] for N, c in zip(Gm.keys(), vec)]))
                    R = self.parent()(wreath_imprimitive_general(N_list, M._group))
                    term += coeff * R
                res += fM * term
            return res


def wreath_imprimitive_general(G_list, H):
    r"""
    H([G[0]] - [G[1]] - ... - [G[m]])
    """
    if len(G_list) != H.degree():
        raise ValueError(f"length of G_list (= {len(G_list)}) must be equal to degree of H (= {H.degree()})")
    gens, dlist = [], []
    dsum = 0
    for G in G_list:
        dlist.append(list(range(dsum + 1, dsum + G.degree() + 1)))
        dsum += G.degree()

    # First find gens of H
    for gen in H.gens():
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
        for gen in G_list[i].gens():
            images = libgap.MappingPermListList(G_list[i].domain().list(),dlist[i])
            gens.append(gen ** images)

    # Finally create group
    return PermutationGroup(gens, domain = range(1, dsum + 1))
