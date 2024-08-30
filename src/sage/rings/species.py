from itertools import accumulate, chain
from sage.categories.graded_algebras_with_basis import GradedAlgebrasWithBasis
from sage.categories.infinite_enumerated_sets import InfiniteEnumeratedSets
from sage.categories.monoids import Monoids
from sage.categories.sets_cat import cartesian_product
from sage.categories.sets_with_grading import SetsWithGrading
from sage.combinat.free_module import CombinatorialFreeModule
from sage.combinat.integer_vector import IntegerVectors
from sage.combinat.partition import Partitions
from sage.groups.perm_gps.constructor import PermutationGroupElement
from sage.groups.perm_gps.permgroup import PermutationGroup, PermutationGroup_generic
from sage.groups.perm_gps.permgroup_named import SymmetricGroup
from sage.libs.gap.libgap import libgap
from sage.misc.cachefunc import cached_method
from sage.monoids.indexed_free_monoid import (IndexedFreeAbelianMonoid,
                                              IndexedFreeAbelianMonoidElement,
                                              IndexedMonoid)
from sage.rings.integer_ring import ZZ
from sage.sets.finite_enumerated_set import FiniteEnumeratedSet
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
        Return the cached element, or add it
        if it doesn't exist.

        ``elm`` must implement the following methods:
            - ``_element_key`` - hashable type for dict lookup.
            - ``__eq__`` - to compare two elements.
            - ``_canonicalize`` - to preprocess the element.

        TESTS::

            sage: P = PolynomialSpecies(ZZ, "X, Y")
            sage: M = P._indices
            sage: G = PermutationGroup([[(1,2),(3,4),(5,6),(7,8,9,10)]]); G
            Permutation Group with generators [(1,2)(3,4)(5,6)(7,8,9,10)]
            sage: A = M(G, {1: [1,2,3,4], 2: [5,6,7,8,9,10]}); A
            {((1,2,3,4)(5,6)(7,8)(9,10),): ({5, 6, 7, 8}, {1, 2, 3, 4, 9, 10})}
            sage: C = M(G, {1: [1,2,5,6], 2: [3,4,7,8,9,10]}); C
            {((1,2,3,4)(5,6)(7,8)(9,10),): ({5, 6, 7, 8}, {1, 2, 3, 4, 9, 10})}
            sage: from sage.rings.species import ElementCache
            sage: E = ElementCache()
            sage: E._cache_get(A)
            {((1,2,3,4)(5,6)(7,8)(9,10),): ({5, 6, 7, 8}, {1, 2, 3, 4, 9, 10})}
            sage: E._cache_get(C)
            {((1,2,3,4)(5,6)(7,8)(9,10),): ({5, 6, 7, 8}, {1, 2, 3, 4, 9, 10})}
        """
        # TODO: Make _canonicalize optional.
        # Possibly the following works:
        # use getattr to check if name exists
        # use callable to check if it is a function
        # if both true, call _canonicalize
        key = elm._element_key()
        if key in self._cache:
            lookup = self._cache[key]
            for other_elm in lookup:
                if elm == other_elm:
                    return other_elm
            elm._canonicalize()
            lookup.append(elm)
        else:
            elm._canonicalize()
            self._cache[key] = [elm]
        return elm


class ConjugacyClassOfDirectlyIndecomposableSubgroups(Element):
    def __init__(self, parent, C):
        r"""
        A conjugacy class of directly indecomposable subgroups.

        TESTS::

            sage: from sage.rings.species import ConjugacyClassesOfDirectlyIndecomposableSubgroups
            sage: C = ConjugacyClassesOfDirectlyIndecomposableSubgroups()
            sage: G = C(PermutationGroup([[(1,2),(3,4)],[(1,2),(5,6)]]))
            sage: TestSuite(G).run()
        """
        Element.__init__(self, parent)
        self._C = C

    def __hash__(self):
        r"""
        Return the hash for the conjugacy class.

        TESTS::

            sage: from sage.rings.species import ConjugacyClassesOfDirectlyIndecomposableSubgroups
            sage: C = ConjugacyClassesOfDirectlyIndecomposableSubgroups()
            sage: G = PermutationGroup([[(1,2),(3,4),(5,6),(7,8,9,10)]]); G
            Permutation Group with generators [(1,2)(3,4)(5,6)(7,8,9,10)]
            sage: H = PermutationGroup([[(1,2,3,4),(5,6),(7,8),(9,10)]]); H
            Permutation Group with generators [(1,2,3,4)(5,6)(7,8)(9,10)]
            sage: hash(C(G)) == hash(C(H))
            True
        """
        return hash(self._C)

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: from sage.rings.species import ConjugacyClassesOfDirectlyIndecomposableSubgroups
            sage: C = ConjugacyClassesOfDirectlyIndecomposableSubgroups()
            sage: G = PermutationGroup([[(1,3),(4,7)], [(2,5),(6,8)], [(1,4),(2,5),(3,7)]])
            sage: C(G)
            ((5,6)(7,8), (1,2)(3,4), (1,3)(2,4)(5,6))
        """
        return f"{self._C.gens()}"

    def _element_key(self):
        r"""
        Return the cache lookup key for ``self``.

        TESTS::

            sage: from sage.rings.species import ConjugacyClassesOfDirectlyIndecomposableSubgroups
            sage: C = ConjugacyClassesOfDirectlyIndecomposableSubgroups()
            sage: G = PermutationGroup([[(1,3),(4,7)], [(2,5),(6,8)], [(1,4),(2,5),(3,7)]])
            sage: C(G)._element_key()
            (8, 8, (2, 2, 4))
        """
        return self._C.degree(), self._C.order(), tuple(len(orbit) for orbit in sorted(self._C.orbits(), key=len))

    @cached_method
    def _canonicalize(self):
        r"""
        Canonicalize this conjugacy class by sorting the orbits by
        length and making them consecutive.

        EXAMPLES::

            sage: from sage.rings.species import ConjugacyClassesOfDirectlyIndecomposableSubgroups
            sage: C = ConjugacyClassesOfDirectlyIndecomposableSubgroups()
            sage: G = PermutationGroup([[(1,3),(4,7)], [(2,5),(6,8)], [(1,4),(2,5),(3,7)]])
            sage: C(G)
            ((5,6)(7,8), (1,2)(3,4), (1,3)(2,4)(5,6))
        """
        if self._C == SymmetricGroup(0):
            return
        sorted_orbits = sorted((sorted(orbit) for orbit in self._C.orbits()),
                               key=len, reverse=True)
        pi = PermutationGroupElement(list(chain.from_iterable(sorted_orbits))).inverse()
        self._C = PermutationGroup(self._C.gens_small(), domain=self._C.domain()).conjugate(pi)

    def __eq__(self, other):
        r"""
        Return whether ``self`` is equal to ``other``.

        ``self`` is equal to ``other`` if they have the same degree (say `n`)
        and are conjugate within `S_n`.

        TESTS::

            sage: from sage.rings.species import ConjugacyClassesOfDirectlyIndecomposableSubgroups
            sage: C = ConjugacyClassesOfDirectlyIndecomposableSubgroups()
            sage: G = PermutationGroup([[(1,2),(3,4),(5,6),(7,8,9,10)]]); G
            Permutation Group with generators [(1,2)(3,4)(5,6)(7,8,9,10)]
            sage: H = PermutationGroup([[(1,2,3,4),(5,6),(7,8),(9,10)]]); H
            Permutation Group with generators [(1,2,3,4)(5,6)(7,8)(9,10)]
            sage: C(G) == C(H)
            True
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

        TESTS::

            sage: from sage.rings.species import ConjugacyClassesOfDirectlyIndecomposableSubgroups
            sage: C = ConjugacyClassesOfDirectlyIndecomposableSubgroups()
            sage: TestSuite(C).run(max_runs=5) # It takes far too long otherwise
        """
        Parent.__init__(self, category=InfiniteEnumeratedSets())
        ElementCache.__init__(self)

    @cached_method
    def an_element(self):
        r"""
        Return an element of ``self``.

        EXAMPLES::

            sage: from sage.rings.species import ConjugacyClassesOfDirectlyIndecomposableSubgroups
            sage: C = ConjugacyClassesOfDirectlyIndecomposableSubgroups()
            sage: C.an_element()
            ((),)
        """
        return self._element_constructor_(SymmetricGroup(1))

    def _element_constructor_(self, x):
        r"""
        ``x`` is an element of ``self`` or a group `H` such that
        `H` is directly indecomposable.

        EXAMPLES::

            sage: from sage.rings.species import ConjugacyClassesOfDirectlyIndecomposableSubgroups
            sage: C = ConjugacyClassesOfDirectlyIndecomposableSubgroups()
            sage: C(PermutationGroup([("a", "b", "c")]))
            ((1,2,3),)
            sage: C(PermutationGroup([(1, 3, 5)], domain=[1,3,5]))
            ((1,2,3),)
            sage: C(PermutationGroup([[(1,3),(4,7)],[(2,5),(6,8)], [(1,4),(2,5),(3,7)]]))
            ((5,6)(7,8), (1,2)(3,4), (1,3)(2,4)(5,6))
        """
        if parent(x) == self:
            return x
        if isinstance(x, PermutationGroup_generic):
            if len(x.disjoint_direct_product_decomposition()) > 1:
                raise ValueError(f"{x} is not directly indecomposable")
            mapping = {v: i for i, v in enumerate(x.domain(), 1)}
            normalized_gens = [[tuple(mapping[x] for x in cyc)
                                for cyc in gen.cycle_tuples()]
                            for gen in x.gens()]
            P = PermutationGroup(gens=normalized_gens)
            # Fix for SymmetricGroup(0)
            if x.degree() == 0:
                P = SymmetricGroup(0)
            elm = self.element_class(self, P)
            return self._cache_get(elm)
        raise ValueError(f"unable to convert {x} to {self}")

    def __iter__(self):
        r"""
        An iterator over all conjugacy classes of directly indecomposable
        subgroups.

        TESTS::

            sage: from sage.rings.species import ConjugacyClassesOfDirectlyIndecomposableSubgroups
            sage: C = ConjugacyClassesOfDirectlyIndecomposableSubgroups()
            sage: iterC = iter(C)
            sage: for i in range(5):
            ....:     print(next(iterC))
            ()
            ((),)
            ((1,2),)
            ((1,2,3),)
            ((2,3), (1,2,3))
        """
        n = 0
        while True:
            for G in SymmetricGroup(n).conjugacy_classes_subgroups():
                if len(G.disjoint_direct_product_decomposition()) <= 1:
                    yield self._element_constructor_(G)
            n += 1

    def __contains__(self, G):
        r"""
        Return if ``G`` is in ``self``.

        TESTS::

            sage: from sage.rings.species import ConjugacyClassesOfDirectlyIndecomposableSubgroups
            sage: C = ConjugacyClassesOfDirectlyIndecomposableSubgroups()
            sage: G = PermutationGroup([[(1,2)], [(3,4)]]); G
            Permutation Group with generators [(3,4), (1,2)]
            sage: G.disjoint_direct_product_decomposition()
            {{1, 2}, {3, 4}}
            sage: G in C
            False
        """
        if parent(G) == self:
            return True
        return len(G.disjoint_direct_product_decomposition()) <= 1

    def _repr_(self):
        r"""
        Return the string representation for ``self``.

        TESTS::

            sage: from sage.rings.species import ConjugacyClassesOfDirectlyIndecomposableSubgroups
            sage: C = ConjugacyClassesOfDirectlyIndecomposableSubgroups(); C
            Infinite set of conjugacy classes of directly indecomposable subgroups
        """
        return "Infinite set of conjugacy classes of directly indecomposable subgroups"

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

        TESTS::

            sage: A = AtomicSpecies("X")
            sage: G = PermutationGroup([[(1,3),(4,7)], [(2,5),(6,8)], [(1,4),(2,5),(3,7)]])
            sage: TestSuite(A(G)).run()
        """
        # Notice that a molecular species (and therefore an atomic species)
        # must be centered on a multicardinality, otherwise it wouldn't be
        # molecular. So it is kind of redundant to say "atomic species
        # centered on a given multicardinality."
        Element.__init__(self, parent)
        self._dis = dis
        self._dompart = domain_partition
        self._mc = tuple(len(v) for v in self._dompart)
        self._tc = sum(self._mc)

    def grade(self):
        r"""
        Return the grade of ``self``.
        """
        return self._mc

    def _element_key(self):
        r"""
        Return the cache lookup key for ``self``.

        TESTS::

            sage: At = AtomicSpecies("X, Y")
            sage: G = PermutationGroup([[(1,2),(3,4),(5,6),(7,8,9,10)]]); G
            Permutation Group with generators [(1,2)(3,4)(5,6)(7,8,9,10)]
            sage: A = At(G, {1: [1,2,3,4], 2: [5,6,7,8,9,10]}); A
            {((1,2,3,4)(5,6)(7,8)(9,10),): ({5, 6, 7, 8}, {1, 2, 3, 4, 9, 10})}
            sage: A._element_key()
            ((4, 6), ((1,2,3,4)(5,6)(7,8)(9,10),))
        """
        return self._mc, self._dis

    @cached_method
    def _canonicalize(self):
        r"""
        Canonicalize this atomic species by sorting the orbits by
        length and making them consecutive.

        EXAMPLES::

            sage: At = AtomicSpecies("X, Y")
            sage: G = PermutationGroup([[(1,2),(3,4),(5,6),(7,8,9,10)]]); G
            Permutation Group with generators [(1,2)(3,4)(5,6)(7,8,9,10)]
            sage: A = At(G, {1: [1,2,3,4], 2: [5,6,7,8,9,10]})
            sage: A._dompart
            (frozenset({5, 6, 7, 8}), frozenset({1, 2, 3, 4, 9, 10}))
            sage: C = At(G, {1: [1,2,5,6], 2: [3,4,7,8,9,10]})
            sage: C._dompart
            (frozenset({5, 6, 7, 8}), frozenset({1, 2, 3, 4, 9, 10}))
        """
        # The canonicalization is done in the element constructor.
        pass

    def __hash__(self):
        r"""
        Return the hash of the atomic species.

        TESTS::

            sage: At = AtomicSpecies("X, Y")
            sage: G = PermutationGroup([[(1,2),(3,4),(5,6),(7,8,9,10)]]); G
            Permutation Group with generators [(1,2)(3,4)(5,6)(7,8,9,10)]
            sage: H = PermutationGroup([[(1,2,3,4),(5,6),(7,8),(9,10)]]); H
            Permutation Group with generators [(1,2,3,4)(5,6)(7,8)(9,10)]
            sage: A = At(G, {1: [1,2,3,4], 2: [5,6,7,8,9,10]}); A
            {((1,2,3,4)(5,6)(7,8)(9,10),): ({5, 6, 7, 8}, {1, 2, 3, 4, 9, 10})}
            sage: B = At(H, {1: [1,2,3,4], 2: [5,6,7,8,9,10]}); B
            {((1,2,3,4)(5,6)(7,8)(9,10),): ({1, 2, 3, 4}, {5, 6, 7, 8, 9, 10})}
            sage: C = At(G, {1: [1,2,5,6], 2: [3,4,7,8,9,10]}); C
            {((1,2,3,4)(5,6)(7,8)(9,10),): ({5, 6, 7, 8}, {1, 2, 3, 4, 9, 10})}
            sage: hash(A) == hash(B)
            True
            sage: hash(A) == hash(C)
            True
        """
        return hash(self._element_key())

    def __eq__(self, other):
        r"""
        Two atomic species are equal if the underlying groups are conjugate,
        and their domain partitions are equal under the conjugating element.

        TESTS::

            sage: At = AtomicSpecies("X, Y")
            sage: G = PermutationGroup([[(1,2),(3,4),(5,6),(7,8,9,10)]]); G
            Permutation Group with generators [(1,2)(3,4)(5,6)(7,8,9,10)]
            sage: H = PermutationGroup([[(1,2,3,4),(5,6),(7,8),(9,10)]]); H
            Permutation Group with generators [(1,2,3,4)(5,6)(7,8)(9,10)]
            sage: A = At(G, {1: [1,2,3,4], 2: [5,6,7,8,9,10]})
            sage: B = At(H, {1: [1,2,3,4], 2: [5,6,7,8,9,10]})
            sage: C = At(G, {1: [1,2,5,6], 2: [3,4,7,8,9,10]})
            sage: A != B
            True
            sage: A == C
            True
        """
        if parent(self) != parent(other):
            return False
        # Check if multicardinalities match
        if self._mc != other._mc:
            return False
        # If they do, construct the mapping between the groups
        selflist, otherlist = [], []
        for i in range(self.parent()._k):
            if len(self._dompart[i]) != len(other._dompart[i]):
                return False
            selflist.extend(sorted(list(self._dompart[i])))
            otherlist.extend(sorted(list(other._dompart[i])))
        mapping = libgap.MappingPermListList(selflist, otherlist)
        G = PermutationGroup(gap_group=libgap.ConjugateGroup(self._dis._C, mapping),
                             domain=self._dis._C.domain())
        # The conjugated group must be exactly equal to the other group
        return G == other._dis._C

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        TESTS::

            sage: At = AtomicSpecies("X, Y")
            sage: G = PermutationGroup([[(1,2),(3,4),(5,6),(7,8,9,10)]]); G
            Permutation Group with generators [(1,2)(3,4)(5,6)(7,8,9,10)]
            sage: A = At(G, {1: [1,2,3,4], 2: [5,6,7,8,9,10]}); A
            {((1,2,3,4)(5,6)(7,8)(9,10),): ({5, 6, 7, 8}, {1, 2, 3, 4, 9, 10})}
            sage: A = At(G, {2: [1,2,3,4,5,6,7,8,9,10]}); A
            {((1,2,3,4)(5,6)(7,8)(9,10),): ({}, {1, 2, 3, 4, 5, 6, 7, 8, 9, 10})}
        """
        if self.parent()._k == 1:
            return "{" + f"{self._dis}" + "}"
        dompart = ', '.join("{" + repr(sorted(b))[1:-1] + "}"
                           for b in self._dompart)
        return "{" + f"{self._dis}: ({dompart})" + "}"


class AtomicSpecies(UniqueRepresentation, Parent, ElementCache):
    @staticmethod
    def __classcall__(cls, names):
        """
        Normalize the arguments.

        TESTS::

            sage: A1 = AtomicSpecies("X")
            sage: A2 = AtomicSpecies("Y")
            sage: A3 = AtomicSpecies("X, Y")
            sage: A4 = AtomicSpecies(["X", "Y"])
            sage: A1 is A2
            False
            sage: A3 is A4
            True
        """
        from sage.structure.category_object import normalize_names
        names = normalize_names(-1, names)
        return super().__classcall__(cls, names)

    def __init__(self, names):
        r"""
        Infinite set of multivariate atomic species.

        INPUT:

        - ``names`` -- an iterable of ``k`` strings

        TESTS::

            sage: At1 = AtomicSpecies(["X"])
            sage: At2 = AtomicSpecies(["X", "Y"])
            sage: TestSuite(At1).run(skip="_test_graded_components")
            sage: TestSuite(At2).run(skip="_test_graded_components")
        """
        category = SetsWithGrading().Infinite()
        Parent.__init__(self, names=names, category=category)
        ElementCache.__init__(self)
        self._k = len(names)
        self._grading_set = IntegerVectors(length=self._k)
        self._renamed = set()  # the degrees that have been renamed already

    @cached_method
    def an_element(self):
        """
        Return an element of ``self``.

        TESTS::

            sage: At1 = AtomicSpecies("X")
            sage: At2 = AtomicSpecies("X, Y")
            sage: At1.an_element()
            E_2
            sage: At2.an_element()
            {((1,2)(3,4),): ({1, 2}, {3, 4})}
        """
        G = PermutationGroup([[(2 * i - 1, 2 * i) for i in range(1, self._k + 1)]])
        m = {i: [2 * i - 1, 2 * i] for i in range(1, self._k + 1)}
        return self._element_constructor_(G, m)

    def _element_constructor_(self, G, pi=None):
        r"""
        Construct the `k`-variate atomic species with the given data.

        INPUT:

        - ``G`` - an element of ``self`` (in this case pi must be ``None``)
          or a permutation group.
        - ``pi`` - a dict mapping sorts to iterables whose union is the domain.
          If `k=1`, `pi` can be omitted.
        """
        if parent(G) == self:
            if pi is not None:
                raise ValueError("cannot reassign sorts to an atomic species")
            return G
        if not isinstance(G, PermutationGroup_generic):
            raise ValueError(f"{G} must be a permutation group")
        if pi is None:
            if self._k == 1:
                pi = {1: G.domain()}
            else:
                raise ValueError("the assignment of sorts to the domain elements must be provided")
        if not set(pi.keys()).issubset(range(1, self._k + 1)):
            raise ValueError(f"keys of pi (={pi.keys()}) must be in the range [1, {self._k}]")
        if sum(len(p) for p in pi.values()) != len(G.domain()) or set(chain.from_iterable(pi.values())) != set(G.domain()):
            raise ValueError(f"values of pi (={pi.values()}) must partition the domain of G (={G.domain()})")
        for orbit in G.orbits():
            if not any(set(orbit).issubset(p) for p in pi.values()):
                raise ValueError(f"All elements of orbit {orbit} must have the same sort")
        dis_elm = ConjugacyClassesOfDirectlyIndecomposableSubgroups()(G)
        mapping = {v: i for i, v in enumerate(G.domain(), 1)}
        mapping2 = PermutationGroupElement([mapping[e] for o in sorted(G.orbits(), key=len, reverse=True)
                                            for e in o]).inverse()
        dpart = [frozenset() for _ in range(self._k)]
        for k, v in pi.items():
            dpart[k - 1] = frozenset(mapping2(mapping[x]) for x in v)
        elm = self._cache_get(self.element_class(self, dis_elm, tuple(dpart)))
        if elm._tc not in self._renamed:
            self._rename(elm._tc)
        return elm

    def _rename(self, n):
        r"""
        Names for common species.

        EXAMPLES::

            sage: P = PolynomialSpecies(ZZ, ["X", "Y"])
            sage: P(SymmetricGroup(4), {1: range(1, 5)})
            E_4(X)
            sage: P(SymmetricGroup(4), {2: range(1, 5)})
            E_4(Y)
            sage: P(CyclicPermutationGroup(4), {1: range(1, 5)})
            C_4(X)
            sage: P(CyclicPermutationGroup(4), {2: range(1, 5)})
            C_4(Y)
            sage: P(DihedralGroup(4), {1: range(1, 5)})
            P_4(X)
            sage: P(DihedralGroup(4), {2: range(1, 5)})
            P_4(Y)
            sage: P(AlternatingGroup(4), {1: range(1, 5)})
            Eo_4(X)
            sage: P(AlternatingGroup(4), {2: range(1, 5)})
            Eo_4(Y)
        """
        from sage.groups.perm_gps.permgroup import PermutationGroup
        from sage.groups.perm_gps.permgroup_named import (AlternatingGroup,
                                                          CyclicPermutationGroup,
                                                          DihedralGroup,
                                                          SymmetricGroup)

        # prevent infinite recursion in self._element_constructor_
        self._renamed.add(n)
        for i in range(self._k):
            if n == 1:
                self(SymmetricGroup(1), {i+1: [1]}).rename(self._names[i])

            if self._k == 1:
                sort = ""
            else:
                sort = f"({self._names[i]})"

            if n >= 2:
                self(SymmetricGroup(n),
                    {i+1: range(1, n+1)}).rename(f"E_{n}" + sort)

            if n >= 3:
                self(CyclicPermutationGroup(n),
                    {i+1: range(1, n+1)}).rename(f"C_{n}" + sort)

            if n >= 4:
                self(DihedralGroup(n),
                    {i+1: range(1, n+1)}).rename(f"P_{n}" + sort)

            if n >= 4:
                self(AlternatingGroup(n),
                    {i+1: range(1, n+1)}).rename(f"Eo_{n}" + sort)

            if n >= 4 and not n % 2:
                gens = [[(i, n-i+1) for i in range(1, n//2 + 1)],
                        [(i, i+1) for i in range(1, n, 2)]]
                self(PermutationGroup(gens),
                    {i+1: range(1, n+1)}).rename(f"Pb_{n}" + sort)

    def __contains__(self, x):
        r"""
        Return if ``x`` is in ``self``.

        TESTS::

            sage: A = AtomicSpecies("X")
            sage: G = PermutationGroup([[(1,2)], [(3,4)]]); G
            Permutation Group with generators [(3,4), (1,2)]
            sage: G.disjoint_direct_product_decomposition()
            {{1, 2}, {3, 4}}
            sage: G in A
            False
        """
        if parent(x) == self:
            return True
        G, pi = None, None
        if isinstance(x, PermutationGroup_generic):
            if self._k == 1:
                G = x
                pi = {1: G.domain()}
            else:
                raise ValueError("the assignment of sorts to the domain elements must be provided")
        else:
            G, pi = x
        if not isinstance(G, PermutationGroup_generic):
            raise ValueError(f"{G} must be a permutation group")
        if not set(pi.keys()).issubset(range(1, self._k + 1)):
            raise ValueError(f"keys of pi (={pi.keys()}) must be in the range [1, {self._k}]")
        if sum(len(p) for p in pi.values()) != len(G.domain()) or set(chain.from_iterable(pi.values())) != set(G.domain()):
            raise ValueError(f"values of pi (={pi.values()}) must partition the domain of G (={G.domain()})")
        for orbit in G.orbits():
            if not any(set(orbit).issubset(p) for p in pi.values()):
                raise ValueError(f"All elements of orbit {orbit} must have the same sort")
        return len(G.disjoint_direct_product_decomposition()) <= 1

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        TESTS::

            sage: AtomicSpecies("X")
            Atomic species in X
            sage: AtomicSpecies("X, Y")
            Atomic species in X, Y
        """
        if len(self._names) == 1:
            return f"Atomic species in {self._names[0]}"
        return f"Atomic species in {', '.join(self._names)}"

    Element = AtomicSpeciesElement


class MolecularSpecies(IndexedFreeAbelianMonoid, ElementCache):
    @staticmethod
    def __classcall__(cls, indices, prefix=None, **kwds):
        return super(IndexedMonoid, cls).__classcall__(cls, indices, prefix, **kwds)

    def __init__(self, indices, prefix=None, **kwds):
        r"""
        Infinite set of multivariate molecular species.

        INPUT:

        - ``indices`` -- the underlying set of atomic species indexing the monoid

        TESTS::

            sage: P1 = PolynomialSpecies(ZZ, "X")
            sage: P2 = PolynomialSpecies(ZZ, ["X", "Y"])
            sage: TestSuite(P1._indices).run(skip="_test_graded_components")
            sage: TestSuite(P2._indices).run(skip="_test_graded_components")
        """
        category = Monoids() & SetsWithGrading().Infinite()
        IndexedFreeAbelianMonoid.__init__(self, indices, prefix=prefix, category=category, **kwds)
        ElementCache.__init__(self)
        self._k = indices._k

    def _project(self, G, pi, part):
        r"""
        Project `G` onto a subset ``part`` of its domain.

        ``part`` must be a union of cycles, but this is not checked.

        TESTS::

            sage: P = PolynomialSpecies(ZZ, ["X", "Y"])
            sage: M = P._indices
            sage: G = PermutationGroup([[(1,2),(3,4)], [(5,6)]]); G
            Permutation Group with generators [(5,6), (1,2)(3,4)]
            sage: parts = G.disjoint_direct_product_decomposition(); parts
            {{1, 2, 3, 4}, {5, 6}}
            sage: pi = {1: [1,2,3,4], 2: [5,6]}
            sage: M._project(G, pi, parts[0])
            (Permutation Group with generators [(1,2)(3,4)], {1: [1, 2, 3, 4]})
            sage: M._project(G, pi, parts[1])
            (Permutation Group with generators [(5,6)], {2: [5, 6]})
        """
        restricted_gens = [[cyc for cyc in gen.cycle_tuples() if cyc[0] in part]
                           for gen in G.gens()]
        restricted_gens = [gen for gen in restricted_gens if gen]
        mapping = dict()
        for k, v in pi.items():
            es = [e for e in v if e in part]
            if es:
                mapping[k] = es
        return PermutationGroup(gens=restricted_gens, domain=part), mapping

    def _element_constructor_(self, G, pi=None):
        r"""
        Construct the `k`-variate molecular species with the given data.

        INPUT:

        - ``G`` - an element of ``self`` (in this case pi must be ``None``)
          or a permutation group, or a pair ``(X, a)`` consisting of a
          finite set and an action.
        - ``pi`` - a dict mapping sorts to iterables whose union is the
          domain of ``G`` (if ``G`` is a permutation group) or `X` (if ``G``)
          is a pair ``(X, a)``. If `k=1`, `pi` can be omitted.

        If `G = (X, a)`, then `X` should be a finite set and `a` a transitive
        action of `G` on `X`.

        TESTS::

            sage: P = PolynomialSpecies(ZZ, ["X", "Y"])
            sage: M = P._indices
            sage: M(CyclicPermutationGroup(4), {1: [1,2], 2: [3,4]})
            Traceback (most recent call last):
            ...
            ValueError: All elements of orbit (1, 2, 3, 4) must have the same sort
        """
        if parent(G) == self:
            if pi is not None:
                raise ValueError("cannot reassign sorts to a molecular species")
            return G
        if isinstance(G, PermutationGroup_generic):
            if pi is None:
                if self._k == 1:
                    pi = {1: G.domain()}
                else:
                    raise ValueError("the assignment of sorts to the domain elements must be provided")
            domain_partition = G.disjoint_direct_product_decomposition()
            elm = self.one()
            for part in domain_partition:
                elm *= self.gen(self._project(G, pi, part))
            return elm
        # Assume G is a tuple (X, a)
        X, a = G
        if pi is None:
            if self._k == 1:
                pi = {1: X}
            else:
                raise ValueError("the assignment of sorts to the domain elements must be provided")
        L = [None for _ in range(self._k)]
        for k, v in pi.items():
            L[k - 1] = list(v)
        # Create group
        # TODO: Is this correct?
        S = SymmetricGroup(list(chain.from_iterable(L))).young_subgroup([len(v) for v in L])
        H = PermutationGroup(S.gens(), action=a, domain=X)
        if len(H.orbits()) > 1:
            # Then it is not transitive
            raise ValueError("Action is not transitive")
        stabG = PermutationGroup([g for g in S.gens() if a(g, H.orbits()[0][0]) == H.orbits()[0][0]], domain=S.domain())
        return self(stabG, pi)

    @cached_method
    def one(self):
        r"""
        Return the one of this monoid.

        EXAMPLES::

            sage: P = PolynomialSpecies(ZZ, "X, Y")
            sage: P._indices.one()
            1
        """
        elm = super().one()
        elm._group = SymmetricGroup(0)
        elm._dompart = tuple()
        elm._mc = tuple(0 for _ in range(self._k))
        elm._tc = 0
        return elm

    def gen(self, x):
        r"""
        Create the molecular species from an atomic species.

        EXAMPLES::

            sage: P = PolynomialSpecies(ZZ, "X, Y")
            sage: At = AtomicSpecies("X, Y")
            sage: G = PermutationGroup([[(1,2),(3,4)]]); G
            Permutation Group with generators [(1,2)(3,4)]
            sage: pi = {1: [1,2], 2: [3,4]}
            sage: E2XY = At(G, pi)
            sage: At(G, pi).rename("E_2(XY)"); E2XY
            E_2(XY)
            sage: m = P._indices.gen((G, pi)); m
            E_2(XY)
            sage: type(m)
            <class 'sage.rings.species.MolecularSpecies_with_category.element_class'>

            sage: M = P._indices
            sage: m = M(SymmetricGroup(6).young_subgroup([2, 2, 2]), {1: [1,2], 2: [3,4,5,6]})
            sage: list(m)
            [(E_2(X), 1), (E_2(Y), 2)]
        """
        if x not in self._indices:
            raise IndexError(f"{x} is not in the index set")
        if isinstance(x, (PermutationGroup_generic, AtomicSpecies.Element)):
            at = self._indices(x)
        else:
            at = self._indices(x[0], x[1])
        elm = self._cache_get(self.element_class(self, {at: ZZ.one()}))
        if elm._group is None:
            elm._group = at._dis._C
            elm._dompart = at._dompart
            elm._mc = at._mc
            elm._tc = at._tc
        return elm

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        TESTS::

            sage: PolynomialSpecies(ZZ, "X")._indices
            Molecular species in X
            sage: PolynomialSpecies(ZZ, "X, Y")._indices
            Molecular species in X, Y
        """
        if len(self._indices._names) == 1:
            return f"Molecular species in {self._indices._names[0]}"
        return f"Molecular species in {', '.join(self._indices._names)}"

    class Element(IndexedFreeAbelianMonoidElement):
        def __init__(self, F, x):
            r"""
            Initialize a molecular species with no group information.
            """
            super().__init__(F, x)
            self._group = None
            self._dompart = None
            self._mc = None
            self._tc = None

        def _assign_group_info(self, other):
            r"""
            Assign the group info of ``other`` to ``self``.
            """
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
            if elm1._tc == 0:
                self._assign_group_info(elm2)
                return
            if elm2._tc == 0:
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
            for gen in gens2:
                gens.append([tuple(elm1._tc + k for k in cyc) for cyc in gen.cycle_tuples()])
            self._group = PermutationGroup(gens, domain=range(1, elm1._tc + elm2._tc + 1))
            self._dompart = list(elm1._dompart)
            for i in range(elm2.parent()._k):
                self._dompart[i] = frozenset(list(self._dompart[i]) + [elm1._tc + e for e in elm2._dompart[i]])
            self._dompart = tuple(self._dompart)

        def __floordiv__(self, elt):
            raise NotImplementedError("Cannot cancel in this monoid")

        def _mul_(self, other):
            r"""
            Multiply ``self`` by ``other``.

            TESTS::

                sage: P = PolynomialSpecies(ZZ, "X, Y")
                sage: M = P._indices
                sage: E2X = M(SymmetricGroup(2), {1: [1,2]}); E2X
                E_2(X)
                sage: (E2X._group, E2X._dompart, E2X._mc, E2X._tc)
                (Permutation Group with generators [(1,2)],
                (frozenset({1, 2}), frozenset()),
                (2, 0),
                2)
                sage: E2Y = M(SymmetricGroup(2), {2: [1,2]}); E2Y
                E_2(Y)
                sage: (E2Y._group, E2Y._dompart, E2Y._mc, E2Y._tc)
                (Permutation Group with generators [(1,2)],
                (frozenset(), frozenset({1, 2})),
                (0, 2),
                2)
                sage: E2XE2Y = E2X * E2Y; E2XE2Y
                E_2(X)*E_2(Y)
                sage: (E2XE2Y._group, E2XE2Y._dompart, E2XE2Y._mc, E2XE2Y._tc)
                (Permutation Group with generators [(3,4), (1,2)],
                (frozenset({1, 2}), frozenset({3, 4})),
                (2, 2),
                4)
            """
            res = super()._mul_(other)
            elm = self.parent()._cache_get(res)
            if elm._group is None:
                elm._elmmul(self, other)
                elm._canonicalize()
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

            TESTS::

                sage: P = PolynomialSpecies(ZZ, "X")
                sage: M = P._indices
                sage: E2 = M(SymmetricGroup(2), {1: [1,2]})
                sage: (E2._group, E2._dompart, E2._mc, E2._tc)
                (Permutation Group with generators [(1,2)], (frozenset({1, 2}),), (2,), 2)
                sage: E2_3 = E2 ^ 3
                sage: (E2_3._group, E2_3._dompart, E2_3._mc, E2_3._tc)
                (Permutation Group with generators [(5,6), (3,4), (1,2)],
                (frozenset({1, 2, 3, 4, 5, 6}),),
                (6,),
                6)
            """
            res = super().__pow__(n)
            elm = self.parent()._cache_get(res)
            if elm._group is None:
                elm._elmexp(self, n)
                elm._canonicalize()
            return elm

        def _element_key(self):
            r"""
            Return the cache lookup key for ``self``.

            TESTS::

                sage: P = PolynomialSpecies(ZZ, "X")
                sage: M = P._indices
                sage: E2 = M(SymmetricGroup(2), {1: [1,2]})
                sage: E2._element_key()
                E_2
            """
            return self

        @cached_method
        def _canonicalize(self):
            r"""
            Canonicalize this molecular species by sorting the orbits by
            length and making them consecutive.

            EXAMPLES::

                sage: P = PolynomialSpecies(ZZ, "X, Y")
                sage: M = P._indices
                sage: G = PermutationGroup([[(1,2),(3,4),(5,6),(7,8,9,10)]]); G
                Permutation Group with generators [(1,2)(3,4)(5,6)(7,8,9,10)]
                sage: A = M(G, {1: [1,2,3,4], 2: [5,6,7,8,9,10]})
                sage: A._dompart
                (frozenset({5, 6, 7, 8}), frozenset({1, 2, 3, 4, 9, 10}))
                sage: C = M(G, {1: [1,2,5,6], 2: [3,4,7,8,9,10]})
                sage: C._dompart
                (frozenset({5, 6, 7, 8}), frozenset({1, 2, 3, 4, 9, 10}))
            """
            if self._group is None or self._group == SymmetricGroup(0):
                return
            sorted_orbits = sorted([sorted(orbit) for orbit in self._group.orbits()], key=len, reverse=True)
            pi = PermutationGroupElement(list(chain.from_iterable(sorted_orbits))).inverse()
            self._group = self._group.conjugate(pi)
            self._dompart = tuple(frozenset(pi(k) for k in v) for v in self._dompart)

        def grade(self):
            r"""
            Return the grade of ``self``.

            EXAMPLES::

                sage: P = PolynomialSpecies(ZZ, "X, Y")
                sage: M = P._indices
                sage: G = PermutationGroup([[(1,2),(3,4),(5,6),(7,8,9,10)]]); G
                Permutation Group with generators [(1,2)(3,4)(5,6)(7,8,9,10)]
                sage: A = M(G, {1: [1,2,3,4], 2: [5,6,7,8,9,10]}); A
                {((1,2,3,4)(5,6)(7,8)(9,10),): ({5, 6, 7, 8}, {1, 2, 3, 4, 9, 10})}
                sage: A.grade()
                (4, 6)
            """
            return self._mc

        def domain(self):
            r"""
            Return the domain of ``self``.

            EXAMPLES::

                sage: P = PolynomialSpecies(ZZ, "X, Y")
                sage: M = P._indices
                sage: G = PermutationGroup([[(1,2),(3,4),(5,6),(7,8,9,10)]]); G
                Permutation Group with generators [(1,2)(3,4)(5,6)(7,8,9,10)]
                sage: A = M(G, {1: [1,2,3,4], 2: [5,6,7,8,9,10]}); A
                {((1,2,3,4)(5,6)(7,8)(9,10),): ({5, 6, 7, 8}, {1, 2, 3, 4, 9, 10})}
                sage: A.domain()
                {1, 2, 3, 4, 5, 6, 7, 8, 9, 10}
            """
            return FiniteEnumeratedSet(range(1, self._tc + 1))

        def _compose_with_singletons(self, base_ring, names, args):
            r"""
            Compute the inner sum of exercise 2.6.16 of BLL book.

            INPUT:

                - ``base_ring``, the base ring of the result

                - ``names``, the (flat) list of names of the result

                - ``args``, the sequence of compositions, each of
                  which sums to the corresponding cardinality of
                  ``self``. The number of ``args`` is equal to the
                  arity of ``self``.

            EXAMPLES::

                sage: P = PolynomialSpecies(ZZ, "X")
                sage: M = P._indices
                sage: C4 = M(CyclicPermutationGroup(4))
                sage: C4._compose_with_singletons(ZZ, "X, Y", [[2, 2]]) # X^2Y^2 + C2(XY)
                E_2(XY) + X^2*Y^2

                sage: P = PolynomialSpecies(ZZ, ["X", "Y"])
                sage: M = P._indices
                sage: F = M(PermutationGroup([[(1,2,3), (4,5,6)]]), {1: [1,2,3], 2: [4,5,6]})
                sage: F
                {((1,2,3)(4,5,6),): ({1, 2, 3}, {4, 5, 6})}
                sage: F._compose_with_singletons(ZZ, "X1, X2, X3, Y1, Y2", [[1, 1, 1], [2, 1]])
                6*X1*X2*X3*Y1^2*Y2

            TESTS::

                sage: P = PolynomialSpecies(ZZ, "X")
                sage: M = P._indices
                sage: M.one()._compose_with_singletons(ZZ, "X", [[1]])
                1
            """
            # TODO: No checks are performed right now, must be added.
            # Checks: all args in Compositions, sums must match cardinalities.

            # Create group of the composition
            Pn = PolynomialSpecies(base_ring, names)
            res = Pn.zero()
            # conjugate self._group so that [1..k] is sort 1, [k+1,..] is sort 2, so on
            conj = PermutationGroupElement(list(chain.from_iterable(self._dompart))).inverse()
            G = libgap.ConjugateGroup(self._group, conj)

            comp = list(chain.from_iterable(args))
            dpart = {i + 1: range(x - comp[i] + 1, x + 1) for i, x in enumerate(accumulate(comp))}
            # Create the double coset representatives.
            S_down = SymmetricGroup(sum(comp)).young_subgroup(comp)
            S_up = SymmetricGroup(self._tc).young_subgroup(self._mc)
            taus = libgap.DoubleCosetRepsAndSizes(S_up, S_down, G)
            # Sum over double coset representatives.
            for tau, _ in taus:
                H = libgap.Intersection(libgap.ConjugateGroup(G, tau), S_down)
                grp = PermutationGroup(gap_group=H, domain=self.domain())
                res += Pn(grp, dpart)
            return res

        def __call__(self, *args):
            r"""
            Substitute `M_1,\dots, M_k` into ``self``.

            The arguments must all have the same parent and must all
            be molecular.  The number of arguments must be equal to
            the arity of ``self``.

            The result is a molecular species, whose parent is the
            same as those of the arguments.

            EXAMPLES::

                sage: P = PolynomialSpecies(ZZ, ["X"])
                sage: M = P._indices
                sage: X = M(SymmetricGroup(1))
                sage: E2 = M(SymmetricGroup(2))
                sage: E2(X)
                E_2
                sage: X(E2)
                E_2
                sage: E2(E2)
                P_4

                sage: P = PolynomialSpecies(ZZ, ["X","Y"])
                sage: M = P._indices
                sage: X = M(SymmetricGroup(1), {1:[1]})
                sage: Y = M(SymmetricGroup(1), {2:[1]})
                sage: (X*Y)(X, Y^2)
                X*Y^2

            A multivariate example::

                sage: M1 = PolynomialSpecies(QQ, "X")._indices
                sage: M2 = PolynomialSpecies(QQ, "X, Y")._indices
                sage: C3 = M1(CyclicPermutationGroup(3))
                sage: X = M2(SymmetricGroup(1), {1: [1]})
                sage: Y = M2(SymmetricGroup(1), {2: [1]})
                sage: C3(X*Y)
                {((1,2,3)(4,5,6),): ({1, 2, 3}, {4, 5, 6})}

            TESTS::

                sage: P = PolynomialSpecies(QQ, ["X"])
                sage: P.one()()
                Traceback (most recent call last):
                ...
                ValueError: number of args must match arity of self
                """
            if len(args) != self.parent()._k:
                raise ValueError("number of args must match arity of self")
            if len(set(arg.parent() for arg in args)) > 1:
                raise ValueError("all args must have the same parent")
            if not all(isinstance(arg, MolecularSpecies.Element) for arg in args):
                raise ValueError("all args must be molecular species")

            gens = []

            # TODO: What happens if in F(G), G has a constant part? E(1+X)?
            Mlist = [None for _ in range(self._tc)]
            for i, v in enumerate(self._dompart):
                for k in v:
                    Mlist[k - 1] = args[i]
            starts = list(accumulate([M._tc for M in Mlist], initial=0))

            # gens from self
            for gen in self._group.gens():
                newgen = []
                for cyc in gen.cycle_tuples():
                    for k in range(1, Mlist[cyc[0] - 1]._tc + 1):
                        newgen.append(tuple(k + starts[i - 1] for i in cyc))
                gens.append(newgen)

            # gens from M_i and dompart
            dpart = {i: [] for i in range(1, args[0].parent()._k + 1)}
            for start, M in zip(starts, Mlist):
                for i, v in enumerate(M._dompart, 1):
                    dpart[i].extend([start + k for k in v])
                for gen in M._group.gens():
                    gens.append([tuple(start + k for k in cyc) for cyc in gen.cycle_tuples()])

            return args[0].parent()(PermutationGroup(gens, domain=range(1, starts[-1] + 1)), dpart)


class PolynomialSpecies(CombinatorialFreeModule):
    def __classcall__(cls, base_ring, names):
        r"""
        Normalize the arguments.

        TESTS::

            sage: P1 = PolynomialSpecies(ZZ, "X, Y")
            sage: P2 = PolynomialSpecies(ZZ, "X, Y")
            sage: P3 = PolynomialSpecies(ZZ, ["X", "Z"])
            sage: P1 is P2
            True
            sage: P1 is P3
            False
        """
        from sage.structure.category_object import normalize_names
        names = normalize_names(-1, names)
        return super().__classcall__(cls, base_ring, names)

    def __init__(self, base_ring, names):
        r"""
        Ring of `k`-variate polynomial (virtual) species.

        TESTS::

            sage: P = PolynomialSpecies(ZZ, "X")
            sage: TestSuite(P).run()
            sage: P2 = PolynomialSpecies(ZZ, "X, Y")
            sage: TestSuite(P2).run()
        """
        # should we pass a category to basis_keys?
        A = AtomicSpecies(names)
        basis_keys = MolecularSpecies(A, prefix='', bracket=False)
        category = GradedAlgebrasWithBasis(base_ring).Commutative()
        CombinatorialFreeModule.__init__(self, base_ring,
                                         basis_keys=basis_keys,
                                         category=category,
                                         element_class=self.Element,
                                         prefix='', bracket=False)
        self._k = len(names)

    def degree_on_basis(self, m):
        r"""
        Return the degree of the molecular species indexed by ``m``.

        EXAMPLES::

            sage: P = PolynomialSpecies(ZZ, ["X", "Y"])
            sage: E4X = P(SymmetricGroup(4), {1: range(1, 5)}); E4X
            E_4(X)
            sage: E4Y = P(SymmetricGroup(4), {2: range(1, 5)}); E4Y
            E_4(Y)
            sage: P.degree_on_basis(E4X.support()[0])
            4
            sage: P.degree_on_basis(E4Y.support()[0])
            4
        """
        return m._tc

    def _element_constructor_(self, G, pi=None):
        r"""
        Construct the `k`-variate polynomial species with the given data.

        INPUT:

        - ``G`` - an element of ``self`` (in this case pi must be ``None``)
          or a permutation group.
        - ``pi`` - a dict mapping sorts to iterables whose union is the domain.
          If `k=1`, `pi` can be omitted.

        If `G = (X, a)`, then `X` should be a finite set and `a` an action of
        `G` on `X`.

        EXAMPLES::

            sage: P = PolynomialSpecies(ZZ, ["X", "Y"])
            sage: P(SymmetricGroup(4).young_subgroup([2, 2]), {1: [1,2], 2: [3,4]})
            E_2(X)*E_2(Y)

            sage: X = SetPartitions(4, 2)
            sage: a = lambda g, x: SetPartition([[g(e) for e in b] for b in x])
            sage: P((X, a), {1: [1,2], 2: [3,4]})
            X^2*E_2(Y) + X^2*Y^2 + E_2(X)*Y^2 + E_2(X)*E_2(Y)
        """
        if parent(G) == self:
            if pi is not None:
                raise ValueError("cannot reassign sorts to a polynomial species")
            return G
        if isinstance(G, PermutationGroup_generic):
            if pi is None:
                if self._k == 1:
                    return self._from_dict({self._indices(G): ZZ.one()})
                raise ValueError("the assignment of sorts to the domain elements must be provided")
            return self._from_dict({self._indices(G, pi): ZZ.one()})
        # Assume G is a tuple (X, a)
        X, a = G
        if pi is None:
            if self._k == 1:
                pi = {1: X}
            else:
                raise ValueError("the assignment of sorts to the domain elements must be provided")
        # Make iteration over values of pi deterministic
        L = [None for _ in range(self._k)]
        for k, v in pi.items():
            L[k - 1] = list(v)
        # Create group
        # TODO: Is this correct?
        S = SymmetricGroup(list(chain.from_iterable(L))).young_subgroup([len(v) for v in L])
        H = PermutationGroup(S.gens(), action=a, domain=X)
        res = self.zero()
        for orbit in H.orbits():
            stabG = PermutationGroup([g for g in S.gens() if a(g, orbit[0]) == orbit[0]], domain=S.domain())
            res += self._from_dict({self._indices(stabG, pi): ZZ.one()})
        return res

    @cached_method
    def one_basis(self):
        r"""
        Returns SymmetricGroup(0), which indexes the one of this algebra,
        as per :meth:`AlgebrasWithBasis.ParentMethods.one_basis`.

        EXAMPLES::

            sage: P = PolynomialSpecies(ZZ, "X")
            sage: P.one_basis()
            1
            sage: P2 = PolynomialSpecies(ZZ, "X, Y")
            sage: P2.one_basis()
            1
        """
        return self._indices.one()

    @cached_method
    def an_element(self):
        """
        Return an element of ``self``.

            sage: P = PolynomialSpecies(ZZ, "X")
            sage: P.an_element()
            1
            sage: P2 = PolynomialSpecies(ZZ, "X, Y")
            sage: P2.an_element()
            1
        """
        return self.one()

    def product_on_basis(self, H, K):
        r"""
        Return the product of the basis elements indexed by `H` and `K`.

        EXAMPLES::

            sage: P = PolynomialSpecies(ZZ, "X")
            sage: L1 = [P(H) for H in SymmetricGroup(3).conjugacy_classes_subgroups()]
            sage: L2 = [P(H) for H in SymmetricGroup(2).conjugacy_classes_subgroups()]
            sage: matrix([[F * G for F in L1] for G in L2])
            [    X^5 X^3*E_2 C_3*X^2 E_3*X^2]
            [X^3*E_2 X*E_2^2 C_3*E_2 E_3*E_2]

        TESTS::

            sage: P = PolynomialSpecies(ZZ, "X")
            sage: X = P(SymmetricGroup(1))
            sage: type(list(X^2)[0][1])
            <class 'sage.rings.integer.Integer'>
        """
        return self._from_dict({H * K: ZZ(1)})

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: PolynomialSpecies(ZZ, "X")
            Polynomial species in X over Integer Ring
            sage: PolynomialSpecies(ZZ, "X, Y")
            Polynomial species in X, Y over Integer Ring
        """
        names = self._indices._indices._names
        if len(names) == 1:
            return f"Polynomial species in {names[0]} over {self.base_ring()}"
        return f"Polynomial species in {', '.join(names)} over {self.base_ring()}"

    @cached_method
    def powersum(self, s, n):
        r"""
        Return `P_n(X_s)`.

        EXAMPLES::

            sage: P = PolynomialSpecies(ZZ, "X")
            sage: P.powersum(1, 4)
            4*E_4 - 4*X*E_3 + 4*X^2*E_2 - X^4 - 2*E_2^2
        """
        assert n in ZZ and n > 0
        if n == 1:
            return self(SymmetricGroup(1), {s: [1]})
        return (ZZ(n) * self(SymmetricGroup(n), {s: range(1, n+1)})
                - sum(self(SymmetricGroup(i), {s: range(1, i+1)})
                      * self.powersum(s, n-i)
                      for i in range(1, n)))

    def exponential(self, multiplicities, degrees):
        r"""
        Return `E(\sum_i m_i X_i)` in the specified degrees.

        EXAMPLES::

            sage: P = PolynomialSpecies(QQ, ["X"])
            sage: P.exponential([3/2], [7])  # random
            3/2*E_7 + 3/4*X*E_6 - 3/16*X^2*E_5 + 3/32*X^3*E_4 - 15/256*E_3*X^4
             + 21/512*X^5*E_2 - 9/2048*X^7 - 15/128*X^3*E_2^2 - 3/8*E_2*E_4*X
             + 3/32*X*E_2^3 - 3/16*X*E_3^2 + 3/4*E_2*E_5 - 3/16*E_3*E_2^2
             + 3/4*E_3*E_4 + 9/32*E_3*E_2*X^2

        We support weights::

            sage: R.<q> = QQ[]
            sage: P = PolynomialSpecies(R, ["X"])
            sage: P.exponential([1], [2])
            E_2

            sage: P.exponential([1+q], [3])
            (q^3+1)*E_3 + (q^2+q)*X*E_2

            sage: P.exponential([1-q], [2])
            (-q^2+1)*E_2 + (q^2-q)*X^2

            sage: P = PolynomialSpecies(R, ["X", "Y"])
            sage: P.exponential([1, q], [2, 2])
            q^2*E_2(X)*E_2(Y)

        TESTS::

            sage: from sage.rings.lazy_species import LazySpecies
            sage: L = LazySpecies(QQ, "X")
            sage: E = L(lambda n: SymmetricGroup(n))
            sage: P = PolynomialSpecies(QQ, ["X"])

            sage: c = 3/2; all((E^c)[i] == P.exponential([c], [i]) for i in range(6))
            True

            sage: c = -5/3; all((E^c)[i] == P.exponential([c], [i]) for i in range(6))
            True

            sage: c = 0; all((E^c)[i] == P.exponential([c], [i]) for i in range(6))
            True

            sage: c = 1; all((E^c)[i] == P.exponential([c], [i]) for i in range(6))
            True

            sage: c = -1; all((E^c)[i] == P.exponential([c], [i]) for i in range(6))
            True

            sage: P = PolynomialSpecies(QQ, ["X"])
            sage: P.exponential([1], [0]).parent()
            Polynomial species in X over Rational Field
        """
        def stretch(c, k):
            r"""
            Return c
            """
            if callable(c):
                B = self.base_ring()
                return c(*[g ** k for g in B.gens() if g != B.one()])
            return c

        def factor(s, c, d):
            r"""
            Return `E(c X_s)_d`.

            We use Proposition 2 in Labelle, New combinatorial
            computational methods arising from pseudo-singletons.
            """
            return self.sum(~ mu.centralizer_size()
                            * self.prod(stretch(c, k)
                                        * self.powersum(s, k) for k in mu)
                            for mu in Partitions(d))

        return self.prod(factor(s+1, multiplicities[s], degrees[s])
                         for s in range(self._k))

    class Element(CombinatorialFreeModule.Element):
        def is_constant(self):
            """
            Return ``True`` if this is a constant polynomial species.

            EXAMPLES::

                sage: P = PolynomialSpecies(ZZ, ["X", "Y"])
                sage: X = P(SymmetricGroup(1), {1: [1]})
                sage: X.is_constant()
                False
                sage: (3*P.one()).is_constant()
                True
                sage: P(0).is_constant()
                True
                sage: (1 + X).is_constant()
                False
            """
            return self.is_zero() or not self.maximal_degree()

        def homogeneous_degree(self):
            """

            ..TODO::

               This implementation should not be necessary.

            EXAMPLES::

                sage: P = PolynomialSpecies(ZZ, ["X"])
                sage: C3 = P(CyclicPermutationGroup(3))
                sage: X = P(SymmetricGroup(1))
                sage: E2 = P(SymmetricGroup(2))
                sage: (E2*X + C3).homogeneous_degree()
                3
            """
            if not self.support():
                raise ValueError("the zero element does not have a well-defined degree")
            if not self.is_homogeneous():
                raise ValueError("element is not homogeneous")
            return self.parent().degree_on_basis(self.support()[0])

        def is_virtual(self):
            r"""
            Return if ``self`` is a virtual species.

            TESTS::

                sage: P = PolynomialSpecies(ZZ, ["X", "Y"])
                sage: X = P(SymmetricGroup(1), {1: [1]})
                sage: Y = P(SymmetricGroup(1), {2: [1]})
                sage: V = 2 * X - 3 * Y; V
                2*X - 3*Y
                sage: V.is_virtual()
                True
                sage: (X * Y).is_virtual()
                False
            """
            return any(x < 0 for x in self.coefficients(sort=False))

        def is_molecular(self):
            r"""
            Return if ``self`` is a molecular species.

            TESTS::

                sage: P = PolynomialSpecies(ZZ, ["X", "Y"])
                sage: X = P(SymmetricGroup(1), {1: [1]})
                sage: Y = P(SymmetricGroup(1), {2: [1]})
                sage: V = 2 * X - 3 * Y; V
                2*X - 3*Y
                sage: V.is_molecular()
                False
                sage: (2 * X).is_molecular()
                False
                sage: (X * Y).is_molecular()
                True
            """
            return len(self.coefficients(sort=False)) == 1 and self.coefficients(sort=False)[0] == 1

        def is_atomic(self):
            r"""
            Return if ``self`` is an atomic species.

            TESTS::

                sage: P = PolynomialSpecies(ZZ, ["X", "Y"])
                sage: X = P(SymmetricGroup(1), {1: [1]})
                sage: Y = P(SymmetricGroup(1), {2: [1]})
                sage: V = 2 * X - 3 * Y; V
                2*X - 3*Y
                sage: V.is_atomic()
                False
                sage: (2 * X).is_atomic()
                False
                sage: (X * Y).is_atomic()
                False
                sage: Y.is_atomic()
                True
            """
            return self.is_molecular() and len(self.support()[0]) == 1

        def hadamard_product(self, other):
            r"""
            Compute the hadamard product of ``self`` and ``other``.

            EXAMPLES:

            Exercise 2.1.9 from the BLL book::

                sage: P = PolynomialSpecies(ZZ, ["X"])
                sage: C3 = P(CyclicPermutationGroup(3))
                sage: X = P(SymmetricGroup(1))
                sage: E2 = P(SymmetricGroup(2))
                sage: C3.hadamard_product(C3)
                2*C_3
                sage: (X^3).hadamard_product(C3)
                2*X^3
                sage: (X*E2).hadamard_product(X*E2)
                X*E_2 + X^3

            TESTS::

                sage: C3.hadamard_product(-C3)
                -2*C_3
            """
            P = self.parent()
            if P is not other.parent():
                raise ValueError("the factors of a Hadamard product must have the same parent")

            res = P.zero()
            # we should first collect matching multicardinalities.
            for L, c in self:
                S = SymmetricGroup(L._tc).young_subgroup(L._mc)
                # conjugate L and R to match S
                conj_L = PermutationGroupElement(list(chain.from_iterable(L._dompart))).inverse()
                G = libgap.ConjugateGroup(L._group, conj_L)
                dpart = {i + 1: range(x - L._mc[i] + 1, x + 1) for i, x in enumerate(accumulate(L._mc))}
                for R, d in other:
                    if L._mc != R._mc:
                        continue
                    conj_R = PermutationGroupElement(list(chain.from_iterable(R._dompart))).inverse()
                    H = libgap.ConjugateGroup(R._group, conj_R)
                    taus = libgap.DoubleCosetRepsAndSizes(S, G, H)
                    # loop over representatives
                    new = P.zero()
                    for tau, _ in taus:
                        F = libgap.Intersection(libgap.ConjugateGroup(H, tau), G)
                        new += P(PermutationGroup(gap_group=F, domain=L.domain()), dpart)
                    res += c * d * new

            return res

        def _compose_with_weighted_singletons(self, names, multiplicities, degrees):
            r"""
            Compute the composition with
            `(\sum_j m_{1,j} X_{1,j}, \sum_j m_{2,j} X_{2,j}, \dots)`
            in the specified degrees.

            The `k`-sort species ``self`` should be homogeneous.

            INPUT:

                - ``names``, the (flat) list of names of the result

                - ``multiplicities``, a (flat) list of constants

                - ``degrees``, a `k`-tuple of compositions `c_1,
                  \dots, c_k`, such that the size of `c_i` is the
                  degree of self in sort `i`.

            EXAMPLES:

            Equation (2.5.41)::

                sage: P = PolynomialSpecies(QQ, ["X"])
                sage: E2 = P(SymmetricGroup(2))
                sage: E2._compose_with_weighted_singletons(["X"], [-1], [[2]])
                -E_2 + X^2

                sage: C4 = P(CyclicPermutationGroup(4))
                sage: C4._compose_with_weighted_singletons(["X"], [-1], [[4]])
                -C_4 + {((1,2)(3,4),): ({1, 2, 3, 4})}

            Exercise (2.5.17)::

                sage: C4._compose_with_weighted_singletons(["X", "Y"], [1, 1], [[2, 2]])
                E_2(XY) + X^2*Y^2
                sage: C4._compose_with_weighted_singletons(["X", "Y"], [1, 1], [[3, 1]])
                X^3*Y
                sage: C4._compose_with_weighted_singletons(["X", "Y"], [1, 1], [[4, 0]])
                C_4(X)

            Auger et al., Equation (4.60)::

                sage: C4._compose_with_weighted_singletons(["X", "Y"], [1, -1], [[2, 2]])
                -E_2(XY) + 2*X^2*Y^2

            TESTS::

                sage: (C4+E2^2)._compose_with_weighted_singletons(["X"], [-1], [[4]])
                -C_4 + {((1,2)(3,4),): ({1, 2, 3, 4})} + E_2^2 - 2*X^2*E_2 + X^4
            """
            P = self.parent()
            if not self.support():
                return P.zero()
            if not self.is_homogeneous():
                raise ValueError("element is not homogeneous")

            left = sum(c * M._compose_with_singletons(P.base_ring(),
                                                      names,
                                                      degrees)
                       for M, c in self)
            P = left.parent()
            right = P.exponential(multiplicities,
                                  list(chain.from_iterable(degrees)))
            return left.hadamard_product(right)

        def __call__(self, *args):
            """

            EXAMPLES::

                sage: P = PolynomialSpecies(QQ, ["X"])
                sage: X = P(SymmetricGroup(1))
                sage: E2 = P(SymmetricGroup(2))
                sage: E2(-X)
                -E_2 + X^2
                sage: E2(X^2)
                {((1,2)(3,4),): ({1, 2, 3, 4})}

                sage: E2(X + X^2)
                E_2 + X^3 + {((1,2)(3,4),): ({1, 2, 3, 4})}

                sage: P2 = PolynomialSpecies(QQ, ["X", "Y"])
                sage: X = P2(SymmetricGroup(1), {1:[1]})
                sage: Y = P2(SymmetricGroup(1), {2:[1]})
                sage: E2(X + Y)
                E_2(Y) + X*Y + E_2(X)

                sage: E2(X*Y)(E2(X), E2(Y))
                {((7,8), (5,6), (3,4), (1,2), (1,3)(2,4)(5,7)(6,8)): ({1, 2, 3, 4}, {5, 6, 7, 8})}

                sage: R.<q> = QQ[]
                sage: P = PolynomialSpecies(R, ["X"])
                sage: X = P(SymmetricGroup(1))
                sage: E2 = P(SymmetricGroup(2))
                sage: E2(q*X)
                q^2*E_2

            """
            P = self.parent()
            if len(args) != P._k:
                raise ValueError("number of args must match arity of self")
            if len(set(arg.parent() for arg in args)) > 1:
                raise ValueError("all args must have the same parent")

            P0 = args[0].parent()
            if not self.support():
                return P0.zero()

            args = [sorted(g, key=lambda x: x[0]._mc) for g in args]
            multiplicities = list(chain.from_iterable([[c for _, c in g] for g in args]))
            molecules = list(chain.from_iterable([[M for M, _ in g] for g in args]))
            F_degrees = sorted(set(M._mc for M, _ in self))
            names = ["X%s" % i for i in range(sum(len(arg) for arg in args))]

            result = P0.zero()
            for n in F_degrees:
                F = P.sum_of_terms((M, c) for M, c in self if M._mc == n)
                for degrees in cartesian_product([IntegerVectors(n_i, length=len(arg))
                                                  for n_i, arg in zip(n, args)]):
                    # each degree is a weak composition of the degree of F in sort i
                    FX = F._compose_with_weighted_singletons(names,
                                                             multiplicities,
                                                             degrees)
                    FG = [(M(*molecules), c) for M, c in FX]
                    result += P0.sum_of_terms(FG)
            return result
