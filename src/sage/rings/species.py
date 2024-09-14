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
                                              IndexedFreeAbelianMonoidElement)
from sage.rings.integer_ring import ZZ
from sage.structure.element import Element, parent
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.modules.free_module_element import vector


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

            sage: A = AtomicSpecies("X, Y")
            sage: G = PermutationGroup([[(1,2),(3,4),(5,6),(7,8,9,10)]]); G
            Permutation Group with generators [(1,2)(3,4)(5,6)(7,8,9,10)]
            sage: a = A(G, {1: [1,2,3,4], 2: [5,6,7,8,9,10]}); a
            {((1,2,3,4)(5,6)(7,8)(9,10),): ({5, 6, 7, 8}, {1, 2, 3, 4, 9, 10})}
            sage: b = A(G, {1: [1,2,5,6], 2: [3,4,7,8,9,10]}); b
            {((1,2,3,4)(5,6)(7,8)(9,10),): ({5, 6, 7, 8}, {1, 2, 3, 4, 9, 10})}
            sage: a is b
            True
            sage: from sage.rings.species import ElementCache
            sage: E = ElementCache()
            sage: E._cache_get(a)
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
            sage: C(G)  # indirect doctest
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
                and (d := self._C.degree()) == other._C.degree()
                and (self._C == other._C
                     or SymmetricGroup(d).are_conjugate(self._C, other._C)))


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

        EXAMPLES::

            sage: A = AtomicSpecies("X, Y")
            sage: G = PermutationGroup([[(1,2),(3,4),(5,6),(7,8,9,10)]])
            sage: a = A(G, {1: [1,2,3,4], 2: [5,6,7,8,9,10]})
            sage: a.grade()
            [4, 6]
        """
        S = self.parent().grading_set()
        return S(self._mc)

    def _element_key(self):
        r"""
        Return the cache lookup key for ``self``.

        TESTS::

            sage: A = AtomicSpecies("X, Y")
            sage: G = PermutationGroup([[(1,2),(3,4),(5,6),(7,8,9,10)]]); G
            Permutation Group with generators [(1,2)(3,4)(5,6)(7,8,9,10)]
            sage: a = A(G, {1: [1,2,3,4], 2: [5,6,7,8,9,10]}); a
            {((1,2,3,4)(5,6)(7,8)(9,10),): ({5, 6, 7, 8}, {1, 2, 3, 4, 9, 10})}
            sage: a._element_key()
            ((4, 6), ((1,2,3,4)(5,6)(7,8)(9,10),))
        """
        return self._mc, self._dis

    @cached_method
    def _canonicalize(self):
        r"""
        Canonicalize this atomic species by sorting the orbits by
        length and making them consecutive.

        EXAMPLES::

            sage: A = AtomicSpecies("X, Y")
            sage: G = PermutationGroup([[(1,2),(3,4),(5,6),(7,8,9,10)]]); G
            Permutation Group with generators [(1,2)(3,4)(5,6)(7,8,9,10)]
            sage: f = A(G, {1: [1,2,3,4], 2: [5,6,7,8,9,10]}); f
            {((1,2,3,4)(5,6)(7,8)(9,10),): ({5, 6, 7, 8}, {1, 2, 3, 4, 9, 10})}
            sage: A(G, {1: [1,2,5,6], 2: [3,4,7,8,9,10]}) is f  # indirect doctest
            True
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
        for i in range(self.parent()._arity):
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

            sage: A = AtomicSpecies("X, Y")
            sage: G = PermutationGroup([[(1,2),(3,4),(5,6),(7,8,9,10)]]); G
            Permutation Group with generators [(1,2)(3,4)(5,6)(7,8,9,10)]
            sage: a = A(G, {1: [1,2,3,4], 2: [5,6,7,8,9,10]}); a
            {((1,2,3,4)(5,6)(7,8)(9,10),): ({5, 6, 7, 8}, {1, 2, 3, 4, 9, 10})}
            sage: a = A(G, {2: [1,2,3,4,5,6,7,8,9,10]}); a
            {((1,2,3,4)(5,6)(7,8)(9,10),): ({}, {1, 2, 3, 4, 5, 6, 7, 8, 9, 10})}
        """
        if self.parent()._arity == 1:
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

        - ``names`` -- an iterable of ``k`` strings for the sorts of the species

        TESTS::

            sage: At1 = AtomicSpecies(["X"])
            sage: At2 = AtomicSpecies(["X", "Y"])
            sage: TestSuite(At1).run(skip="_test_graded_components")
            sage: TestSuite(At2).run(skip="_test_graded_components")
        """
        category = SetsWithGrading().Infinite()
        Parent.__init__(self, names=names, category=category)
        ElementCache.__init__(self)
        self._arity = len(names)
        self._renamed = set()  # the degrees that have been renamed already

    def grading_set(self):
        r"""
        Return the set of non-negative integer vectors, whose length is
        the arity of ``self``.

        EXAMPLES::

            sage: AtomicSpecies(["X"]).grading_set()
            Integer vectors of length 1
        """
        return IntegerVectors(length=self._arity)

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
        G = PermutationGroup([[(2 * i - 1, 2 * i) for i in range(1, self._arity + 1)]])
        m = {i: [2 * i - 1, 2 * i] for i in range(1, self._arity + 1)}
        return self._element_constructor_(G, m)

    def _element_constructor_(self, G, pi=None):
        r"""
        Construct the `k`-variate atomic species with the given data.

        INPUT:

        - ``G`` - an element of ``self`` (in this case pi must be ``None``)
          or a permutation group.
        - ``pi`` - a dict mapping sorts to iterables whose union is the domain.
          If `k=1`, `pi` can be omitted.


        EXAMPLES::

            sage: A = AtomicSpecies("X, Y")
            sage: A(DihedralGroup(5), {1: [1,2,3,4,5]})
            P_5(X)

            sage: G = PermutationGroup([[(1,2),(3,4,5,6)]])
            sage: a = A(G, {1: [1,2], 2: [3,4,5,6]}); a
            {((1,2,3,4)(5,6),): ({5, 6}, {1, 2, 3, 4})}
        """
        if parent(G) == self:
            if pi is not None:
                raise ValueError("cannot reassign sorts to an atomic species")
            return G
        if not isinstance(G, PermutationGroup_generic):
            raise ValueError(f"{G} must be a permutation group")
        if pi is None:
            if self._arity == 1:
                pi = {1: G.domain()}
            else:
                raise ValueError("the assignment of sorts to the domain elements must be provided")
        if not set(pi.keys()).issubset(range(1, self._arity + 1)):
            raise ValueError(f"keys of pi (={pi.keys()}) must be in the range [1, {self._arity}]")
        if sum(len(p) for p in pi.values()) != len(G.domain()) or set(chain.from_iterable(pi.values())) != set(G.domain()):
            raise ValueError(f"values of pi (={pi.values()}) must partition the domain of G (={G.domain()})")
        for orbit in G.orbits():
            if not any(set(orbit).issubset(p) for p in pi.values()):
                raise ValueError(f"All elements of orbit {orbit} must have the same sort")
        dis_elm = ConjugacyClassesOfDirectlyIndecomposableSubgroups()(G)
        mapping = {v: i for i, v in enumerate(G.domain(), 1)}
        mapping2 = PermutationGroupElement([mapping[e] for o in sorted(G.orbits(), key=len, reverse=True)
                                            for e in o]).inverse()
        dpart = [frozenset() for _ in range(self._arity)]
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

            sage: A = AtomicSpecies(["X", "Y"])
            sage: A(SymmetricGroup(4), {1: range(1, 5)})  # indirect doctest
            E_4(X)
            sage: A(SymmetricGroup(4), {2: range(1, 5)})
            E_4(Y)
            sage: A(CyclicPermutationGroup(4), {1: range(1, 5)})
            C_4(X)
            sage: A(CyclicPermutationGroup(4), {2: range(1, 5)})
            C_4(Y)
            sage: A(DihedralGroup(4), {1: range(1, 5)})
            P_4(X)
            sage: A(DihedralGroup(4), {2: range(1, 5)})
            P_4(Y)
            sage: A(AlternatingGroup(4), {1: range(1, 5)})
            Eo_4(X)
            sage: A(AlternatingGroup(4), {2: range(1, 5)})
            Eo_4(Y)
        """
        from sage.groups.perm_gps.permgroup import PermutationGroup
        from sage.groups.perm_gps.permgroup_named import (AlternatingGroup,
                                                          CyclicPermutationGroup,
                                                          DihedralGroup,
                                                          SymmetricGroup)

        # prevent infinite recursion in self._element_constructor_
        self._renamed.add(n)
        for i in range(self._arity):
            if n == 1:
                self(SymmetricGroup(1), {i+1: [1]}).rename(self._names[i])

            if self._arity == 1:
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
            if self._arity == 1:
                G = x
                pi = {1: G.domain()}
            else:
                raise ValueError("the assignment of sorts to the domain elements must be provided")
        else:
            G, pi = x
        if not isinstance(G, PermutationGroup_generic):
            raise ValueError(f"{G} must be a permutation group")
        if not set(pi.keys()).issubset(range(1, self._arity + 1)):
            raise ValueError(f"keys of pi (={pi.keys()}) must be in the range [1, {self._arity}]")
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


class MolecularSpecies(IndexedFreeAbelianMonoid):
    """
    The set of (multivariate) molecular species.
    """
    @staticmethod
    def __classcall__(cls, *args, **kwds):
        """
        Normalize the arguments.

        EXAMPLES::

            sage: from sage.rings.species import MolecularSpecies
            sage: MolecularSpecies("X,Y") is MolecularSpecies(["X", "Y"])
            True

            sage: MolecularSpecies("X,Y") == MolecularSpecies(["X", "Z"])
            False
        """
        if isinstance(args[0], AtomicSpecies):
            indices = args[0]
        else:
            assert "names" not in kwds or kwds["names"] is None
            indices = AtomicSpecies(args[0])
        category = Monoids().Commutative() & SetsWithGrading().Infinite()
        return super().__classcall__(cls, indices,
                                     prefix='', bracket=False,
                                     category=category)

    def __init__(self, *args, **kwds):
        r"""
        Initialize the class of (multivariate) molecular species.

        INPUT:

        - ``names`` -- an iterable of ``k`` strings for the sorts of
          the species

        EXAMPLES::

            sage: from sage.rings.species import MolecularSpecies
            sage: M = MolecularSpecies("X,Y")
            sage: G = PermutationGroup([[(1,2),(3,4)], [(5,6)]])
            sage: M(G, {1: [5,6], 2: [1,2,3,4]})
            {((1,2)(3,4),): ({}, {1, 2, 3, 4})}*E_2(X)

        TESTS::

            sage: M1 = MolecularSpecies("X")
            sage: TestSuite(M1).run(skip="_test_graded_components")
            sage: M2 = MolecularSpecies(["X", "Y"])
            sage: TestSuite(M2).run(skip="_test_graded_components")
        """
        IndexedFreeAbelianMonoid.__init__(self, *args, **kwds)
        self._arity = args[0]._arity

    def grading_set(self):
        r"""
        Return the set of non-negative integer vectors, whose length is
        the arity of ``self``.

        EXAMPLES::

            sage: from sage.rings.species import MolecularSpecies
            sage: MolecularSpecies(["X", "Y"]).grading_set()
            Integer vectors of length 2
        """
        return IntegerVectors(length=self._arity)

    def graded_component(self, grade):
        """
        Return the set of molecular species with given multicardinality.

        The default implementation just calls the method :meth:`subset()`
        with the first argument ``grade``.

        EXAMPLES::

            sage: from sage.rings.species import MolecularSpecies
            sage: M = MolecularSpecies(["X", "Y"])
            sage: M.graded_component([3,2])  # random
            {E_3(X)*Y^2, X^3*Y^2, X*E_2(X)*E_2(Y), X^3*E_2(Y),
             {((1,2,3), (1,3)(4,5)): ({1, 2, 3}, {4, 5})},
             X*{((1,2)(3,4),): ({1, 2}, {3, 4})}, X*E_2(X)*Y^2, E_3(X)*E_2(Y),
             C_3(X)*Y^2, C_3(X)*E_2(Y)}
        """
        from sage.sets.set import Set
        assert len(grade) == self._arity
        n = sum(grade)
        S = SymmetricGroup(n).young_subgroup(grade)
        dom = S.domain()
        dom_part = {i+1: dom[sum(grade[:i]): sum(grade[:i+1])] for i in range(len(grade))}
        return Set([self(G, dom_part) for G in S.conjugacy_classes_subgroups()])

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

        EXAMPLES:

        Create a molecular species given a permutation group::

            sage: from sage.rings.species import MolecularSpecies
            sage: M = MolecularSpecies(["X", "Y"])
            sage: G = PermutationGroup([[(1,2)], [(3,4)]])
            sage: M(G, {1: [1,2], 2: [3,4]})
            E_2(X)*E_2(Y)

        TESTS::

            sage: M(CyclicPermutationGroup(4), {1: [1,2], 2: [3,4]})
            Traceback (most recent call last):
            ...
            ValueError: All elements of orbit (1, 2, 3, 4) must have the same sort

            sage: G = PermutationGroup([[(2,3),(4,5)]], domain=[2,3,4,5])
            sage: M(G, {1:[2, 3], 2:[4,5]})
            {((1,2)(3,4),): ({1, 2}, {3, 4})}
        """
        if parent(G) == self:
            if pi is not None:
                raise ValueError("cannot reassign sorts to a molecular species")
            return G

        if isinstance(G, PermutationGroup_generic):
            if pi is None:
                if self._arity == 1:
                    pi = {1: G.domain()}
                else:
                    raise ValueError("the assignment of sorts to the domain elements must be provided")
            domain = [e for p in pi.values() for e in p]
            if len(domain) != len(set(domain)):
                raise ValueError("each domain element must have exactly one sort")
            if set(G.domain()) != set(domain):
                raise ValueError("each element of the domain of the group must have one sort")
            domain_partition = G.disjoint_direct_product_decomposition()
            elm = self.one()
            for part in domain_partition:
                elm *= self.gen(self._project(G, pi, part))
            return elm

        # Assume G is a tuple (X, a)
        X, a = G
        if pi is None:
            if self._arity == 1:
                pi = {1: X}
            else:
                raise ValueError("the assignment of sorts to the domain elements must be provided")
        L = [None for _ in range(self._arity)]
        for k, v in pi.items():
            L[k - 1] = list(v)
        # Create group
        # TODO: Is this correct?
        S = SymmetricGroup(list(chain.from_iterable(L))).young_subgroup([len(v) for v in L])
        H = PermutationGroup(S.gens(), action=a, domain=X)
        if len(H.orbits()) > 1:
            # Then it is not transitive
            raise ValueError("Action is not transitive")
        stabG = PermutationGroup([g for g in S.gens()
                                  if a(g, H.orbits()[0][0]) == H.orbits()[0][0]],
                                 domain=S.domain())
        return self(stabG, pi)

    def gen(self, x):
        r"""
        The molecular species given by an atomic species.

        EXAMPLES::

            sage: from sage.rings.species import MolecularSpecies
            sage: M = MolecularSpecies("X, Y")

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

            sage: P = PolynomialSpecies(ZZ, ["X"])
            sage: M = P._indices
            sage: A = AtomicSpecies("X")
            sage: a = A(CyclicPermutationGroup(4))
            sage: M.gen(a)
            C_4

            sage: P = PolynomialSpecies(ZZ, ["X", "Y"])
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
        return self.element_class(self, {at: ZZ.one()})

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        TESTS::

            sage: from sage.rings.species import MolecularSpecies
            sage: MolecularSpecies("X")
            Molecular species in X
            sage: MolecularSpecies("A, B")
            Molecular species in A, B
        """
        if len(self._indices._names) == 1:
            return f"Molecular species in {self._indices._names[0]}"
        return f"Molecular species in {', '.join(self._indices._names)}"

    class Element(IndexedFreeAbelianMonoidElement):
        def __init__(self, parent, x):
            r"""
            Initialize a molecular species.

            INPUT:

            - ``x``, a dictionary mapping atomic species to exponents

            EXAMPLES::

                sage: from sage.rings.species import MolecularSpecies
                sage: M = MolecularSpecies("X")
                sage: M(CyclicPermutationGroup(3))  # indirect doctest
                C_3
            """
            super().__init__(parent, x)

        @cached_method
        def group_and_partition(self):
            """
            Return the (transitive) permutation group corresponding to ``self``.

            EXAMPLES::

                sage: from sage.rings.species import MolecularSpecies
                sage: M = MolecularSpecies("X,Y")
                sage: G = PermutationGroup([[(1,2),(3,4)], [(5,6)]])
                sage: A = M(G, {1: [5,6], 2: [1,2,3,4]})
                sage: A.group_and_partition()
                (Permutation Group with generators [(5,6), (1,2)(3,4)],
                 (frozenset({5, 6}), frozenset({1, 2, 3, 4})))

            TESTS::

                sage: B = M(PermutationGroup([(1,2,3)]), {1: [1,2,3]})
                sage: B.group_and_partition()
                (Permutation Group with generators [(1,2,3)],
                 (frozenset({1, 2, 3}), frozenset()))

                sage: (A*B).group_and_partition()
                (Permutation Group with generators [(7,8,9), (5,6), (1,2)(3,4)],
                 (frozenset({5, 6, 7, 8, 9}), frozenset({1, 2, 3, 4})))

                sage: C = M(PermutationGroup([(2,3)]), {1: [1], 2: [2,3]})
                sage: C.group_and_partition()
                (Permutation Group with generators [(2,3)],
                 (frozenset({1}), frozenset({2, 3})))

                sage: (C^3).group_and_partition()
                (Permutation Group with generators [(8,9), (6,7), (4,5)],
                 (frozenset({1, 2, 3}), frozenset({4, 5, 6, 7, 8, 9})))

                sage: M = MolecularSpecies("X")
                sage: F = M(SymmetricGroup(1)) * M(SymmetricGroup(2))
                sage: F.group_and_partition()
                (Permutation Group with generators [(2,3)], (frozenset({1, 2, 3}),))
            """
            def shift_gens(gens, n):
                """
                Given a list of generators ``gens``, increase every element of the
                domain by ``n``.
                """
                return tuple([tuple([tuple([n + e for e in cyc])
                                     for cyc in gen.cycle_tuples()])
                              for gen in gens])

            factors = list(self)
            if not factors:
                k = self.parent()._arity
                return SymmetricGroup(0), tuple([frozenset()]*k)

            if len(factors) == 1:
                A, n = factors[0]
                if n == 1:
                    a = list(A._monomial)[0]  # as atomic species
                    return a._dis._C, a._dompart

                if n % 2 == 1:
                    a = list(A._monomial)[0] # as atomic species
                    b, b_dompart = (A ** (n-1)).group_and_partition()
                    gens = a._dis._C.gens() + shift_gens(b.gens(), a._tc)
                    new_dompart = tuple([frozenset(list(p_a) + [a._tc + e for e in p_b])
                                         for p_a, p_b in zip(a._dompart, b_dompart)])
                else:
                    f, f_dompart = (A ** (n // 2)).group_and_partition()
                    tc = sum(len(p) for p in f_dompart)
                    gens = f.gens() + shift_gens(f.gens(), tc)
                    new_dompart = tuple([frozenset(list(p) + [tc + e for e in p])
                                         for p in f_dompart])

                G = PermutationGroup(gens)
                return G, new_dompart

            f_dompart_list = [(A ** n).group_and_partition() for A, n in factors]
            f_list = [f for f, _ in f_dompart_list]
            dompart_list = [f_dompart for _, f_dompart in f_dompart_list]
            tc_list = list(accumulate([sum(len(p) for p in f_dompart)
                                       for f_dompart in dompart_list],
                                      initial=0))
            gens = [gen
                    for f, tc in zip(f_list, tc_list)
                    for gen in shift_gens(f.gens(), tc) if gen]  # gen is a tuple
            G = PermutationGroup(gens)
            new_dompart = tuple([frozenset(chain(*[[tc + e for e in p]
                                                   for p, tc in zip(f_dompart, tc_list)]))
                                 for f_dompart in zip(*dompart_list)])

            return G, new_dompart

        @cached_method
        def grade(self):
            r"""
            Return the grade of ``self``.

            EXAMPLES::

                sage: from sage.rings.species import MolecularSpecies
                sage: M = MolecularSpecies("X, Y")
                sage: G = PermutationGroup([[(1,2),(3,4),(5,6),(7,8,9,10)]]); G
                Permutation Group with generators [(1,2)(3,4)(5,6)(7,8,9,10)]
                sage: a = M(G, {1: [1,2,3,4], 2: [5,6,7,8,9,10]}); a
                {((1,2,3,4)(5,6)(7,8)(9,10),): ({5, 6, 7, 8}, {1, 2, 3, 4, 9, 10})}
                sage: a.grade()
                [4, 6]

            TESTS::

                sage: M.one().grade()
                [0, 0]
            """
            P = self.parent()
            S = P.grading_set()
            mons = self._monomial
            if not mons:
                return S([0] * P._arity)
            mc = sum(n * vector(a._mc) for a, n in mons.items())
            return S(mc)

        def _compose_with_singletons(self, base_ring, names, args):
            r"""
            Return the inner sum of Exercise 2.6.16 in [BLL1998]_,
            generalized to the case of arbitrary arity.

            INPUT:

                - ``base_ring``, the base ring of the result

                - ``names``, the (flat) list of names of the result

                - ``args``, the sequence of `k` compositions, each of
                  which sums to the corresponding cardinality of
                  ``self``, where `k` is the arity of ``self``.

            OUTPUT:

            - the polynomial species
              `self(X_1 + \dots + X_{m_1}, Y_1 + \dots + Y_{m_2}, \dots)`,
              where `m_i` is the number of parts of the `i`-th composition.

            EXAMPLES::

                sage: from sage.rings.species import MolecularSpecies
                sage: M = MolecularSpecies("X")
                sage: C4 = M(CyclicPermutationGroup(4))
                sage: C4._compose_with_singletons(ZZ, "X, Y", [[2, 2]]) # X^2Y^2 + C2(XY)
                E_2(XY) + X^2*Y^2

                sage: M = MolecularSpecies(["X", "Y"])
                sage: F = M(PermutationGroup([[(1,2,3), (4,5,6)]]), {1: [1,2,3], 2: [4,5,6]})
                sage: F
                {((1,2,3)(4,5,6),): ({1, 2, 3}, {4, 5, 6})}
                sage: F._compose_with_singletons(ZZ, "X1, X2, X3, Y1, Y2", [[1, 1, 1], [2, 1]])
                6*X1*X2*X3*Y1^2*Y2

            TESTS::

                sage: M = MolecularSpecies("X")
                sage: O = M.one()
                sage: O._compose_with_singletons(ZZ, "X", [[]])
                1

                sage: F = M(SymmetricGroup(1)) * M(SymmetricGroup(2))
                sage: F._compose_with_singletons(QQ, ["T", "S"], [[2, 1]])
                T^2*S + E_2(T)*S
                sage: F._compose_with_singletons(QQ, ["T", "S"], [[1, 2]])
                T*E_2(S) + T*S^2
            """
            # TODO: No checks are performed right now, must be added.
            # Checks: all args in Compositions, sums must match cardinalities.

            # Create group of the composition
            # conjugate self.group() so that [1..k] is sort 1, [k+1,..] is sort 2, so on
            G, dompart = self.group_and_partition()
            conj = PermutationGroupElement(list(chain.from_iterable(dompart))).inverse()
            G = libgap.ConjugateGroup(G, conj)

            comp = list(chain.from_iterable(args))
            dpart = {i + 1: range(x - comp[i] + 1, x + 1) for i, x in enumerate(accumulate(comp))}
            # Create the double coset representatives.
            S_down = SymmetricGroup(sum(comp)).young_subgroup(comp)
            mc = self.grade()
            tc = sum(mc)
            S_up = SymmetricGroup(tc).young_subgroup(mc)
            taus = libgap.DoubleCosetRepsAndSizes(S_up, S_down, G)
            # sum over double coset representatives.
            Pn = PolynomialSpecies(base_ring, names)
            res = Pn.zero()
            for tau, _ in taus:
                H = libgap.Intersection(libgap.ConjugateGroup(G, tau ** -1), S_down)
                grp = PermutationGroup(gap_group=H, domain=range(1, tc + 1))
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

                sage: from sage.rings.species import MolecularSpecies
                sage: M = MolecularSpecies("X")
                sage: X = M(SymmetricGroup(1))
                sage: E2 = M(SymmetricGroup(2))
                sage: E2(X)
                E_2
                sage: X(E2)
                E_2
                sage: E2(E2)
                P_4

                sage: M = MolecularSpecies(["X","Y"])
                sage: X = M(SymmetricGroup(1), {1:[1]})
                sage: Y = M(SymmetricGroup(1), {2:[1]})
                sage: (X*Y)(X, Y^2)
                X*Y^2

            A multivariate example::

                sage: M1 = MolecularSpecies("X")
                sage: M2 = MolecularSpecies("X, Y")
                sage: C3 = M1(CyclicPermutationGroup(3))
                sage: X = M2(SymmetricGroup(1), {1: [1]})
                sage: Y = M2(SymmetricGroup(1), {2: [1]})
                sage: C3(X*Y)
                {((1,2,3)(4,5,6),): ({1, 2, 3}, {4, 5, 6})}

            TESTS::

                sage: M = MolecularSpecies("X")
                sage: M.one()()
                Traceback (most recent call last):
                ...
                ValueError: number of args must match arity of self
            """
            if len(args) != self.parent()._arity:
                raise ValueError("number of args must match arity of self")
            if len(set(arg.parent() for arg in args)) > 1:
                raise ValueError("all args must have the same parent")
            if not all(isinstance(arg, MolecularSpecies.Element) for arg in args):
                raise ValueError("all args must be molecular species")

            # TODO: What happens if in F(G), G has a constant part? E(1+X)?
            Mlist = [None for _ in range(sum(self.grade()))]
            G, dompart = self.group_and_partition()
            for i, v in enumerate(dompart):
                for k in v:
                    Mlist[k - 1] = args[i]
            starts = list(accumulate([sum(M.grade()) for M in Mlist],
                                     initial=0))

            # gens from self
            gens = []
            for gen in G.gens():
                newgen = []
                for cyc in gen.cycle_tuples():
                    for k in range(1, sum(Mlist[cyc[0] - 1].grade()) + 1):
                        newgen.append(tuple(k + starts[i - 1] for i in cyc))
                gens.append(newgen)

            # gens from M_i and dompart
            P = args[0].parent()
            dpart = {i: [] for i in range(1, P._arity + 1)}
            for start, M in zip(starts, Mlist):
                K, dompart = M.group_and_partition()
                for i, v in enumerate(dompart, 1):
                    dpart[i].extend([start + k for k in v])
                for gen in K.gens():
                    gens.append([tuple(start + k for k in cyc)
                                 for cyc in gen.cycle_tuples()])

            return P(PermutationGroup(gens,
                                      domain=range(1, starts[-1] + 1)), dpart)


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
        category = GradedAlgebrasWithBasis(base_ring).Commutative()
        CombinatorialFreeModule.__init__(self, base_ring,
                                         basis_keys=MolecularSpecies(names),
                                         category=category,
                                         element_class=self.Element,
                                         prefix='', bracket=False)
        self._arity = len(names)

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
        return sum(m.grade())

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
                if self._arity == 1:
                    return self._from_dict({self._indices(G): ZZ.one()})
                raise ValueError("the assignment of sorts to the domain elements must be provided")
            return self._from_dict({self._indices(G, pi): ZZ.one()})
        # Assume G is a tuple (X, a)
        X, a = G
        if pi is None:
            if self._arity == 1:
                pi = {1: X}
            else:
                raise ValueError("the assignment of sorts to the domain elements must be provided")
        # Make iteration over values of pi deterministic
        L = [None for _ in range(self._arity)]
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
            sage: matrix([[F * G for F in L1] for G in L2])  # indirect doctest
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
                         for s in range(self._arity))

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

            Exercise 2.1.9 from [BLL1998]_::

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
                mc = L.grade()
                tc = sum(mc)
                S = SymmetricGroup(tc).young_subgroup(mc)
                # conjugate L and R to match S
                G, dompart = L.group_and_partition()
                g = list(chain.from_iterable(dompart))
                conj_L = PermutationGroupElement(g).inverse()
                G = libgap.ConjugateGroup(G, conj_L)
                dpart = {i + 1: range(x - mc[i] + 1, x + 1)
                         for i, x in enumerate(accumulate(mc))}
                for R, d in other:
                    if mc != R.grade():
                        continue
                    G_R, dompart_R = R.group_and_partition()
                    g = list(chain.from_iterable(dompart_R))
                    conj_R = PermutationGroupElement(g).inverse()
                    H = libgap.ConjugateGroup(G_R, conj_R)
                    taus = libgap.DoubleCosetRepsAndSizes(S, G, H)
                    # loop over representatives
                    new = P.zero()
                    for tau, _ in taus:
                        F = libgap.Intersection(libgap.ConjugateGroup(H, tau), G)
                        new += P(PermutationGroup(gap_group=F, domain=range(1, tc + 1)),
                                 dpart)
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

            Equation (2.5.41) in [BLL1998]_::

                sage: P = PolynomialSpecies(QQ, ["X"])
                sage: E2 = P(SymmetricGroup(2))
                sage: E2._compose_with_weighted_singletons(["X"], [-1], [[2]])
                -E_2 + X^2

                sage: C4 = P(CyclicPermutationGroup(4))
                sage: C4._compose_with_weighted_singletons(["X"], [-1], [[4]])
                -C_4 + {((1,2)(3,4),)}

            Exercise (2.5.17) in [BLL1998]_::

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
                -C_4 + {((1,2)(3,4),)} + E_2^2 - 2*X^2*E_2 + X^4

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
                {((1,2)(3,4),)}

                sage: E2(X + X^2)
                E_2 + X^3 + {((1,2)(3,4),)}

                sage: P2 = PolynomialSpecies(QQ, ["X", "Y"])
                sage: X = P2(SymmetricGroup(1), {1:[1]})
                sage: Y = P2(SymmetricGroup(1), {2:[1]})
                sage: E2(X + Y)
                E_2(Y) + Y*X + E_2(X)

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
            if len(args) != P._arity:
                raise ValueError("number of args must match arity of self")
            if len(set(arg.parent() for arg in args)) > 1:
                raise ValueError("all args must have the same parent")

            P0 = args[0].parent()
            if not self.support():
                return P0.zero()

            args = [sorted(g, key=lambda x: x[0].grade()) for g in args]
            multiplicities = list(chain.from_iterable([[c for _, c in g] for g in args]))
            molecules = list(chain.from_iterable([[M for M, _ in g] for g in args]))
            F_degrees = sorted(set(M.grade() for M, _ in self))
            names = ["X%s" % i for i in range(sum(len(arg) for arg in args))]

            result = P0.zero()
            for n in F_degrees:
                F = P.sum_of_terms((M, c) for M, c in self if M.grade() == n)
                for degrees in cartesian_product([IntegerVectors(n_i, length=len(arg))
                                                  for n_i, arg in zip(n, args)]):
                    # each degree is a weak composition of the degree of F in sort i
                    FX = F._compose_with_weighted_singletons(names,
                                                             multiplicities,
                                                             degrees)
                    FG = [(M(*molecules), c) for M, c in FX]
                    result += P0.sum_of_terms(FG)
            return result
