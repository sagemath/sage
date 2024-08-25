from itertools import accumulate, chain

from sage.categories.graded_algebras_with_basis import GradedAlgebrasWithBasis
from sage.categories.infinite_enumerated_sets import InfiniteEnumeratedSets
from sage.categories.monoids import Monoids
from sage.categories.sets_with_grading import SetsWithGrading
from sage.categories.cartesian_product import cartesian_product
from sage.combinat.free_module import CombinatorialFreeModule
from sage.combinat.integer_vector import IntegerVectors
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
        Additionally, if a method ``_canonicalize`` is implemented, it is used to preprocess the element.
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
        sorted_orbits = sorted([sorted(orbit) for orbit in self._C.orbits()], key=len, reverse=True)
        pi = PermutationGroupElement(list(chain.from_iterable(sorted_orbits))).inverse()
        self._C = PermutationGroup(gap_group=libgap.ConjugateGroup(self._C, pi))

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

    def _element_key(self):
        r"""
        Return the cache lookup key for ``self``.
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
            sage: H = PermutationGroup([[(1,2,3,4),(5,6),(7,8),(9,10)]]); H
            Permutation Group with generators [(1,2,3,4)(5,6)(7,8)(9,10)]
            sage: A = At(G, {1: [1,2,3,4], 2: [5,6,7,8,9,10]})
            sage: A._dompart
            ((5, 6, 7, 8), (9, 10, 1, 2, 3, 4))
            sage: C = At(G, {1: [1,2,5,6], 2: [3,4,7,8,9,10]})
            sage: C._dompart
            ((5, 6, 7, 8), (9, 10, 1, 2, 3, 4))
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
            {((1,2,3,4)(5,6)(7,8)(9,10),): ((5, 6, 7, 8), (9, 10, 1, 2, 3, 4))}
            sage: B = At(H, {1: [1,2,3,4], 2: [5,6,7,8,9,10]}); B
            {((1,2,3,4)(5,6)(7,8)(9,10),): ((1, 2, 3, 4), (5, 6, 7, 8, 9, 10))}
            sage: C = At(G, {1: [1,2,5,6], 2: [3,4,7,8,9,10]}); C
            {((1,2,3,4)(5,6)(7,8)(9,10),): ((5, 6, 7, 8), (9, 10, 1, 2, 3, 4))}
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
            selflist.extend(self._dompart[i])
            otherlist.extend(other._dompart[i])
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
            {((1,2,3,4)(5,6)(7,8)(9,10),): ((5, 6, 7, 8), (9, 10, 1, 2, 3, 4))}
        """
        return "{" + f"{self._dis}: {self._dompart}" + "}"


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
            {((1,2)(3,4),): ((1, 2), (3, 4))}
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
        dpart = [tuple() for _ in range(self._k)]
        for k, v in pi.items():
            dpart[k - 1] = tuple(mapping2(mapping[x]) for x in v)
        elm = self._cache_get(self.element_class(self, dis_elm, tuple(dpart)))
        if elm._tc not in self._renamed:
            self._rename(elm._tc)
        return elm

    def _rename(self, n):
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
    def __classcall__(cls, indices, prefix, **kwds):
        return super(IndexedMonoid, cls).__classcall__(cls, indices, prefix, **kwds)

    def __init__(self, indices, prefix, **kwds):
        category = Monoids() & InfiniteEnumeratedSets()
        IndexedFreeAbelianMonoid.__init__(self, indices, prefix=prefix, category=category, **kwds)
        ElementCache.__init__(self)
        self._k = indices._k

    def _project(self, G, pi, part):
        r"""
        Project `G` onto a subset ``part`` of its domain.

        ``part`` must be a union of cycles, but this is not checked.
        """
        restricted_gens = [[cyc for cyc in gen.cycle_tuples() if cyc[0] in part] for gen in G.gens()]
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
        """
        if x not in self._indices:
            raise IndexError(f"{x} is not in the index set")
        at = None
        if isinstance(x, PermutationGroup_generic):
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

    class Element(IndexedFreeAbelianMonoidElement):
        def __init__(self, F, x):
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
            for gen in gens2:
                gens.append([tuple(elm1._tc + k for k in cyc) for cyc in gen.cycle_tuples()])
            self._group = PermutationGroup(gens, domain=range(1, elm1._tc + elm2._tc + 1))
            self._dompart = list(elm1._dompart)
            for i in range(elm2.parent()._k):
                self._dompart[i] = tuple(list(self._dompart[i]) + [elm1._tc + e for e in elm2._dompart[i]])
            self._dompart = tuple(self._dompart)

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
            """
            return self

        @cached_method
        def _canonicalize(self):
            r"""
            Canonicalize this molecular species by sorting the orbits by
            length and making them consecutive.

            EXAMPLES::

                sage: P = PolynomialSpecies(ZZ, "X, Y")
                sage: G = PermutationGroup([[(1,2),(3,4),(5,6),(7,8,9,10)]]); G
                Permutation Group with generators [(1,2)(3,4)(5,6)(7,8,9,10)]
                sage: H = PermutationGroup([[(1,2,3,4),(5,6),(7,8),(9,10)]]); H
                Permutation Group with generators [(1,2,3,4)(5,6)(7,8)(9,10)]
                sage: A = P(G, {1: [1,2,3,4], 2: [5,6,7,8,9,10]})
                sage: A.support()[0]._dompart
                ((5, 6, 7, 8), (9, 10, 1, 2, 3, 4))
                sage: C = P(G, {1: [1,2,5,6], 2: [3,4,7,8,9,10]})
                sage: C.support()[0]._dompart
                ((5, 6, 7, 8), (9, 10, 1, 2, 3, 4))
            """
            if self._group is None or self._group == SymmetricGroup(0):
                return
            sorted_orbits = sorted([sorted(orbit) for orbit in self._group.orbits()], key=len, reverse=True)
            pi = PermutationGroupElement(list(chain.from_iterable(sorted_orbits))).inverse()
            self._group = PermutationGroup(gap_group=libgap.ConjugateGroup(self._group, pi))
            self._dompart = tuple(tuple(pi(k) for k in v) for v in self._dompart)

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

        def hadamard_product(self, other):
            r"""
            Compute the hadamard product of ``self`` and ``other``.

            EXAMPLES:

            Exercise 2.1.9 from the BLL book::

                sage: P = PolynomialSpecies(ZZ, ["X"])
                sage: M = P._indices
                sage: C3 = M(CyclicPermutationGroup(3))
                sage: X = M(SymmetricGroup(1))
                sage: E2 = M(SymmetricGroup(2))
                sage: C3.hadamard_product(C3)
                2*C_3
                sage: (X^3).hadamard_product(C3)
                2*X^3
                sage: (X*E2).hadamard_product(X*E2)
                X*E_2 + X^3
            """
            P = self.parent()
            if P is not other.parent():
                raise ValueError("the factors of a Hadamard product must be the same")
            Pn = PolynomialSpecies(ZZ, P._indices._names)

            if self._mc != other._mc:
                return Pn.zero()
            # create S
            S = SymmetricGroup(self._tc).young_subgroup(self._mc)
            # conjugate self and other to match S
            conj_self = PermutationGroupElement(list(chain.from_iterable(self._dompart))).inverse()
            conj_other = PermutationGroupElement(list(chain.from_iterable(other._dompart))).inverse()
            G = libgap.ConjugateGroup(self._group, conj_self)
            H = libgap.ConjugateGroup(other._group, conj_other)
            # create dompart
            dpart = {i + 1: range(x - self._mc[i] + 1, x + 1) for i, x in enumerate(accumulate(self._mc))}
            # create double coset representatives
            taus = libgap.DoubleCosetRepsAndSizes(S, G, H)
            # loop over representatives
            res = Pn.zero()
            for tau, _ in taus:
                F = libgap.Intersection(libgap.ConjugateGroup(H, tau), G)
                res += Pn(PermutationGroup(gap_group=F, domain=self.domain()), dpart)
            return res

        def functorial_composition(self, other):
            r"""
            Return the functorial composition of ``self`` with ``other``.
            Currrent this only works for univariate molecular species only.

            TESTS::

                sage: P = PolynomialSpecies(ZZ, ["X"])
                sage: M = P._indices
                sage: E = [None] + [M(SymmetricGroup(n)) for n in range(1, 4)]
                sage: E[2].functorial_composition(E[2])
                E_4
                sage: E[1].functorial_composition(E[1])
                X
                sage: E[1].functorial_composition(E[2])
                E_2
                sage: E[2].functorial_composition(E[1])
                X
                sage: E[3].functorial_composition(E[2])
                E_8
                sage: E[2].functorial_composition(E[3])
                E_9
            """
            if not isinstance(other, MolecularSpecies.Element):
                raise ValueError(f"{other} must be a molecular species")
            if other.parent()._k > 1:
                raise ValueError(f"{other} must be univariate")
            if self.parent()._k > 1:
                raise ValueError(f"{self} must be univariate")

            gens = []
            dpart = {1: range(1, other.grade() ** self.grade() + 1)}

            def to_number(l):
                B = ZZ.one() # ZZ.one or 1?
                R = 0 # Just 0 is fine?
                for e in reversed(l):
                    R += e * B
                    B *= other.grade()
                return R

            # it is a base other.grade() representation
            for pre in cartesian_product([range(other.grade())]*self.grade()):
                pre_n = to_number(pre) + 1
                # gens from self
                for gen in self._group.gens():
                    # gen swaps around the numbers in pre
                    im = gen(list(pre))
                    im_n = to_number(im) + 1
                    if pre_n == im_n:
                        continue
                    gens.append([tuple([pre_n, im_n])])
            print('after self',gens)

            # gens from other
            for gen in other._group.gens():
                for pre in cartesian_product([range(other.grade())]*self.grade()):
                    pre_n = to_number(pre) + 1
                    for i in range(self.grade()):
                        im_n = pre_n - pre[-i-1] * other.grade() ** i + (gen(pre[-i-1] + 1) - 1) * other.grade() ** i
                        if pre_n == im_n:
                            continue
                        gens.append([tuple([pre_n, im_n])])
            print('after others',gens)

            G = PermutationGroup(gens,domain=range(1, other.grade() ** self.grade() + 1))
            print('G',G)
            print('dompart',dpart)
            return other.parent()(G, dpart)

        def inner_sum(self, base_ring, names, *args):
            r"""
            Compute the inner sum of exercise 2.6.16 of BLL book.

            INPUT:

                - ``base_ring``, the base ring of the result

                - ``names``, the names of the result

                - ``args``, the sequence of compositions, each of
                  which sums to the corresponding cardinality of
                  ``self``. The number of ``args`` is equal to the
                  arity of ``self``.

            EXAMPLES::

                sage: P = PolynomialSpecies(ZZ, "X")
                sage: M = P._indices
                sage: C4 = M(CyclicPermutationGroup(4))
                sage: C4.inner_sum(ZZ, "X, Y", [2, 2]) # X^2Y^2 + C2(XY)
                {((1,2)(3,4),): ((1, 2), (3, 4))} + X^2*Y^2

            """
            # TODO: No checks are performed right now, must be added.
            # Checks: all args in compositions, sums must match cardinalities.

            # Create group of the composition
            Pn = PolynomialSpecies(base_ring, names)
            res = Pn.zero()
            # conjugate self._group so that [1..k] is sort 1, [k+1,..] is sort 2, so on
            conj = PermutationGroupElement(list(chain.from_iterable(self._dompart))).inverse()
            G = libgap.ConjugateGroup(self._group, conj)

            comp = list(chain.from_iterable(args))
            # Create the double coset representatives.
            S_down = SymmetricGroup(sum(comp)).young_subgroup(comp)
            S_up = SymmetricGroup(self._tc).young_subgroup(self._mc)
            taus = libgap.DoubleCosetRepsAndSizes(S_up, S_down, G)
            # Sum over double coset representatives.
            for tau, _ in taus:
                H = libgap.Intersection(libgap.ConjugateGroup(G, tau), S_down)
                grp = PermutationGroup(gap_group=H, domain=self.domain())
                dpart = {i + 1: range(x - comp[i] + 1, x + 1) for i, x in enumerate(accumulate(comp))}
                res += Pn(grp, dpart)
            return res

        def __call__(self, *args):
            r"""
            Substitute M_1...M_k into self.
            M_i must all have same arity, same multicardinality,
            and must be molecular.

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
            """
            if len(args) != self.parent()._k:
                raise ValueError("number of args must match arity of self")
            if not all(isinstance(arg, MolecularSpecies.Element) for arg in args):
                raise ValueError("all args must be molecular species")
            if len(set(arg.parent()._k for arg in args)) > 1:
                raise ValueError("all args must have same arity")
            if len(set(arg._mc for arg in args)) > 1:
                raise ValueError("all args must have same multicardinality")

            gens = []

            # TODO: What happens if in F(G), G has a constant part? E(1+X)?
            Mlist = [None for _ in range(self.grade())]
            for i, v in enumerate(self._dompart):
                for k in v:
                    Mlist[k - 1] = args[i]
            starts = list(accumulate([M.grade() for M in Mlist], initial=0))

            # gens from self
            for gen in self._group.gens():
                newgen = []
                for cyc in gen.cycle_tuples():
                    for k in range(1, Mlist[cyc[0] - 1].grade() + 1):
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
        """
        return m.grade()

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
        """
        return self._from_dict({H * K: 1})

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
            """
            return self.is_zero() or not self.degree()
        
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
