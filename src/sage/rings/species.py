from itertools import accumulate, chain

from sage.categories.graded_algebras_with_basis import GradedAlgebrasWithBasis
from sage.categories.infinite_enumerated_sets import InfiniteEnumeratedSets
from sage.categories.modules import Modules
from sage.categories.monoids import Monoids
from sage.categories.sets_cat import cartesian_product
from sage.categories.sets_with_grading import SetsWithGrading
from sage.categories.tensor import tensor
from sage.combinat.free_module import CombinatorialFreeModule
from sage.combinat.integer_vector import IntegerVectors
from sage.combinat.partition import Partitions, _Partitions
from sage.combinat.sf.sf import SymmetricFunctions
from sage.groups.perm_gps.constructor import PermutationGroupElement
from sage.groups.perm_gps.permgroup import PermutationGroup, PermutationGroup_generic
from sage.groups.perm_gps.permgroup_named import SymmetricGroup
from sage.libs.gap.libgap import libgap
from sage.misc.cachefunc import cached_method
from sage.misc.fast_methods import WithEqualityById
from sage.misc.inherit_comparison import InheritComparisonClasscallMetaclass
from sage.modules.free_module_element import vector
from sage.monoids.indexed_free_monoid import (IndexedFreeAbelianMonoid,
                                              IndexedFreeAbelianMonoidElement)
from sage.rings.rational_field import QQ
from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.sets.set import Set
from sage.structure.category_object import normalize_names
from sage.structure.element import Element, parent
from sage.structure.factorization import Factorization
from sage.structure.parent import Parent
from sage.structure.richcmp import op_LT, op_LE, op_EQ, op_NE, op_GT, op_GE
from sage.structure.unique_representation import (UniqueRepresentation,
                                                  WithPicklingByInitArgs)

GAP_FAIL = libgap.eval('fail')


class ConjugacyClassOfDirectlyIndecomposableSubgroups(Element, WithEqualityById,
                                                      WithPicklingByInitArgs,
                                                      metaclass=InheritComparisonClasscallMetaclass):
    r"""
    A directly indecomposable conjugacy class of subgroups of a
    symmetric group.

    Two conjugacy classes of subgroups are equal if they have the
    same degree (say `n`) and are conjugate within `S_n`.
    """
    @staticmethod
    def _invariant(C):
        r"""
        Return an invariant of the permutation group ``C`` used for
        hashing and as key in the cache.

        EXAMPLES::

            sage: from sage.rings.species import ConjugacyClassOfDirectlyIndecomposableSubgroups
            sage: C = ConjugacyClassOfDirectlyIndecomposableSubgroups
            sage: G = PermutationGroup([[(1,2),(3,4)]])
            sage: H = PermutationGroup([[(1,3),(2,4)], [(1,2),(3,4)]])
            sage: K = PermutationGroup([[(1,2,3,4)]])
            sage: C._invariant(G)
            (2, (2, 2))
            sage: C._invariant(H)
            (4, (4,))
            sage: C._invariant(K)
            (4, (4,))
        """
        return C.cardinality(), tuple(sorted(len(o) for o in C.orbits()))

    @staticmethod
    def __classcall__(cls, parent, C):
        """
        Normalize the input for unique representation.

        TESTS::

            sage: from sage.rings.species import ConjugacyClassesOfDirectlyIndecomposableSubgroups
            sage: C = ConjugacyClassesOfDirectlyIndecomposableSubgroups()

        Check that the domain is standardized::

            sage: G = PermutationGroup([("a", "b", "c"), ("d", "a"), ("d", "b"), ("d", "c")])
            sage: C(G)
            ((2,4,3), (1,2))

            sage: G = PermutationGroup([[(1,2),(3,4),(5,6),(7,8,9,10)]]); G
            Permutation Group with generators [(1,2)(3,4)(5,6)(7,8,9,10)]
            sage: H = PermutationGroup([[(1,2,3,4),(5,6),(7,8),(9,10)]]); H
            Permutation Group with generators [(1,2,3,4)(5,6)(7,8)(9,10)]
            sage: C(G) is C(H)
            True

        Two different transitive groups with the same cardinality::

            sage: a = C(PermutationGroup([[(1,3),(2,4)], [(1,4),(2,3)]]))
            sage: b = C(PermutationGroup([[(1,3,2,4)]]))
            sage: a == b
            False
        """
        def is_equal(elm):
            return SymmetricGroup(elm._C.degree()).are_conjugate(elm._C, C)

        def standardize(C):
            if not C.degree():
                return SymmetricGroup(0)
            # a directly indecomposable group is transitive, so we
            # can use the gap group without worrying about the domain
            G = C.gap()
            sorted_orbits = sorted(G.Orbits().sage(), key=len, reverse=True)
            pi = PermutationGroupElement([e for o in sorted_orbits for e in o])
            G = PermutationGroup(G.SmallGeneratingSet().sage())
            return G.conjugate(pi.inverse())

        key = cls._invariant(C)
        if key in parent._cache:
            lookup = parent._cache[key]
            for elm in lookup:
                if is_equal(elm):
                    return elm
        else:
            lookup = parent._cache[key] = []

        elm = WithPicklingByInitArgs.__classcall__(cls, parent, standardize(C))
        lookup.append(elm)
        return elm

    def __hash__(self):
        """
        Return a hash of ``self``.

        EXAMPLES::

            sage: from sage.rings.species import ConjugacyClassesOfDirectlyIndecomposableSubgroups
            sage: C = ConjugacyClassesOfDirectlyIndecomposableSubgroups()
            sage: G = PermutationGroup([[(1,2),(3,4)]])
            sage: H = PermutationGroup([[(1,3),(2,4)], [(1,2),(3,4)]])
            sage: K = PermutationGroup([[(1,2,3,4)]])
            sage: hash(C(G)) == hash(C(H))
            False
            sage: hash(C(H)) == hash(C(K))
            True
            sage: C(H) == C(K)
            False
        """
        return hash(self._invariant(self._C))

    def __init__(self, parent, C):
        r"""
        Initialize a conjugacy class of directly indecomposable subgroups
        of a symmetric group.

        TESTS::

            sage: from sage.rings.species import ConjugacyClassesOfDirectlyIndecomposableSubgroups
            sage: C = ConjugacyClassesOfDirectlyIndecomposableSubgroups()
            sage: G = C(PermutationGroup([[(1,2),(3,4)],[(1,2),(5,6)]]))
            sage: TestSuite(G).run()
        """
        Element.__init__(self, parent)
        self._C = C

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
        return repr(self._C.gens())


class ConjugacyClassesOfDirectlyIndecomposableSubgroups(UniqueRepresentation, Parent):
    def __init__(self):
        r"""
        Conjugacy classes of directly indecomposable subgroups of a symmetric group.

        TESTS::

            sage: from sage.rings.species import ConjugacyClassesOfDirectlyIndecomposableSubgroups
            sage: C = ConjugacyClassesOfDirectlyIndecomposableSubgroups()
            sage: TestSuite(C).run(max_runs=5) # it takes far too long otherwise
        """
        Parent.__init__(self, category=InfiniteEnumeratedSets())
        self._cache = dict()

    def _element_constructor_(self, x):
        r"""
        Construct a conjugacy class from ``x``.

        INPUT:

        - ``x``, an element of ``self`` or a directly indecomposable
          permutation group.

        EXAMPLES::

            sage: from sage.rings.species import ConjugacyClassesOfDirectlyIndecomposableSubgroups
            sage: C = ConjugacyClassesOfDirectlyIndecomposableSubgroups()
            sage: C(PermutationGroup([("a", "b", "c")]))
            ((1,2,3),)
            sage: C(PermutationGroup([(1, 3, 5)], domain=[1,3,5]))
            ((1,2,3),)
            sage: C(PermutationGroup([[(1,3),(4,7)],[(2,5),(6,8)], [(1,4),(2,5),(3,7)]]))
            ((5,6)(7,8), (1,2)(3,4), (1,3)(2,4)(5,6))

        TESTS:

        Providing a group that decomposes as a direct product raises
        an error::

            sage: C(PermutationGroup([[(1,2,3),(4,5)]]))
            Traceback (most recent call last):
            ...
            ValueError: Permutation Group with generators [(1,2,3)(4,5)] is not directly indecomposable

        Check some trivial cases::

            sage: C(0)
            Traceback (most recent call last):
            ...
            ValueError: unable to convert 0 to Set of conjugacy classes of
             directly indecomposable subgroups

            sage: C(groups.presentation.KleinFour())
            Traceback (most recent call last):
            ...
            ValueError: unable to convert Finitely presented group
             < a, b | a^2, b^2, a^-1*b^-1*a*b > to Set of conjugacy classes of
             directly indecomposable subgroups
        """
        if isinstance(x, PermutationGroup_generic):
            if len(x.disjoint_direct_product_decomposition()) > 1:
                raise ValueError(f"{x} is not directly indecomposable")
            return self.element_class(self, x)
        raise ValueError(f"unable to convert {x} to {self}")

    def _repr_(self):
        r"""
        Return the string representation for ``self``.

        TESTS::

            sage: from sage.rings.species import ConjugacyClassesOfDirectlyIndecomposableSubgroups
            sage: C = ConjugacyClassesOfDirectlyIndecomposableSubgroups(); C
            Set of conjugacy classes of directly indecomposable subgroups
        """
        return "Set of conjugacy classes of directly indecomposable subgroups"

    def __iter__(self):
        r"""
        An iterator over all conjugacy classes of directly indecomposable
        subgroups.

        EXAMPLES::

            sage: from sage.rings.species import ConjugacyClassesOfDirectlyIndecomposableSubgroups
            sage: C = ConjugacyClassesOfDirectlyIndecomposableSubgroups()
            sage: it = iter(C)
            sage: for i in range(5):
            ....:     print(next(it))
            ()
            ((),)
            ((1,2),)
            ((1,2,3),)
            ((2,3), (1,2,3))

        TESTS::

            sage: it = iter(C)
            sage: l = [next(it) for _ in range(100)]  # long time
            sage: len(set(l)) == len(l)  # long time
            True
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

    Element = ConjugacyClassOfDirectlyIndecomposableSubgroups


class AtomicSpeciesElement(Element, WithEqualityById,
                           WithPicklingByInitArgs,
                           metaclass=InheritComparisonClasscallMetaclass):
    r"""
    An atomic species.

    Two atomic species are equal if the underlying groups are
    conjugate, and their domain partitions are equal under the
    conjugating element.
    """
    @staticmethod
    def __classcall__(cls, parent, dis, domain_partition):
        r"""
        Normalize the input for unique representation.

        TESTS::

            sage: from sage.rings.species import AtomicSpecies
            sage: A = AtomicSpecies("X, Y")

        Check that the domain is irrelevant::

            sage: G = PermutationGroup([[("a", "b", "c", "d"), ("e", "f")]])
            sage: a = A(G, {0: "abcd", 1: "ef"}); a  # random
            {((1,2,3,4)(5,6),): ({1, 2, 3, 4}, {5, 6})}
            sage: H = PermutationGroup([[(1,2,3,4), (5,6)]])
            sage: a is A(H, {0: [1,2,3,4], 1: [5,6]})
            True

        The advantage of the unique representation is that we can
        rename the species::

            sage: a.rename("CD(X,Y)"); a
            CD(X,Y)

        We create two different atomic species `a` and `b` with the
        same multicardinality and the same underlying permutation
        group::

            sage: G = PermutationGroup([[(1,2),(3,4),(5,6),(7,8,9,10)]]); G
            Permutation Group with generators [(1,2)(3,4)(5,6)(7,8,9,10)]
            sage: H = PermutationGroup([[(1,2,3,4),(5,6),(7,8),(9,10)]]); H
            Permutation Group with generators [(1,2,3,4)(5,6)(7,8)(9,10)]
            sage: a = A(G, {0: [1,2,3,4], 1: [5,6,7,8,9,10]}); a
            {((1,2,3,4)(5,6)(7,8)(9,10),): ({5, 6, 7, 8}, {1, 2, 3, 4, 9, 10})}
            sage: b = A(H, {0: [1,2,3,4], 1: [5,6,7,8,9,10]}); b
            {((1,2,3,4)(5,6)(7,8)(9,10),): ({1, 2, 3, 4}, {5, 6, 7, 8, 9, 10})}
            sage: c = A(G, {0: [1,2,5,6], 1: [3,4,7,8,9,10]}); c
            {((1,2,3,4)(5,6)(7,8)(9,10),): ({5, 6, 7, 8}, {1, 2, 3, 4, 9, 10})}
            sage: a == b
            False
            sage: a is c
            True
        """
        mc = tuple(len(v) for v in domain_partition)
        domain = list(chain(*map(sorted, domain_partition)))

        def is_equal(elm):
            # mc and dis match because they are part of the key, we
            # construct the mapping between the groups
            elm_domain = list(chain(*map(sorted, elm._dompart)))
            mapping = libgap.MappingPermListList(elm_domain, domain)
            G = PermutationGroup(gap_group=libgap.ConjugateGroup(elm._dis._C,
                                                                 mapping),
                                 domain=elm._dis._C.domain())
            # The conjugated group must be exactly equal to the other group
            return G == dis._C

        key = mc, dis
        if key in parent._cache:
            lookup = parent._cache[key]
            for elm in lookup:
                if is_equal(elm):
                    return elm
        else:
            lookup = parent._cache[key] = []

        elm = WithPicklingByInitArgs.__classcall__(cls, parent, dis, domain_partition)
        lookup.append(elm)
        return elm

    def __init__(self, parent, dis, domain_partition):
        r"""
        Initialize an atomic species.

        INPUT:

        - ``dis`` -- :class:`ConjugacyClassOfDirectlyIndecomposableSubgroups`
        - ``domain_partition`` -- ``dict`` representing the
          assignment of each element of the domain of ``dis`` to
          a "variable"

        TESTS::

            sage: from sage.rings.species import AtomicSpecies
            sage: A = AtomicSpecies("X")
            sage: G = PermutationGroup([[(1,3),(4,7)], [(2,5),(6,8)], [(1,4),(2,5),(3,7)]])
            sage: TestSuite(A(G)).run()

            sage: loads(dumps(A(G))) is A(G)
            True
        """
        Element.__init__(self, parent)
        self._dis = dis
        self._dompart = domain_partition
        self._mc = tuple(len(v) for v in self._dompart)
        self._tc = sum(self._mc)

    def __hash__(self):
        """
        Return a hash of ``self``.

        EXAMPLES::

            sage: from sage.rings.species import AtomicSpecies
            sage: A = AtomicSpecies("X")
            sage: G = PermutationGroup([[(1,2),(3,4)]])
            sage: H = PermutationGroup([[(1,3),(2,4)], [(1,2),(3,4)]])
            sage: K = PermutationGroup([[(1,2,3,4)]])
            sage: hash(A(G)) == hash(A(H))
            False
            sage: hash(A(H)) == hash(A(K))
            True
            sage: A(H) == A(K)
            False
        """
        return hash((self._dis, self._dompart))

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        TESTS::

            sage: from sage.rings.species import AtomicSpecies
            sage: A = AtomicSpecies("X, Y")
            sage: G = PermutationGroup([[(1,2),(3,4),(5,6),(7,8,9,10)]]); G
            Permutation Group with generators [(1,2)(3,4)(5,6)(7,8,9,10)]
            sage: a = A(G, {0: [1,2,3,4], 1: [5,6,7,8,9,10]}); a
            {((1,2,3,4)(5,6)(7,8)(9,10),): ({5, 6, 7, 8}, {1, 2, 3, 4, 9, 10})}
            sage: a = A(G, {1: [1,2,3,4,5,6,7,8,9,10]}); a
            {((1,2,3,4)(5,6)(7,8)(9,10),): ({}, {1, 2, 3, 4, 5, 6, 7, 8, 9, 10})}
        """
        if self.parent()._arity == 1:
            return "{" + f"{self._dis}" + "}"
        dompart = ', '.join("{" + repr(sorted(b))[1:-1] + "}"
                            for b in self._dompart)
        return "{" + f"{self._dis}: ({dompart})" + "}"

    def grade(self):
        r"""
        Return the grade of ``self``.

        EXAMPLES::

            sage: from sage.rings.species import AtomicSpecies
            sage: A = AtomicSpecies("X, Y")
            sage: G = PermutationGroup([[(1,2),(3,4),(5,6),(7,8,9,10)]])
            sage: a = A(G, {0: [1,2,3,4], 1: [5,6,7,8,9,10]})
            sage: a.grade()
            [4, 6]
        """
        S = self.parent().grading_set()
        return S(self._mc)

    def _richcmp_(self, other, op):
        r"""
        Compare ``self`` with ``other`` with respect to the comparison
        operator ``op``.

        ``self`` is less than or equal to ``other`` if it is
        conjugate to a subgroup of ``other`` in the parent group.

        EXAMPLES::

            sage: from sage.rings.species import AtomicSpecies
            sage: A = AtomicSpecies("X")
            sage: A(DihedralGroup(4)) < A(CyclicPermutationGroup(4))
            True

        We create the poset of atomic species of degree four::

            sage: P = Poset([A.subset(4), lambda b, c: b <= c])
            sage: len(P.cover_relations())
            7
            sage: P.cover_relations()  # random
            [[E_4, Eo_4],
             [E_4, P_4],
             [Eo_4, Pb_4],
             [P_4, C_4],
             [P_4, Pb_4],
             [C_4, {((1,2)(3,4),)}],
             [Pb_4, {((1,2)(3,4),)}]]

        TESTS::

            sage: A = AtomicSpecies("X")
            sage: [(a, b) for a, b in Subsets(A.subset(4), 2) if (a < b) != (b > a)]
            []
            sage: [(a, b) for a, b in Subsets(A.subset(4), 2) if (a <= b) != (b >= a)]
            []
            sage: A = AtomicSpecies("X, Y")
            sage: [(a, b) for a, b in Subsets(A.subset(3), 2) if (a < b) != (b > a)]
            []
        """
        if op is op_EQ:
            return self is other
        if op is op_NE:
            return self is not other
        if op is op_LE:
            return self is other or self < other
        if op is op_GE:
            return other <= self
        if op is op_GT:
            return other < self
        if op is op_LT:
            # the arities match because the parents are equal
            if self._mc != other._mc:
                # X should come before Y
                return (sum(self._mc) < sum(other._mc)
                        or (sum(self._mc) == sum(other._mc)
                            and self._mc > other._mc))
            S = SymmetricGroup(sum(self._mc)).young_subgroup(self._mc)
            # conjugate self and other to match S
            g = list(chain.from_iterable(self._dompart))
            conj_self = PermutationGroupElement(g).inverse()
            G = libgap.ConjugateGroup(self._dis._C, conj_self)
            h = list(chain.from_iterable(other._dompart))
            conj_other = PermutationGroupElement(h).inverse()
            H = libgap.ConjugateGroup(other._dis._C, conj_other)
            return GAP_FAIL != libgap.ContainedConjugates(S, G, H, True)


class AtomicSpecies(UniqueRepresentation, Parent):
    """
    The class of atomic species.
    """
    @staticmethod
    def __classcall__(cls, names):
        """
        Normalize the arguments.

        TESTS::

            sage: from sage.rings.species import AtomicSpecies
            sage: A1 = AtomicSpecies("X")
            sage: A2 = AtomicSpecies("Y")
            sage: A3 = AtomicSpecies("X, Y")
            sage: A4 = AtomicSpecies(["X", "Y"])
            sage: A1 == A2
            False
            sage: A3 is A4
            True
        """
        names = normalize_names(-1, names)
        return super().__classcall__(cls, names)

    def __init__(self, names):
        r"""
        Initialize the class of (multivariate) atomic species.

        INPUT:

        - ``names`` -- an iterable of ``k`` strings for the sorts of
          the species

        TESTS:

        We have to exclude `_test_graded_components`, because
        :meth:`~sage.combinat.integer_vector.IntegerVectors.some_elements`
        yields degrees that are too large::

            sage: from sage.rings.species import AtomicSpecies
            sage: A1 = AtomicSpecies(["X"])
            sage: A2 = AtomicSpecies(["X", "Y"])
            sage: TestSuite(A1).run(skip="_test_graded_components")
            sage: TestSuite(A2).run(skip="_test_graded_components")

        """
        category = SetsWithGrading().Infinite()
        Parent.__init__(self, names=names, category=category)
        self._arity = len(names)
        self._cache = dict()
        self._renamed = set()  # the degrees that have been renamed already

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        TESTS::

            sage: from sage.rings.species import AtomicSpecies
            sage: AtomicSpecies("X")
            Atomic species in X
            sage: AtomicSpecies("X, Y")
            Atomic species in X, Y
        """
        if len(self._names) == 1:
            return f"Atomic species in {self._names[0]}"
        return f"Atomic species in {', '.join(self._names)}"

    def _element_constructor_(self, G, pi=None):
        r"""
        Construct the `k`-variate atomic species with the given data.

        INPUT:

        - ``G`` -- element of ``self`` (in this case ``pi`` must be ``None``)
          or permutation group
        - ``pi`` -- `k`-tuple or list of iterables or a dict mapping
          sorts to iterables whose union is the domain; if `k=1`,
          ``pi`` can be omitted

        EXAMPLES::

            sage: from sage.rings.species import AtomicSpecies
            sage: A = AtomicSpecies("X, Y")
            sage: A(DihedralGroup(5), {0: [1,2,3,4,5]})
            P_5(X)

            sage: G = PermutationGroup([[(1,2),(3,4,5,6)]])
            sage: A(G, {0: [1,2], 1: [3,4,5,6]})  # random
            {((1,2,3,4)(5,6),): ({5, 6}, {1, 2, 3, 4})}

            sage: A(G, ([1,2], [3,4,5,6]))  # random
            {((1,2,3,4)(5,6),): ({5, 6}, {1, 2, 3, 4})}

        TESTS::

            sage: G = PermutationGroup([[(1,3),(5,7)]], domain=[1,3,5,7])
            sage: A(G, ([1,3], [5,7]))
            {((1,2)(3,4),): ({1, 2}, {3, 4})}

        Test that errors are raised on some possible misuses::

            sage: A = AtomicSpecies("X, Y")
            sage: G = PermutationGroup([[(1,2), (3,4,5,6)]])
            sage: a = A(G, {0:[1,2], 1:[3,4,5,6]})
            sage: A(a) == a
            True
            sage: A(a, {0:[3,4,5,6], 1:[1,2]})
            Traceback (most recent call last):
            ...
            ValueError: cannot reassign sorts to an atomic species

            sage: A(G)
            Traceback (most recent call last):
            ...
            ValueError: the assignment of sorts to the domain elements must be
             provided

            sage: A(G, {0:[1,2], 1:[]})
            Traceback (most recent call last):
            ...
            ValueError: values of pi (=dict_values([[1, 2], []])) must
             partition the domain of G (={1, 2, 3, 4, 5, 6})

            sage: A(G, {1:[1,2], 2:[3,4,5,6]})
            Traceback (most recent call last):
            ...
            ValueError: keys of pi (=dict_keys([1, 2])) must be in range(2)

            sage: A(G, {0:[1], 1:[2,3,4,5,6]})
            Traceback (most recent call last):
            ...
            ValueError: All elements of orbit (1, 2) must have the same sort

            sage: A(0)
            Traceback (most recent call last):
            ...
            ValueError: 0 must be a permutation group
        """
        if parent(G) == self:
            # pi cannot be None because of framework
            raise ValueError("cannot reassign sorts to an atomic species")
        if not isinstance(G, PermutationGroup_generic):
            raise ValueError(f"{G} must be a permutation group")
        if pi is None:
            if self._arity == 1:
                pi = {0: G.domain()}
            else:
                raise ValueError("the assignment of sorts to the domain elements must be provided")
        elif not isinstance(pi, dict):
            pi = {i: v for i, v in enumerate(pi)}
        if not set(pi.keys()).issubset(range(self._arity)):
            raise ValueError(f"keys of pi (={pi.keys()}) must be in range({self._arity})")
        if sum(len(p) for p in pi.values()) != len(G.domain()) or set(chain.from_iterable(pi.values())) != set(G.domain()):
            raise ValueError(f"values of pi (={pi.values()}) must partition the domain of G (={G.domain()})")
        for orbit in G.orbits():
            if not any(set(orbit).issubset(p) for p in pi.values()):
                raise ValueError(f"All elements of orbit {orbit} must have the same sort")
        # TODO: perhaps move to AtomicSpeciesElement.__classcall__
        dis_elm = ConjugacyClassesOfDirectlyIndecomposableSubgroups()(G)
        mapping = {v: i for i, v in enumerate(G.domain(), 1)}
        sorted_orbits = sorted(G.orbits(), key=len, reverse=True)
        mapping2 = PermutationGroupElement([mapping[e] for o in sorted_orbits for e in o]).inverse()
        dompart = [frozenset(mapping2(mapping[x]) for x in pi.get(k, []))
                   for k in range(self._arity)]
        elm = self.element_class(self, dis_elm, tuple(dompart))
        if elm._tc not in self._renamed:
            self._rename(elm._tc)
        return elm

    def _rename(self, n):
        r"""
        Names for common species.

        EXAMPLES::

            sage: from sage.rings.species import AtomicSpecies
            sage: A = AtomicSpecies(["X", "Y"])
            sage: A(SymmetricGroup(4), {0: range(1, 5)})  # indirect doctest
            E_4(X)
            sage: A(SymmetricGroup(4), {1: range(1, 5)})
            E_4(Y)
            sage: A(CyclicPermutationGroup(4), {0: range(1, 5)})
            C_4(X)
            sage: A(CyclicPermutationGroup(4), {1: range(1, 5)})
            C_4(Y)
            sage: A(DihedralGroup(4), {0: range(1, 5)})
            P_4(X)
            sage: A(DihedralGroup(4), {1: range(1, 5)})
            P_4(Y)
            sage: A(AlternatingGroup(4), {0: range(1, 5)})
            Eo_4(X)
            sage: A(AlternatingGroup(4), {1: range(1, 5)})
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
                self(SymmetricGroup(1), {i: [1]}).rename(self._names[i])

            if self._arity == 1:
                sort = ""
            else:
                sort = f"({self._names[i]})"

            if n >= 2:
                self(SymmetricGroup(n),
                     {i: range(1, n+1)}).rename(f"E_{n}" + sort)

            if n >= 3:
                self(CyclicPermutationGroup(n),
                     {i: range(1, n+1)}).rename(f"C_{n}" + sort)

            if n >= 4:
                self(DihedralGroup(n),
                     {i: range(1, n+1)}).rename(f"P_{n}" + sort)

            if n >= 4:
                self(AlternatingGroup(n),
                     {i: range(1, n+1)}).rename(f"Eo_{n}" + sort)

            if n >= 4 and not n % 2:
                gens = [[(i, n-i+1) for i in range(1, n//2 + 1)],
                        [(i, i+1) for i in range(1, n, 2)]]
                self(PermutationGroup(gens),
                     {i: range(1, n+1)}).rename(f"Pb_{n}" + sort)

    def __contains__(self, x):
        r"""
        Return whether ``x`` is in ``self``.

        TESTS::

            sage: from sage.rings.species import AtomicSpecies
            sage: A = AtomicSpecies("X")
            sage: G = PermutationGroup([[(1,2)], [(3,4)]]); G
            Permutation Group with generators [(3,4), (1,2)]
            sage: G.disjoint_direct_product_decomposition()
            {{1, 2}, {3, 4}}
            sage: G in A
            False

        For convenience, directly indecomposable permutation groups
        are regarded as being in `AtomicSpecies`::

            sage: G = PermutationGroup([(1,2)])
            sage: G in AtomicSpecies("X")
            True
            sage: A(G) in A
            True
            sage: G in AtomicSpecies("X, Y")
            False
            sage: (G, {0: [1,2]}) in AtomicSpecies("X, Y")
            True
            sage: (G, {3: [1,2]}) in AtomicSpecies("X, Y")
            False
            sage: (G, {0: [1]}) in AtomicSpecies("X, Y")
            False
            sage: (G, {0: [1], 1: [2]}) in AtomicSpecies("X, Y")
            False
            sage: (0, {0: []}) in AtomicSpecies("X, Y")
            False
        """
        if parent(x) == self:
            return True
        if isinstance(x, PermutationGroup_generic):
            if self._arity == 1:
                G = x
                pi = {0: G.domain()}
            else:
                return False
        else:
            G, pi = x
        if not isinstance(G, PermutationGroup_generic):
            return False
        if not set(pi.keys()).issubset(range(self._arity)):
            return False
        if (sum(len(p) for p in pi.values()) != len(G.domain())
            or set(chain.from_iterable(pi.values())) != set(G.domain())):
            return False
        for orbit in G.orbits():
            if not any(set(orbit).issubset(p) for p in pi.values()):
                return False
        return len(G.disjoint_direct_product_decomposition()) <= 1

    def grading_set(self):
        r"""
        Return the set of non-negative integer vectors, whose length is
        the arity of ``self``.

        EXAMPLES::

            sage: from sage.rings.species import AtomicSpecies
            sage: AtomicSpecies(["X"]).grading_set()
            Integer vectors of length 1
        """
        return IntegerVectors(length=self._arity)

    def subset(self, size):
        """
        Return the set of atomic species with given total cardinality.

        EXAMPLES::

            sage: from sage.rings.species import AtomicSpecies
            sage: A = AtomicSpecies(["X", "Y"])
            sage: sorted(A.subset(3))
            [E_3(X), C_3(X), E_3(Y), C_3(Y)]
        """
        result = Set()
        for grade in IntegerVectors(size, length=self._arity):
            result = result.union(self.graded_component(grade))
        return result

    def graded_component(self, grade):
        """
        Return the set of atomic species with given multicardinality.

        EXAMPLES::

            sage: from sage.rings.species import AtomicSpecies
            sage: A = AtomicSpecies(["X", "Y"])
            sage: len(A.graded_component([2, 3]))
            1
            sage: A.graded_component([2, 3])  # random
            {{((2,3)(4,5), (1,2,3)): ({4, 5}, {1, 2, 3})}}
        """
        assert len(grade) == self._arity
        n = sum(grade)
        S = SymmetricGroup(n).young_subgroup(grade)
        dom = S.domain()
        dompart = {i: dom[sum(grade[:i]): sum(grade[:i+1])] for i in range(len(grade))}
        return Set([self(G, dompart) for G in S.conjugacy_classes_subgroups()
                    if len(G.disjoint_direct_product_decomposition()) <= 1])

    @cached_method
    def an_element(self):
        """
        Return an element of ``self``.

        TESTS::

            sage: from sage.rings.species import AtomicSpecies
            sage: A = AtomicSpecies("X")
            sage: A.an_element()
            E_2

            sage: A = AtomicSpecies("X, Y")
            sage: a = A.an_element(); a
            {((1,2)(3,4),): ({1, 2}, {3, 4})}

            sage: a.rename("E_2(XY)")
            sage: a
            E_2(XY)
        """
        G = PermutationGroup([[(2 * i + 1, 2 * i + 2) for i in range(self._arity)]])
        m = {i: [2 * i + 1, 2 * i + 2] for i in range(self._arity)}
        return self._element_constructor_(G, m)

    Element = AtomicSpeciesElement


def _stabilizer_subgroups(G, X, a):
    r"""
    Return subgroups conjugate to the stabilizer subgroups of the
    given (left) group action.

    INPUT:

    - ``G`` -- the acting group
    - ``X`` -- the set ``G`` is acting on
    - ``a`` -- the (left) action

    EXAMPLES::

        sage: from sage.rings.species import _stabilizer_subgroups
        sage: S = SymmetricGroup(4)
        sage: X = SetPartitions(S.degree(), [2,2])
        sage: a = lambda g, x: SetPartition([[g(e) for e in b] for b in x])
        sage: _stabilizer_subgroups(S, X, a)
        [Permutation Group with generators [(1,2), (1,3)(2,4)]]

        sage: S = SymmetricGroup(8)
        sage: X = SetPartitions(S.degree(), [3,3,2])
        sage: _stabilizer_subgroups(S, X, a)
        [Permutation Group with generators [(7,8), (6,7), (4,5), (1,3)(2,6)(4,7)(5,8), (1,3)]]

        sage: S = SymmetricGroup(4)
        sage: X = SetPartitions(S.degree(), 2)
        sage: _stabilizer_subgroups(S, X, a)
        [Permutation Group with generators [(1,4), (1,3,4)],
         Permutation Group with generators [(1,3)(2,4), (1,4)]]

    Let us keep 3 and 6 in separate blocks, and check that the
    returned subgroups have the proper domains::

        sage: S = SymmetricGroup([1,2,4,5,3,6]).young_subgroup([4, 2])
        sage: X = [pi for pi in SetPartitions(6, [3,3]) if all(sum(1 for e in b if e % 3) == 2 for b in pi)]
        sage: _stabilizer_subgroups(S, X, a)
        [Permutation Group with generators [(1,2), (4,5), (1,4)(2,5)(3,6)]]

    TESTS::

        sage: _stabilizer_subgroups(SymmetricGroup(2), [1], lambda pi, H: H)
        [Permutation Group with generators [(1,2)]]
    """
    from sage.combinat.cyclic_sieving_phenomenon import orbit_decomposition
    to_gap = {x: i for i, x in enumerate(X, 1)}

    g_orbits = [orbit_decomposition(list(to_gap), lambda x: a(g, x))
                for g in G.gens()]

    gens = [PermutationGroupElement([tuple([to_gap[x] for x in o])
                                     for o in g_orbit])
            for g_orbit in g_orbits]
    result = []
    M = set(range(1, len(to_gap) + 1))
    while M:
        p = M.pop()
        OS = libgap.OrbitStabilizer(G, p, G.gens(), gens)
        result.append(PermutationGroup(gap_group=OS["stabilizer"], domain=G.domain()))
        M.difference_update(OS["orbit"].sage())
    return result


class MolecularSpecies(IndexedFreeAbelianMonoid):
    """
    The monoid of (multivariate) molecular species.
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
            sage: M(G, {0: [5,6], 1: [1,2,3,4]})
            E_2(X)*{((1,2)(3,4),): ({}, {1, 2, 3, 4})}

        TESTS:

        We have to exclude `_test_graded_components`, because
        :meth:`~sage.combinat.integer_vector.IntegerVectors.some_elements`
        yields degrees that are too large::

            sage: M1 = MolecularSpecies("X")
            sage: TestSuite(M1).run(skip="_test_graded_components")
            sage: M2 = MolecularSpecies(["X", "Y"])
            sage: TestSuite(M2).run(skip="_test_graded_components")
        """
        IndexedFreeAbelianMonoid.__init__(self, *args, **kwds)
        self._arity = args[0]._arity

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

    def _element_constructor_(self, G, pi=None):
        r"""
        Construct the `k`-variate molecular species with the given data.

        INPUT:

        - ``G`` -- element of ``self`` (in this case ``pi`` must be ``None``)
          permutation group, or pair ``(X, a)`` consisting of a
          finite set and a transitive action
        - ``pi`` -- ``dict`` mapping sorts to iterables whose union is the
          domain of ``G`` (if ``G`` is a permutation group) or `X` (if ``G``
          is a pair ``(X, a)``); if `k=1`, ``pi`` can be omitted

        EXAMPLES:

        Create a molecular species given a permutation group::

            sage: from sage.rings.species import MolecularSpecies
            sage: M = MolecularSpecies(["X", "Y"])
            sage: G = PermutationGroup([[(1,2)], [(3,4)]])
            sage: M(G, {0: [1,2], 1: [3,4]})
            E_2(X)*E_2(Y)

        Create a molecular species given an action::

            sage: M = MolecularSpecies("X")
            sage: a = lambda g, x: SetPartition([[g(e) for e in b] for b in x])
            sage: X = SetPartitions(4, [2, 2])
            sage: M((X, a), {0: X.base_set()})
            P_4

            sage: X = SetPartitions(8, [4, 2, 2])
            sage: M((X, a), {0: X.base_set()})
            E_4*P_4

        TESTS::

            sage: M = MolecularSpecies(["X", "Y"])
            sage: M(CyclicPermutationGroup(4), {0: [1,2], 1: [3,4]})
            Traceback (most recent call last):
            ...
            ValueError: All elements of orbit (1, 2, 3, 4) must have the same sort

            sage: G = PermutationGroup([[(2,3),(4,5)]], domain=[2,3,4,5])
            sage: M(G, {0: [2, 3], 1: [4, 5]})
            E_2(XY)

            sage: X = SetPartitions(4, 2)
            sage: a = lambda g, x: SetPartition([[g(e) for e in b] for b in x])
            sage: M((X, a), {0: [1,2], 1: [3,4]})
            Traceback (most recent call last):
            ...
            ValueError: Action is not transitive

            sage: G = PermutationGroup([[(1,3),(5,7)]], domain=[1,3,5,7])
            sage: M(G, ([1,3], [5,7]))
            E_2(XY)

            sage: G = PermutationGroup([[(1,2), (3,4,5,6)]])
            sage: a = M(G, {0:[1,2], 1:[3,4,5,6]})
            sage: M(a) == a
            True
            sage: M(a, {0:[3,4,5,6], 1:[1,2]})
            Traceback (most recent call last):
            ...
            ValueError: cannot reassign sorts to a molecular species

            sage: M(G)
            Traceback (most recent call last):
            ...
            ValueError: the assignment of sorts to the domain elements must be
             provided

            sage: M(G, {0:[1,2], 1:[]})
            Traceback (most recent call last):
            ...
            ValueError: values of pi (=dict_values([[1, 2], []])) must
             partition the domain of G (={1, 2, 3, 4, 5, 6})

            sage: M(G, {1:[1,2], 2:[3,4,5,6]})
            Traceback (most recent call last):
            ...
            ValueError: keys of pi (=dict_keys([1, 2])) must be in range(2)

            sage: M(0)
            Traceback (most recent call last):
            ...
            ValueError: 0 must be a permutation group or a pair specifying a
             group action on the given domain pi=None
        """
        if parent(G) == self:
            # pi cannot be None because of framework
            raise ValueError("cannot reassign sorts to a molecular species")

        if isinstance(G, PermutationGroup_generic):
            if pi is None:
                if self._arity == 1:
                    pi = {0: G.domain()}
                else:
                    raise ValueError("the assignment of sorts to the domain elements must be provided")
            elif not isinstance(pi, dict):
                pi = {i: v for i, v in enumerate(pi)}
            domain = [e for p in pi.values() for e in p]
            if len(domain) != len(set(domain)) or set(G.domain()) != set(domain):
                raise ValueError(f"values of pi (={pi.values()}) must partition the domain of G (={G.domain()})")
            components = G.disjoint_direct_product_decomposition()
            elm = self.one()
            for component in components:
                gens = [[cyc for cyc in gen.cycle_tuples() if cyc[0] in component]
                        for gen in G.gens()]
                H = PermutationGroup([gen for gen in gens if gen],
                                     domain=component)
                dompart = {k: [e for e in v if e in component]
                           for k, v in pi.items()}
                elm *= self.gen(self._indices(H, dompart))
            return elm

        try:
            X, a = G
        except TypeError:
            raise ValueError(f"{G} must be a permutation group or a pair specifying a group action on the given domain pi={pi}")
        L = [len(pi.get(i, [])) for i in range(self._arity)]
        S = SymmetricGroup(sum(L)).young_subgroup(L)
        H = _stabilizer_subgroups(S, X, a)
        if len(H) > 1:
            raise ValueError("Action is not transitive")
        return self(H[0], pi)

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

    def subset(self, size):
        """
        Return the set of molecular species with given total cardinality.

        EXAMPLES::

            sage: from sage.rings.species import MolecularSpecies
            sage: M = MolecularSpecies(["X", "Y"])
            sage: M.subset(3)  # random
            {X*E_2(Y), X*Y^2, C_3(Y), E_3(X), Y^3, Y*E_2(Y), C_3(X), X^2*Y,
             E_3(Y), E_2(X)*Y, X*E_2(X), X^3}
        """
        result = Set()
        for grade in IntegerVectors(size, length=self._arity):
            result = result.union(self.graded_component(grade))
        return result

    def graded_component(self, grade):
        """
        Return the set of molecular species with given multicardinality.

        EXAMPLES::

            sage: from sage.rings.species import MolecularSpecies
            sage: M = MolecularSpecies(["X", "Y"])
            sage: M.graded_component([3, 2])  # random
            {E_3(X)*Y^2, X^3*Y^2, X*E_2(X)*E_2(Y), X^3*E_2(Y),
             {((1,2,3), (1,3)(4,5)): ({1, 2, 3}, {4, 5})},
             X*{((1,2)(3,4),): ({1, 2}, {3, 4})}, X*E_2(X)*Y^2, E_3(X)*E_2(Y),
             C_3(X)*Y^2, C_3(X)*E_2(Y)}
        """
        assert len(grade) == self._arity
        n = sum(grade)
        S = SymmetricGroup(n).young_subgroup(grade)
        dom = S.domain()
        dompart = {i: dom[sum(grade[:i]): sum(grade[:i+1])] for i in range(len(grade))}
        return Set([self(G, dompart) for G in S.conjugacy_classes_subgroups()])

    class Element(IndexedFreeAbelianMonoidElement):
        def __init__(self, parent, x):
            r"""
            Initialize a molecular species.

            INPUT:

            - ``x`` -- a dictionary mapping atomic species to exponents

            EXAMPLES::

                sage: from sage.rings.species import MolecularSpecies
                sage: M = MolecularSpecies("X")
                sage: M(CyclicPermutationGroup(3))  # indirect doctest
                C_3

            TESTS::

                sage: X = M(CyclicPermutationGroup(3))
                sage: C3 = M(CyclicPermutationGroup(3))
                sage: TestSuite(X*C3).run()
            """
            super().__init__(parent, x)

        @cached_method
        def grade(self):
            r"""
            Return the grade of ``self``.

            EXAMPLES::

                sage: from sage.rings.species import MolecularSpecies
                sage: M = MolecularSpecies("X, Y")
                sage: G = PermutationGroup([[(1,2),(3,4),(5,6),(7,8,9,10)]]); G
                Permutation Group with generators [(1,2)(3,4)(5,6)(7,8,9,10)]
                sage: a = M(G, {0: [1,2,3,4], 1: [5,6,7,8,9,10]}); a
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

        def _richcmp_(self, other, op):
            r"""
            Compare ``self`` with ``other`` with respect to the
            comparison operator ``op``.

            ``self`` is less than or equal to ``other`` if it is
            conjugate to a subgroup of ``other`` in the parent group.

            EXAMPLES::

                sage: from sage.rings.species import MolecularSpecies
                sage: M = MolecularSpecies("X")
                sage: M(DihedralGroup(4)) < M(CyclicPermutationGroup(4))
                True

            We create the lattice of molecular species of degree four::

                sage: P = Poset([M.subset(4), lambda b, c: b <= c])
                sage: len(P.cover_relations())
                17
                sage: P.cover_relations()  # random
                [[E_4, P_4],
                 [E_4, Eo_4],
                 [E_4, X*E_3],
                 [P_4, C_4],
                 [P_4, E_2^2],
                 [P_4, Pb_4],
                 [C_4, {((1,2)(3,4),)}],
                 [E_2^2, {((1,2)(3,4),)}],
                 [E_2^2, X^2*E_2],
                 [Eo_4, Pb_4],
                 [Eo_4, X*C_3],
                 [Pb_4, {((1,2)(3,4),)}],
                 [{((1,2)(3,4),)}, X^4],
                 [X*E_3, X^2*E_2],
                 [X*E_3, X*C_3],
                 [X^2*E_2, X^4],
                 [X*C_3, X^4]]

            TESTS::

                sage: [(a, b) for a, b in Subsets(M.subset(4), 2) if (a < b) != (b > a)]
                []
                sage: [(a, b) for a, b in Subsets(M.subset(4), 2) if (a <= b) != (b >= a)]
                []

                sage: M = MolecularSpecies("S, T")
                sage: S = M(SymmetricGroup(1), {0: [1]})
                sage: T = M(SymmetricGroup(1), {1: [1]})
                sage: E2S = M(SymmetricGroup(2), {0: [1,2]})
                sage: T * E2S < S^2 * T
                True
            """
            if op is op_EQ or op is op_NE:
                return super()._richcmp_(other, op)
            if op is op_LE:
                return self == other or self < other
            if op is op_GE:
                return other <= self
            if op is op_GT:
                return other < self
            if op is op_LT:
                # the arities match because the parents are equal
                if self.grade() != other.grade():
                    # X should come before Y
                    return (sum(self.grade()) < sum(other.grade())
                            or (sum(self.grade()) == sum(other.grade())
                                and self.grade() > other.grade()))

                S = SymmetricGroup(sum(self.grade())).young_subgroup(self.grade())
                # conjugate self and other to match S
                G, G_dompart = self.group_and_partition()
                g = list(chain.from_iterable(G_dompart))
                conj_self = PermutationGroupElement(g).inverse()
                G = libgap.ConjugateGroup(G, conj_self)
                H, H_dompart = other.group_and_partition()
                h = list(chain.from_iterable(H_dompart))
                conj_other = PermutationGroupElement(h).inverse()
                H = libgap.ConjugateGroup(H, conj_other)
                return GAP_FAIL != libgap.ContainedConjugates(S, G, H, True)

        @cached_method
        def group_and_partition(self):
            r"""
            Return the (transitive) permutation group
            corresponding to ``self``, together with the partition of
            the domain into sorts.

            EXAMPLES::

                sage: from sage.rings.species import MolecularSpecies
                sage: M = MolecularSpecies("X,Y")
                sage: G = PermutationGroup([[(1,2),(3,4)], [(5,6)]])
                sage: A = M(G, {0: [5,6], 1: [1,2,3,4]}); A
                E_2(X)*{((1,2)(3,4),): ({}, {1, 2, 3, 4})}
                sage: A.group_and_partition()
                (Permutation Group with generators [(3,4)(5,6), (1,2)],
                 (frozenset({1, 2}), frozenset({3, 4, 5, 6})))

            Note that we cannot rely on blocks of the partition being
            consecutive::

                sage: A = M(PermutationGroup([(3,4)]), {0:[1,3,4], 1:[2]})
                sage: A
                X*Y*E_2(X)
                sage: A.group_and_partition()[1]
                (frozenset({1, 3, 4}), frozenset({2}))

            TESTS::

                sage: B = M(PermutationGroup([(1,2,3)]), {0: [1,2,3]}); B
                C_3(X)
                sage: B.group_and_partition()
                (Permutation Group with generators [(1,2,3)],
                 (frozenset({1, 2, 3}), frozenset()))

                sage: G = PermutationGroup([[(1,2),(3,4)], [(5,6)]])
                sage: A = M(G, {0: [5,6], 1: [1,2,3,4]})
                sage: A * B
                E_2(X)*C_3(X)*{((1,2)(3,4),): ({}, {1, 2, 3, 4})}
                sage: (A*B).group_and_partition()
                (Permutation Group with generators [(6,7)(8,9), (3,4,5), (1,2)],
                 (frozenset({1, 2, 3, 4, 5}), frozenset({6, 7, 8, 9})))

                sage: C = M(PermutationGroup([(2,3)]), {0: [1], 1: [2,3]}); C
                X*E_2(Y)
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

                sage: F = M(PermutationGroup([(1,2),(3,)]))
                sage: F.group_and_partition()[0].domain()
                {1, 2, 3}

                sage: F = M(AlternatingGroup(2))
                sage: F.group_and_partition()[0].domain()
                {1, 2}
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
                    # split off a single monomial
                    a = list(A._monomial)[0]  # as atomic species
                    b, b_dompart = (A ** (n-1)).group_and_partition()
                    gens = a._dis._C.gens() + shift_gens(b.gens(), a._tc)
                    new_dompart = tuple([frozenset(list(p_a) + [a._tc + e for e in p_b])
                                         for p_a, p_b in zip(a._dompart, b_dompart)])
                    domain = range(1, n * a._tc + 1)
                else:
                    f, f_dompart = (A ** (n // 2)).group_and_partition()
                    tc = sum(len(p) for p in f_dompart)
                    gens = f.gens() + shift_gens(f.gens(), tc)
                    new_dompart = tuple([frozenset(list(p) + [tc + e for e in p])
                                         for p in f_dompart])
                    domain = range(1, 2 * tc + 1)

                G = PermutationGroup(gens, domain=domain)
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
            G = PermutationGroup(gens, domain=range(1, tc_list[-1]+1))
            new_dompart = tuple([frozenset(chain(*[[tc + e for e in p]
                                                   for p, tc in zip(f_dompart, tc_list)]))
                                 for f_dompart in zip(*dompart_list)])

            return G, new_dompart

        def cycle_index(self, parent=None):
            r"""
            Return the cycle index of ``self``.

            This is essentially a variant of
            :meth:`~sage.categories.finite_permutation_groups.FinitePermutationGroups.ParentMethods.cycle_index`
            for subgroups of a Young subgroup of the symmetric group.

            EXAMPLES::

                sage: from sage.rings.species import MolecularSpecies
                sage: M = MolecularSpecies("X,Y")
                sage: G = PermutationGroup([[(1,2),(3,4)], [(5,6)]])
                sage: A = M(G, {0: [5,6], 1: [1,2,3,4]})
                sage: A.cycle_index()
                1/4*p[1, 1] # p[1, 1, 1, 1] + 1/4*p[1, 1] # p[2, 2] + 1/4*p[2] # p[1, 1, 1, 1] + 1/4*p[2] # p[2, 2]

            TESTS:

            Check that we support different parents::

                sage: F = CombinatorialFreeModule(QQ, Partitions())
                sage: P = A.cycle_index(parent=tensor([F, F]))
                sage: P
                1/4*B[[1, 1]] # B[[1, 1, 1, 1]] + 1/4*B[[1, 1]] # B[[2, 2]] + 1/4*B[[2]] # B[[1, 1, 1, 1]] + 1/4*B[[2]] # B[[2, 2]]
                sage: P.parent() is tensor([F, F])
                True

            This parent should be a module with basis indexed by partitions::

                sage: A.cycle_index(parent=QQ)
                Traceback (most recent call last):
                  ...
                ValueError: `parent` should be a module with basis indexed by partitions
            """
            k = self.parent()._arity
            if parent is None:
                p = SymmetricFunctions(QQ).powersum()
                parent = tensor([p]*k)
            elif parent not in Modules.WithBasis:
                raise ValueError("`parent` should be a module with basis indexed by partitions")
            base_ring = parent.base_ring()
            G, dompart = self.group_and_partition()
            dompart_dict = {}
            for i, s in enumerate(dompart):
                dompart_dict.update({e: i for e in s})

            def cycle_type(pi):
                tuples = pi.cycle_tuples(singletons=True)
                cycle_type = [[] for _ in range(k)]
                for c in tuples:
                    cycle_type[dompart_dict[c[0]]].append(len(c))
                return tuple(_Partitions(sorted(c, reverse=True))
                             for c in cycle_type)

            return (parent.sum_of_terms([cycle_type(C.an_element()),
                                         base_ring(C.cardinality())]
                                        for C in G.conjugacy_classes())
                    / G.cardinality())

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
                sage: X = M(SymmetricGroup(1), {0: [1]})
                sage: Y = M(SymmetricGroup(1), {1: [1]})
                sage: (X*Y)(X, Y^2)
                X*Y^2

            A multivariate example::

                sage: M1 = MolecularSpecies("X")
                sage: M2 = MolecularSpecies("X, Y")
                sage: C3 = M1(CyclicPermutationGroup(3))
                sage: X = M2(SymmetricGroup(1), {0: [1]})
                sage: Y = M2(SymmetricGroup(1), {1: [1]})
                sage: C3(X*Y)
                {((1,2,3)(4,5,6),): ({1, 2, 3}, {4, 5, 6})}

            TESTS::

                sage: M = MolecularSpecies("X")
                sage: M.one()()
                Traceback (most recent call last):
                ...
                ValueError: number of args must match arity of self

                sage: M.one()(2)
                Traceback (most recent call last):
                ...
                ValueError: all args must be molecular species

                sage: M2 = MolecularSpecies("X, Y")
                sage: X2 = M2(SymmetricGroup(1), {0: [1]})
                sage: Y2 = M2(SymmetricGroup(1), {1: [1]})
                sage: X = M(SymmetricGroup(1))
                sage: (X2*Y2)(X2, X)
                Traceback (most recent call last):
                ...
                ValueError: all args must have the same parent
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
            dompart = {i: [] for i in range(P._arity)}
            for start, M in zip(starts, Mlist):
                K, K_dompart = M.group_and_partition()
                for i, v in enumerate(K_dompart):
                    dompart[i].extend([start + k for k in v])
                for gen in K.gens():
                    gens.append([tuple(start + k for k in cyc)
                                 for cyc in gen.cycle_tuples()])

            return P(PermutationGroup(gens, domain=range(1, starts[-1] + 1)),
                     dompart)


class PolynomialSpeciesElement(CombinatorialFreeModule.Element):
    r"""
    A (virtual) polynomial species.

    TESTS::

        sage: from sage.rings.species import PolynomialSpecies
        sage: P = PolynomialSpecies(ZZ, ["X"])
        sage: C3 = P(CyclicPermutationGroup(3))
        sage: X = P(SymmetricGroup(1))
        sage: E2 = P(SymmetricGroup(2))
        sage: (E2*X + C3).homogeneous_degree()
        3

        sage: TestSuite(E2*X + C3).run()
    """
    def is_constant(self):
        """
        Return ``True`` if this is a constant polynomial species.

        EXAMPLES::

            sage: from sage.rings.species import PolynomialSpecies
            sage: P = PolynomialSpecies(ZZ, ["X", "Y"])
            sage: X = P(SymmetricGroup(1), {0: [1]})
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

    def is_virtual(self):
        r"""
        Return if ``self`` is a virtual species.

        TESTS::

            sage: from sage.rings.species import PolynomialSpecies
            sage: P = PolynomialSpecies(ZZ, ["X", "Y"])
            sage: X = P(SymmetricGroup(1), {0: [1]})
            sage: Y = P(SymmetricGroup(1), {1: [1]})
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

            sage: from sage.rings.species import PolynomialSpecies
            sage: P = PolynomialSpecies(ZZ, ["X", "Y"])
            sage: X = P(SymmetricGroup(1), {0: [1]})
            sage: Y = P(SymmetricGroup(1), {1: [1]})
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

            sage: from sage.rings.species import PolynomialSpecies
            sage: P = PolynomialSpecies(ZZ, ["X", "Y"])
            sage: X = P(SymmetricGroup(1), {0: [1]})
            sage: Y = P(SymmetricGroup(1), {1: [1]})
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

            sage: from sage.rings.species import PolynomialSpecies
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

            sage: C3.hadamard_product(E2^2)
            0

        TESTS::

            sage: C3.hadamard_product(-C3)
            -2*C_3

            sage: Q = PolynomialSpecies(ZZ, ["Y"])
            sage: P.one().hadamard_product(Q.one())
            Traceback (most recent call last):
            ...
            ValueError: the factors of a Hadamard product must have the same parent
        """
        P = self.parent()
        if P is not other.parent():
            raise ValueError("the factors of a Hadamard product must have the same parent")

        result = P.zero()
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
            dompart = {i: range(x - mc[i] + 1, x + 1)
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
                             dompart)
                result += c * d * new

        return result

    def _compose_with_singletons(self, names, args):
        r"""
        Return the inner sum of Exercise 2.6.16 in [BLL1998]_,
        generalized to the case of arbitrary arity.

        The `k`-sort species ``self`` should be homogeneous.

        INPUT:

            - ``names``, the (flat) list of names of the result

            - ``args``, the sequence of `k` compositions, each of
              which sums to the corresponding degree of ``self``,
              where `k` is the arity of ``self``.

        OUTPUT:

        - the polynomial species `self(X_1 + \dots + X_{m_1}, Y_1
          + \dots + Y_{m_2}, \dots)`, where `m_i` is the number
          of parts of the `i`-th composition, restricted to the
          degrees given by ``args``.

        EXAMPLES::

            sage: from sage.rings.species import PolynomialSpecies
            sage: P = PolynomialSpecies(ZZ, "X")
            sage: C4 = P(CyclicPermutationGroup(4))
            sage: C4._compose_with_singletons("X, Y", [[2, 2]])
            E_2(XY) + X^2*Y^2

            sage: P = PolynomialSpecies(ZZ, ["X", "Y"])
            sage: F = P(PermutationGroup([[(1,2,3), (4,5,6)]]), {0: [1,2,3], 1: [4,5,6]})
            sage: F
            {((1,2,3)(4,5,6),): ({1, 2, 3}, {4, 5, 6})}
            sage: F._compose_with_singletons("X1, X2, X3, Y1, Y2", [[1, 1, 1], [2, 1]])
            6*X1*X2*X3*Y1^2*Y2

        TESTS::

            sage: P = PolynomialSpecies(ZZ, "X")
            sage: O = P.one()
            sage: O._compose_with_singletons("X", [[]])
            1

            sage: F = P(SymmetricGroup(1)) * P(SymmetricGroup(2))
            sage: F._compose_with_singletons(["S", "T"], [[2, 1]])
            T*E_2(S) + S^2*T
            sage: F._compose_with_singletons(["S", "T"], [[1, 2]])
            S*E_2(T) + S*T^2
        """
        # TODO: possibly check that all args are compositions,
        # and that sums match cardinalities
        comp = list(chain.from_iterable(args))
        result_dompart = {i: range(x - comp[i] + 1, x + 1)
                          for i, x in enumerate(accumulate(comp))}
        S_down = SymmetricGroup(sum(comp)).young_subgroup(comp)

        Pn = PolynomialSpecies(self.parent().base_ring(), names)
        result = Pn.zero()
        for M, c in self:
            # Create group of the composition
            # conjugate self.group() so that [1..k] is sort 1, [k+1,..] is sort 2, so on
            G, dompart = M.group_and_partition()
            conj = PermutationGroupElement(list(chain.from_iterable(dompart))).inverse()
            G = libgap.ConjugateGroup(G, conj)

            mc = M.grade()
            tc = sum(mc)
            S_up = SymmetricGroup(tc).young_subgroup(mc)
            taus = libgap.DoubleCosetRepsAndSizes(S_up, S_down, G)

            # sum over double coset representatives.
            summand = Pn.zero()
            for tau, _ in taus:
                H = libgap.Intersection(libgap.ConjugateGroup(G, tau.Inverse()),
                                        S_down)
                K = PermutationGroup(gap_group=H, domain=range(1, tc + 1))
                summand += Pn(K, result_dompart)
            result += c*summand
        return result

    def _compose_with_weighted_singletons(self, names, multiplicities, degrees):
        r"""
        Compute the composition with
        `(\sum_j m_{1,j} X_{1,j}, \sum_j m_{2,j} X_{2,j}, \ldots)`
        in the specified degrees.

        The `k`-sort species ``self`` should be homogeneous.

        INPUT:

        - ``names`` -- the (flat) list of names of the result
        - ``multiplicities`` -- a (flat) list of constants
        - ``degrees`` -- a `k`-tuple of compositions `c_1,
          \ldots, c_k`, such that the size of `c_i` is the
          degree of ``self`` in sort `i`

        EXAMPLES:

        Equation (2.5.41) in [BLL1998]_::

            sage: from sage.rings.species import PolynomialSpecies
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

        Equation (4.60) in [ALL2002]_::

            sage: C4._compose_with_weighted_singletons(["X", "Y"], [1, -1], [[2, 2]])
            -E_2(XY) + 2*X^2*Y^2

        TESTS::

            sage: (C4+E2^2)._compose_with_weighted_singletons(["X"], [-1], [[4]])
            -C_4 + E_2^2 + {((1,2)(3,4),)} - 2*X^2*E_2 + X^4
        """
        P = self.parent()
        if not self.support():
            return P.zero()
        if not self.is_homogeneous():
            raise ValueError("element is not homogeneous")

        left = self._compose_with_singletons(names, degrees)
        Pn = left.parent()
        right = Pn.exponential(multiplicities,
                               list(chain.from_iterable(degrees)))
        return left.hadamard_product(right)

    def __call__(self, *args):
        r"""
        Substitute the arguments into ``self``.

        EXAMPLES::

            sage: from sage.rings.species import PolynomialSpecies
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
            sage: X = P2(SymmetricGroup(1), {0: [1]})
            sage: Y = P2(SymmetricGroup(1), {1: [1]})
            sage: E2(X + Y)
            E_2(X) + X*Y + E_2(Y)

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

    def factor(self):
        r"""
        Return the factorization of this species.

        EXAMPLES::

            sage: from sage.rings.species import PolynomialSpecies
            sage: P = PolynomialSpecies(ZZ, ["X"])
            sage: C3 = P(CyclicPermutationGroup(3))
            sage: X = P(SymmetricGroup(1))
            sage: E2 = P(SymmetricGroup(2))
            sage: f = (3*E2*X + C3)*(2*E2 + C3)
            sage: factor(f)
            (2*E_2 + C_3) * (3*X*E_2 + C_3)
        """
        # find the set of atoms and fix an order
        atoms = list(set(a for m in self.monomial_coefficients()
                         for a in m.support()))
        R = PolynomialRing(self.base_ring(), "x", len(atoms))
        var_dict = dict(zip(atoms, R.gens()))
        # create the polynomial
        poly = R.sum(c * R.prod(var_dict[a] ** e for a, e in m.dict().items())
                     for m, c in self)
        factors = poly.factor()
        unit = self.base_ring()(factors.unit())
        if factors.universe() == self.base_ring():
            return Factorization(factors, unit=unit)
        P = self.parent()
        M = P._indices

        def _from_etuple(e):
            return M.element_class(M, {a: i for a, i in zip(atoms, e) if i})

        factors = [(P.sum_of_terms((_from_etuple(mon), c)
                                   for mon, c in factor.monomial_coefficients().items()),
                    exponent)
                   for factor, exponent in factors]
        return Factorization(factors, unit=unit, sort=False)


class PolynomialSpecies(CombinatorialFreeModule):
    def __classcall__(cls, base_ring, names):
        r"""
        Normalize the arguments.

        TESTS::

            sage: from sage.rings.species import PolynomialSpecies
            sage: P1 = PolynomialSpecies(ZZ, "X, Y")
            sage: P2 = PolynomialSpecies(ZZ, "X, Y")
            sage: P3 = PolynomialSpecies(ZZ, ["X", "Z"])
            sage: P1 is P2
            True
            sage: P1 == P3
            False
        """
        from sage.structure.category_object import normalize_names
        names = normalize_names(-1, names)
        return super().__classcall__(cls, base_ring, names)

    def __init__(self, base_ring, names):
        r"""
        Ring of `k`-variate polynomial (virtual) species.

        TESTS::

            sage: from sage.rings.species import PolynomialSpecies
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
                                         prefix='', bracket=False)
        self._arity = len(names)

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: from sage.rings.species import PolynomialSpecies
            sage: PolynomialSpecies(ZZ, "X")
            Polynomial species in X over Integer Ring
            sage: PolynomialSpecies(ZZ, "X, Y")
            Polynomial species in X, Y over Integer Ring
        """
        names = self._indices._indices._names
        if len(names) == 1:
            return f"Polynomial species in {names[0]} over {self.base_ring()}"
        return f"Polynomial species in {', '.join(names)} over {self.base_ring()}"

    def _element_constructor_(self, G, pi=None):
        r"""
        Construct the `k`-variate polynomial species with the given data.

        INPUT:

        - ``G`` -- element of ``self`` (in this case ``pi`` must be ``None``)
          permutation group, or pair ``(X, a)`` consisting of a
          finite set and an action
        - ``pi`` -- ``dict`` mapping sorts to iterables whose union is the
          domain of ``G`` (if ``G`` is a permutation group) or `X` (if ``G``
          is a pair ``(X, a)``); if `k=1`, ``pi`` can be omitted

        EXAMPLES::

            sage: from sage.rings.species import PolynomialSpecies
            sage: P = PolynomialSpecies(ZZ, ["X", "Y"])
            sage: P(SymmetricGroup(4).young_subgroup([2, 2]), ([1,2], [3,4]))
            E_2(X)*E_2(Y)

            sage: X = SetPartitions(4, 2)
            sage: a = lambda g, x: SetPartition([[g(e) for e in b] for b in x])
            sage: P((X, a), {0: [1,2], 1: [3,4]})
            E_2(X)*E_2(Y) + X^2*E_2(Y) + E_2(XY) + Y^2*E_2(X)

            sage: P = PolynomialSpecies(ZZ, ["X"])
            sage: X = SetPartitions(4, 2)
            sage: a = lambda g, x: SetPartition([[g(e) for e in b] for b in x])
            sage: P((X, a), {0: [1,2,3,4]})
            X*E_3 + P_4

        The species of permutation groups::

            sage: P = PolynomialSpecies(ZZ, ["X"])
            sage: n = 4
            sage: S = SymmetricGroup(n)
            sage: X = S.subgroups()
            sage: def act(pi, G):  # WARNING: returning H does not work because of equality problems
            ....:     H = S.subgroup(G.conjugate(pi).gens())
            ....:     return next(K for K in X if K == H)
            sage: P((X, act), {0: range(1, n+1)})
            4*E_4 + 4*P_4 + E_2^2 + 2*X*E_3
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
            elif not isinstance(pi, dict):
                pi = {i: v for i, v in enumerate(pi)}
            return self._from_dict({self._indices(G, pi): ZZ.one()})

        X, a = G
        L = [len(pi.get(i, [])) for i in range(self._arity)]
        S = SymmetricGroup(sum(L)).young_subgroup(L)
        Hs = _stabilizer_subgroups(S, X, a)
        return self.sum_of_terms((self._indices(H, pi), ZZ.one()) for H in Hs)

    def _first_ngens(self, n):
        r"""
        Used by the preparser for ``F.<x> = ...``.

        We do not use the generic implementation of
        :class:`sage.combinat.CombinatorialFreeModule`, because we do
        not want to implement `gens`.

        EXAMPLES::

            sage: from sage.rings.species import PolynomialSpecies
            sage: P.<X, Y> = PolynomialSpecies(QQ)  # indirect doctest
            sage: X + 2*Y
            X + 2*Y
        """
        B = self.basis()
        return tuple(B[i] for grade in IntegerVectors(1, length=self._arity)
                     for i in self._indices.graded_component(grade))

    def change_ring(self, R):
        r"""
        Return the base change of ``self`` to `R`.

        EXAMPLES::

            sage: from sage.rings.species import PolynomialSpecies
            sage: P = PolynomialSpecies(ZZ, ["X", "Y"])
            sage: P.change_ring(QQ)
            Polynomial species in X, Y over Rational Field
        """
        if R is self.base_ring():
            return self
        return PolynomialSpecies(R, self._indices._indices._names)

    def degree_on_basis(self, m):
        r"""
        Return the degree of the molecular species indexed by ``m``.

        EXAMPLES::

            sage: from sage.rings.species import PolynomialSpecies
            sage: P = PolynomialSpecies(ZZ, ["X", "Y"])
            sage: E4X = P(SymmetricGroup(4), {0: range(1, 5)}); E4X
            E_4(X)
            sage: E4Y = P(SymmetricGroup(4), {1: range(1, 5)}); E4Y
            E_4(Y)
            sage: P.degree_on_basis(E4X.support()[0])
            4
            sage: P.degree_on_basis(E4Y.support()[0])
            4
        """
        return sum(m.grade())

    @cached_method
    def one_basis(self):
        r"""
        Return ``SymmetricGroup(0)``, which indexes the one of this algebra,
        as per :meth:`AlgebrasWithBasis.ParentMethods.one_basis`.

        EXAMPLES::

            sage: from sage.rings.species import PolynomialSpecies
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

        EXAMPLES::

            sage: from sage.rings.species import PolynomialSpecies
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

            sage: from sage.rings.species import PolynomialSpecies
            sage: P = PolynomialSpecies(ZZ, "X")
            sage: L1 = [P(H) for H in SymmetricGroup(3).conjugacy_classes_subgroups()]
            sage: L2 = [P(H) for H in SymmetricGroup(2).conjugacy_classes_subgroups()]
            sage: matrix([[F * G for F in L1] for G in L2])  # indirect doctest
            [    X^5 X^3*E_2 X^2*C_3 X^2*E_3]
            [X^3*E_2 X*E_2^2 E_2*C_3 E_2*E_3]

        TESTS::

            sage: P = PolynomialSpecies(ZZ, "X")
            sage: X = P(SymmetricGroup(1))
            sage: type(list(X^2)[0][1])
            <class 'sage.rings.integer.Integer'>
        """
        return self.element_class(self, {H * K: ZZ.one()})

    @cached_method
    def powersum(self, s, n):
        r"""
        Return the combinatorial powersum species `P_n(X_s)`.

        The species `P_n(X)` is introduced in [Labelle2008]_ as the
        coefficient of `t^n/n` in `\log E(tX)`, where `E` is the
        species of sets.

        EXAMPLES::

            sage: from sage.rings.species import PolynomialSpecies
            sage: P = PolynomialSpecies(ZZ, "X")
            sage: P.powersum(0, 4)
            4*E_4 - 4*X*E_3 - 2*E_2^2 + 4*X^2*E_2 - X^4
        """
        if not (n in ZZ and n > 0):
            raise ValueError("n must be a positive integer")
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

            sage: from sage.rings.species import PolynomialSpecies
            sage: P = PolynomialSpecies(QQ, ["X"])
            sage: P.exponential([3/2], [7])
            3/2*E_7 + 3/4*X*E_6 + 3/4*E_2*E_5 - 3/16*X^2*E_5 + 3/4*E_3*E_4
             - 3/8*X*E_2*E_4 + 3/32*X^3*E_4 - 3/16*X*E_3^2 - 3/16*E_2^2*E_3
             + 9/32*X^2*E_2*E_3 - 15/256*X^4*E_3 + 3/32*X*E_2^3
             - 15/128*X^3*E_2^2 + 21/512*X^5*E_2 - 9/2048*X^7

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

            sage: P = PolynomialSpecies(QQ, ["X"])
            sage: P.exponential([1], [0]).parent()
            Polynomial species in X over Rational Field
        """
        def stretch(c, k):
            r"""
            Substitute in ``c`` all variables appearing in the
            base ring with their ``k``-th power.
            """
            if callable(c):
                B = self.base_ring()
                return c(*[g ** k for g in B.gens() if g != B.one()])
            return c

        def factor(s, c, d):
            r"""
            Return `E(c X_s)_d`.

            We use Proposition 2 in [Labelle2008]_.
            """
            return self.sum(~ mu.centralizer_size()
                            * self.prod(stretch(c, k)
                                        * self.powersum(s, k) for k in mu)
                            for mu in Partitions(d))

        return self.prod(factor(s, multiplicities[s], degrees[s])
                         for s in range(self._arity))

    Element = PolynomialSpeciesElement