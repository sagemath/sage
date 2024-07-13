from sage.misc.cachefunc import cached_method
from sage.structure.parent import Parent
from sage.structure.element import Element
from sage.structure.element import parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.rings.integer_ring import ZZ
from sage.categories.sets_cat import cartesian_product
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.categories.sets_with_grading import SetsWithGrading
from sage.categories.graded_algebras_with_basis import GradedAlgebrasWithBasis
from sage.categories.algebras import Algebras
from sage.groups.perm_gps.permgroup import PermutationGroup, PermutationGroup_generic
from sage.groups.perm_gps.permgroup_named import SymmetricGroup
from sage.libs.gap.libgap import libgap
from sage.combinat.free_module import CombinatorialFreeModule

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

class SubgroupStore():
    def __init__(self):
        r"""
        This class caches subgroup information and provides
        helper methods to handle subgroup <-> name associations.

        TESTS::

            sage: P = PolynomialMolecularDecomposition()
            sage: Z3 = CyclicPermutationGroup(3)
            sage: P[Z3]
            PMD[(3, ((1,2,3),))]
            sage: P._indices.set_name(Z3, "C3")
            sage: P[Z3]
            C3
            sage: P._indices.get_name(Z3)
            'C3'
            sage: P._indices.unset_name(Z3)
            sage: P[Z3]
            PMD[(3, ((1,2,3),))]
            sage: P._indices.get_name(Z3)
            '(3, ((1,2,3),))'
        """
        self._cache = dict() # invariant to subgroups
        self._names = dict() # stores subgroup names

    def _group_invariant(self, H):
        r"""
        Returns the set of computed group invariants
        associated with a subgroup H.

        TESTS::

            sage: G = SymmetricGroup(3)
            sage: B = BurnsideRing(G)
            sage: Z3 = CyclicPermutationGroup(3)
            sage: B._indices._group_invariant(Z3)
            3
        """
        return H.order()

    def _normalize(self, H):
        r"""
        Converts a subgroup into its canonical representative
        by finding a representative of its conjugacy class, or
        using the group itself if the conjugacy class didn't exist
        before.

        TESTS::

            sage: P = PolynomialMolecularDecomposition()
            sage: H1 = PermutationGroup([(1,2),(3,4)])
            sage: H2 = PermutationGroup([(1,3),(2,4)])
            sage: P[H1]
            PMD[(4, ((3,4), (1,2)))]
            sage: P[H2]
            PMD[(4, ((3,4), (1,2)))]

            sage: G = SymmetricGroup(4)
            sage: B = BurnsideRing(G)
            sage: H1 = PermutationGroup([(1,2),(3,4)])
            sage: H2 = PermutationGroup([(1,3),(2,4)])
            sage: B[H1]
            B[((3,4), (1,2))]
            sage: B[H2]
            B[((3,4), (1,2))]
        """
        # H is of type self.element_class
        G = H.subgroup_of()
        p = self._group_invariant(H._C)
        if p in self._cache:
            for H0 in self._cache[p]:
                if _is_conjugate(G, H._C, H0._C):
                    return H0
            else:
                g = H._C.gens_small()
                H._C = G.subgroup(g)
                self._cache[p].append(H)
        else:
            g = H._C.gens_small()
            H._C = G.subgroup(g)
            self._cache[p] = [H]
        return H

    def get_name(self, H):
        r"""
        Takes a subgroup as input and returns its associated name, if any.
        Otherwise, the generators are returned. Returns a string.

        TESTS::

            sage: P = PolynomialMolecularDecomposition()
            sage: Z3 = CyclicPermutationGroup(3)
            sage: P[Z3]
            PMD[(3, ((1,2,3),))]
            sage: P._indices.set_name(Z3, "C3")
            sage: P[Z3]
            C3
            sage: P._indices.get_name(Z3)
            'C3'
            sage: P._indices.unset_name(Z3)
        """
        key = self.element_class(self, H)
        G = self._normalize(key)
        name = self._names.get(G, None)
        return name if name else repr(G)
    
    def set_name(self, H, name):
        r"""
        Takes a subgroup as input and sets its name.

        TESTS::

            sage: P = PolynomialMolecularDecomposition()
            sage: Z3 = CyclicPermutationGroup(3)
            sage: P[Z3]
            PMD[(3, ((1,2,3),))]
            sage: P._indices.set_name(Z3, "C3")
            sage: P[Z3]
            C3
            sage: P._indices.get_name(Z3)
            'C3'
            sage: P._indices.unset_name(Z3)
        """
        if not isinstance(name, str):
            raise ValueError("name must be a string")
        key = self.element_class(self, H)
        G = self._normalize(key)
        self._names[G] = name
    
    def unset_name(self, H):
        r"""
        Takes a subgroup as input and removes its name, if any.

        TESTS::

            sage: P = PolynomialMolecularDecomposition()
            sage: Z3 = CyclicPermutationGroup(3)
            sage: P[Z3]
            PMD[(3, ((1,2,3),))]
            sage: P._indices.set_name(Z3, "C3")
            sage: P[Z3]
            C3
            sage: P._indices.get_name(Z3)
            'C3'
            sage: P._indices.unset_name(Z3)
            sage: P[Z3]
            PMD[(3, ((1,2,3),))]
        """
        key = self.element_class(self, H)
        G = self._normalize(key)
        self._names.pop(G, None)

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
            ((1,2,3),)
        """
        name = self.parent()._names.get(self._C, None)
        return name if name else repr(self._C.gens())

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

class ConjugacyClassesOfSubgroups(Parent, SubgroupStore):
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
        SubgroupStore.__init__(self)

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
            ((1,2,3,4),)
        """
        if x.is_subgroup(self._G):
            key = self.element_class(self, x)
            return self._normalize(key)
        raise ValueError(f"unable to convert {x} into {self}: not a subgroup of {self._G}")

    def __iter__(self):
        r"""
        Return iterator over conjugacy classes of subgroups of the group.

        TESTS::

            sage: G = SymmetricGroup(3)
            sage: B = BurnsideRing(G)
            sage: [g for g in B._indices]
            [((),), ((1,2),), ((1,2,3),), ((1,2,3), (1,2))]
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

class ConjugacyClassOfSubgroups_SymmetricGroup(ConjugacyClassOfSubgroups):
    def __init__(self, parent, C):
        r"""
        Initialize the conjugacy class of ``C`` in SymmetricGroup(n).

        TESTS::

            sage: P = PolynomialMolecularDecomposition()
            sage: Z4 = CyclicPermutationGroup(4)
            sage: TestSuite(P(Z4)).run()
        """
        ConjugacyClassOfSubgroups.__init__(self, parent, C)
    
    @cached_method
    def subgroup_of(self):
        r"""
        Return the symmetric group which this conjugacy class
        of subgroups belongs to.

        TESTS::

            sage: from sage.rings.burnside import ConjugacyClassesOfSubgroups_SymmetricGroup
            sage: C = ConjugacyClassesOfSubgroups_SymmetricGroup(3)
            sage: Z3 = CyclicPermutationGroup(3)
            sage: CZ3 = C(Z3)
            sage: CZ3.subgroup_of()
            Symmetric group of order 3! as a permutation group
            sage: Z2 = CyclicPermutationGroup(2)
            sage: CZ2 = C(Z2)
            sage: CZ2.subgroup_of()
            Symmetric group of order 3! as a permutation group
        """
        return SymmetricGroup(self._C.degree())

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        TESTS::

            sage: P = PolynomialMolecularDecomposition()
            sage: H = PermutationGroup([[(1,2),(5,6)],[(3,4)]])
            sage: P._indices(H)
            (6, ((3,4), (1,2)(5,6)))
        """
        return f"({self.grade()}, {super()._repr_()})"

    @cached_method
    def grade(self):
        r"""
        Return the degree of this subgroup (which is the degree
        of the symmetric group it is a subgroup of).

        TESTS::

            sage: P = PolynomialMolecularDecomposition()
            sage: D2 = DiCyclicGroup(2)
            sage: P._indices(D2).grade()
            8
        """
        return self._C.degree()

    def __hash__(self):
        r"""
        Return the hash of the representative of the conjugacy class.

        TESTS::

            sage: P = PolynomialMolecularDecomposition()
            sage: H1 = P(PermutationGroup([(1,2)],domain=[1,2,3]))
            sage: H2 = P(PermutationGroup([(2,3)],domain=[1,2,3]))
            sage: hash(H1) == hash(H2)
            True
        """
        return hash((hash(SymmetricGroup(self.grade())), hash(self._C)))

    def __eq__(self, other):
        r"""
        Return if this element is equal to ``other``.

        Two elements compare equal if they are conjugate subgroups in the parent group.

        TESTS::

            sage: P = PolynomialMolecularDecomposition()
            sage: H1 = PermutationGroup([(1,2)],domain=[1,2,3])
            sage: H2 = PermutationGroup([(2,3)],domain=[1,2,3])
            sage: P[H1] == P[H2]
            True
        """
        return (isinstance(other, ConjugacyClassOfSubgroups_SymmetricGroup)
                and self.grade() == other.grade()
                and _is_conjugate(self.subgroup_of(), self._C, other._C))

    def __le__(self, other):
        r"""
        Return if this element is less than or equal to ``other``.

        If ``self and ``other`` belong to the same symmetric group, then
        ``self`` is less than or equal to ``other`` if it is conjugate to
        a subgroup of ``other`` in the parent group.

        Otherwise, ``self`` is less than ``other`` if the degree of ``self``
        is less than the degree of ``other``.
        """
        return (isinstance(other, ConjugacyClassOfSubgroups_SymmetricGroup)
                and (self.grade() < other.grade() or
                     (GAP_FAIL != libgap.ContainedConjugates(self.subgroup_of(),
                                                            other._C, self._C,
                                                            True))))

    def __lt__(self, other):
        r"""
        Return if this element is less than ``other``.
        """
        return self != other and self <= other

class ConjugacyClassesOfSubgroups_SymmetricGroup(ConjugacyClassesOfSubgroups, SubgroupStore):
    def __init__(self, n):
        r"""
        Initialize the set of conjugacy classes of SymmetricGroup(n).

        INPUT:

        ``n`` -- the degree of the symmetric group.

        TESTS::

            sage: from sage.rings.burnside import ConjugacyClassesOfSubgroups_SymmetricGroup
            sage: C = ConjugacyClassesOfSubgroups_SymmetricGroup(4)
            sage: TestSuite(C).run()
        """
        ConjugacyClassesOfSubgroups.__init__(self, SymmetricGroup(n))

    def _element_constructor_(self, x):
        r"""
        Construct the conjugacy class of subgroups containing ``x``.

        TESTS::

            sage: from sage.rings.burnside import ConjugacyClassesOfSubgroups_SymmetricGroup
            sage: C = ConjugacyClassesOfSubgroups_SymmetricGroup(4)
            sage: Z4 = CyclicPermutationGroup(4)
            sage: C(Z4)
            (4, ((1,2,3,4),))
        """
        if x.is_subgroup(self._G):
            key = self.element_class(self, self._G.subgroup(x))
            return self._normalize(key)
        raise ValueError(f"unable to convert {x} into {self}: not a subgroup of {self._G}")

    Element = ConjugacyClassOfSubgroups_SymmetricGroup

class ConjugacyClassesOfSubgroups_SymmetricGroup_all(UniqueRepresentation, Parent, SubgroupStore):
    def __init__(self):
        r"""
        Initialize the set of conjugacy classes of all symmetric groups.

        This is a graded set graded by the non-negative integers. The homogeneous piece
        with grade n is the set of conjugacy classes of subgroups of SymmetricGroup(n).

        TESTS::

            sage: from sage.rings.burnside import ConjugacyClassesOfSubgroups_SymmetricGroup_all
            sage: C = ConjugacyClassesOfSubgroups_SymmetricGroup_all()
            sage: TestSuite(C).run()
        """
        category = SetsWithGrading().Infinite()
        Parent.__init__(self, category=category)
        SubgroupStore.__init__(self)

    @cached_method
    def an_element(self):
        return self.element_class(self, SymmetricGroup(0))

    def _group_invariant(self, H):
        r"""
        Returns the set of computed group invariants
        associated with a subgroup H.

        TESTS::

            sage: P = PolynomialMolecularDecomposition()
            sage: S3 = SymmetricGroup(3)
            sage: P._indices._group_invariant(S3)
            (6, 3)
        """
        return (H.order(), H.degree())

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        TESTS::

            sage: P = PolynomialMolecularDecomposition()
            sage: P._indices
            Conjugacy classes of subgroups of symmetric groups
        """
        return "Conjugacy classes of subgroups of symmetric groups"

    def subset(self, n):
        r"""
        Returns the conjugacy classes of subgroups of SymmetricGroup(n).

        TESTS::

            sage: P = PolynomialMolecularDecomposition()
            sage: P._indices.subset(3)
            Conjugacy classes of subgroups of Symmetric group of order 3! as a permutation group
        """
        return ConjugacyClassesOfSubgroups_SymmetricGroup(n)

    def _element_constructor_(self, x):
        r"""
        Construct the conjugacy class of subgroups containing ``x``.

        TESTS::

            sage: from sage.rings.burnside import ConjugacyClassesOfSubgroups_SymmetricGroup_all
            sage: C = ConjugacyClassesOfSubgroups_SymmetricGroup_all()
            sage: Z4 = CyclicPermutationGroup(4)
            sage: C(Z4)
            (4, ((1,2,3,4),))
            sage: S3 = SymmetricGroup(3)
            sage: C(S3)
            (3, ((1,2,3), (1,2)))
        """
        if parent(x) == self:
            return x

        G = SymmetricGroup(x.degree())
        if x.is_subgroup(G):
            key = self.element_class(self, x)
            return self._normalize(key)
        raise ValueError(f"unable to convert {x} into {self}: not a subgroup of {G}")

    def __iter__(self):
        r"""
        Return iterator over conjugacy classes of subgroups of symmetric groups.

        TESTS::

            sage: P = PolynomialMolecularDecomposition()
            sage: it = iter(P._indices)
            sage: [next(it) for _ in range(1,7)]
            [(0, ((),)), (1, ((),)), (2, ((),)), (2, ((1,2),)), (3, ((),)), (3, ((2,3),))]
        """
        n = 0
        while True:
            yield from self.subset(n)
            n += 1

    def __contains__(self, H):
        r"""
        Returns if H is a subgroup of a symmetric group.

        TESTS::

            sage: P = PolynomialMolecularDecomposition()
            sage: Z3 = CyclicPermutationGroup(3)
            sage: Z3 in P._indices
            True
            sage: SH = PermutationGroup([[(1,3)]],domain=[1,3])
            sage: SH in P._indices
            False
            sage: SH2 = PermutationGroup([[(1,3)]])
            sage: SH2 in P._indices
            True
        """
        return H in self.subset(H.degree())

    Element = ConjugacyClassOfSubgroups_SymmetricGroup

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
            {48: [((1,2)(4,5,6), (3,5,4,6))], 720: [((1,2,3,4,5,6), (1,2))]}
            sage: b^2
            B[((4,5), (4,6))] + B[((5,6), (3,4), (1,2))] + B[((1,2)(4,5,6), (3,5,4,6))]
            sage: B.basis().keys()._cache
            {6: [((4,5), (4,6))],
            8: [((5,6), (3,4), (1,2))],
            48: [((1,2)(4,5,6), (3,5,4,6))],
            720: [((1,2,3,4,5,6), (1,2))]}

            sage: G = SymmetricGroup(4)
            sage: B = BurnsideRing(G)
            sage: X = Subsets(4, 2)
            sage: b = B.construct_from_action(lambda g, x: X([g(e) for e in x]), X)
            sage: b.tensor(b)
            B[((3,4), (1,2))] # B[((3,4), (1,2))]
            sage: (b.tensor(b))^2
            B[((),)] # B[((),)] + 2*B[((),)] # B[((3,4), (1,2))] + 2*B[((3,4), (1,2))] # B[((),)] + 4*B[((3,4), (1,2))] # B[((3,4), (1,2))]

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
        self._print_options['names'] = self._indices._names
        self._indices.set_name(G, "1")

    def __getitem__(self, H):
        r"""
        Return the basis element indexed by ``H``.

        ``H`` must be a subgroup of the group.

        EXAMPLES::

            sage: G = SymmetricGroup(4)
            sage: B = BurnsideRing(G)
            sage: Z4 = CyclicPermutationGroup(4)
            sage: B[Z4]
            B[((1,2,3,4),)]
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
            B[((3,4), (1,2))]

        Next, we create a group action of `S_4` on itself via conjugation::

            sage: X = G
            sage: a = lambda g, x: g*x*g.inverse()
            sage: B.construct_from_action(a, X)
            B[((3,4), (1,3)(2,4), (1,4)(2,3))] + B[((2,4,3),)] + B[((3,4), (1,2))] + B[((1,2,3,4),)] + 1

        TESTS::

            sage: G = SymmetricGroup(4)
            sage: B = BurnsideRing(G)
            sage: [b.support()[0]._C.order() for b in B.gens()]
            [1, 2, 2, 3, 4, 4, 4, 6, 8, 12, 24]
            sage: sorted((o, len(l)) for o, l in B._indices._cache.items())
            [(1, 1), (2, 2), (3, 1), (4, 3), (6, 1), (8, 1), (12, 1), (24, 1)]

            sage: G = SymmetricGroup(4)
            sage: B = BurnsideRing(G)
            sage: B(-3)
            -3
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
            ((1,2,3,4,5,6,7,8)(9,10,11,12,13,14,15,16), (1,9,5,13)(2,16,6,12)(3,15,7,11)(4,14,8,10))
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
            [            6*B[((),)]             3*B[((),)]             2*B[((),)]               B[((),)]]
            [            3*B[((),)] B[((),)] + B[((1,2),)]               B[((),)]            B[((1,2),)]]
            [            2*B[((),)]               B[((),)]        2*B[((1,2,3),)]          B[((1,2,3),)]]
            [              B[((),)]            B[((1,2),)]          B[((1,2,3),)]                      1]

        TESTS::

            sage: G = SymmetricGroup(3)
            sage: B = BurnsideRing(G)
            sage: Z3 = CyclicPermutationGroup(3)
            sage: Z2 = CyclicPermutationGroup(2)
            sage: from sage.rings.burnside import ConjugacyClassesOfSubgroups
            sage: C = ConjugacyClassesOfSubgroups(G)
            sage: B.product_on_basis(C(Z2), C(Z3))
            B[((),)]
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

class PolynomialMolecularDecomposition(CombinatorialFreeModule):
    def __init__(self, base_ring=ZZ):
        basis_keys = ConjugacyClassesOfSubgroups_SymmetricGroup_all()
        category = GradedAlgebrasWithBasis(base_ring).Commutative()
        CombinatorialFreeModule.__init__(self, base_ring,
                                        basis_keys=basis_keys,
                                        category=category,
                                        prefix="PMD")
        self._print_options['names'] = self._indices._names

    def __getitem__(self, H):
        r"""
        Return the basis element indexed by ``H``.

        ``H`` must be a subgroup of a symmetric group.

        EXAMPLES::

            sage: P = PolynomialMolecularDecomposition()
            sage: Z4 = CyclicPermutationGroup(4)
            sage: P[Z4]
            PMD[(4, ((1,2,3,4),))]
        """
        return self._from_dict({self._indices(H): 1})

    @cached_method
    def one_basis(self):
        r"""
        Returns SymmetricGroup(0), which indexes the one of this algebra,
        as per :meth:`AlgebrasWithBasis.ParentMethods.one_basis`.

        EXAMPLES::

            sage: P = PolynomialMolecularDecomposition()
            sage: P.one_basis()
            (0, ((),))
        """
        return self._indices(SymmetricGroup(0))

    # Remember, a basis element here is a molecular species.
    # When two basis elements are multiplied, you get another
    # molecular species, ie a basis element.
    # Any molecular species centered on cardinality n M_n,
    # is equivalent to [S_n/H] where H is some conjugacy class
    # of subgroups of S_n.

    def product_on_basis(self, H, K):
        r"""
        Return the product of the basis elements indexed by ``H`` and ``K``.

        Let `H` be a subgroup of `\mathfrak{G}_n` and `K` be a subgroup of `\mathfrak{G}_m`.
        Then we define their Cauchy product as the subgroup of `\mathfrak{G}_{n+m}` given by
        `H \ast K = \{h \dot k \vert h \in H_{\{1 \ldots n\}}, K_{\{n+1 \ldots n+m\}}\}` where
        the subscripts denote the domains on which H and K act. Note that this is isomorphic to
        the direct product of `H` and `K`.

        EXAMPLES::

            sage: P = PolynomialMolecularDecomposition()
            sage: matrix([[P.product_on_basis(x,y) for x in P._indices.subset(3)] for y in P._indices.subset(2)])
            [                  PMD[(5, ((),))]                PMD[(5, ((2,3),))]              PMD[(5, ((1,2,3),))]        PMD[(5, ((1,3,2), (2,3)))]]
            [               PMD[(5, ((2,3),))]          PMD[(5, ((4,5), (2,3)))]        PMD[(5, ((1,2,3), (4,5)))] PMD[(5, ((1,3,2), (4,5), (2,3)))]]
        """        
        n, m = H.grade(), K.grade()
        # There is no way to create SymmetricGroup(0) using the
        # PermutationGroup constructor as used here, so a special
        # case has to be added.
        if n+m == 0:
            return self._from_dict({self._indices(SymmetricGroup(0)): 1})
        # We only really need to multiply generators, since we are multiplying
        # permutations acting on disjoint domains.
        H_action = lambda g, e: g(e) if e<=n else e
        H_extend = PermutationGroup(H._C.gens_small(), action=H_action, domain=range(1,n+m+1))
        K_action = lambda g, e: n + g(e-n) if e>n else e
        K_restrict = PermutationGroup(K._C.gens_small(), action=K_action, domain=range(1,n+m+1))
        # We need to add the identity elements to the generating sets
        # to obtain a generating set for H*K.
        G = PermutationGroup(H_extend.gens_small()
                             + K_restrict.gens_small()
                             + [H_extend.identity()], domain=range(1,n+m+1))
        return self._from_dict({self._indices(G): 1})

    def degree_on_basis(self, H):
        r"""
        Return the degree of the basis element indexed by ``H``.

        H is an instance of ConjugacyClassOfSubgroups_SymmetricGroup.

        EXAMPLES::

            sage: P = PolynomialMolecularDecomposition()
            sage: Z4 = CyclicPermutationGroup(4)
            sage: P.degree_on_basis(P._indices(Z4))
            4
        """
        return H.grade()

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: P = PolynomialMolecularDecomposition()
            sage: P
            Polynomial Molecular Decomposition
        """
        return "Polynomial Molecular Decomposition"
