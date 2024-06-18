from sage.misc.cachefunc import cached_method
from sage.structure.parent import Parent
from sage.structure.element import Element
from sage.structure.element import parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.rings.integer_ring import ZZ
from sage.categories.sets_cat import cartesian_product
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.categories.sets_with_grading import SetsWithGrading
from sage.groups.perm_gps.permgroup import PermutationGroup, PermutationGroup_generic
from sage.groups.perm_gps.permgroup_named import SymmetricGroup
from sage.libs.gap.libgap import libgap
from sage.categories.algebras import Algebras
from sage.combinat.free_module import CombinatorialFreeModule
from sage.categories.graded_algebras_with_basis import GradedAlgebrasWithBasis
from sage.categories.filtered_modules_with_basis import FilteredModulesWithBasis
from sage.rings.integer import Integer

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

class SubgroupNameStore():
    def __init__(self, basis_keys):
        r"""
        This class stores subgroup <-> name associations.
        """
        self._key_constructor = basis_keys
        self._cache = dict() # invariant to subgroups
        self._names = dict() # dictionary mapping subgroups to names

    def _group_invariant(self, H):
        return H.order()

    def _normalize(self, H):
        key = self._key_constructor(H)
        p = self._group_invariant(key._C)
        if p in self._cache:
            for H0 in self._cache[p]:
                if _is_conjugate(key.subgroup_of(), key._C, H0):
                    return H0
            else:
                self._cache[p].append(H)
        else:
            self._cache[p] = [H]
        return H

    def get_name(self, H):
        r"""
        Takes a subgroup as input and returns its associated name, if any.
        Otherwise, the generators are returned. Returns a string.
        """
        # If name exists, I want to return that, otherwise I want a default name.
        # But only basis_keys can provide that default name.
        G = self._normalize(H)
        name = self._names.get(G, None)
        return name if name else self._key_constructor(H).default_name()
    
    def set_name(self, H, name):
        r"""
        Takes a subgroup as input and sets its name.
        """
        if not isinstance(name, str):
            raise ValueError("name must be a string")
        G = self._normalize(H)
        self._names[G] = name
    
    def unset_name(self, H):
        r"""
        Takes a subgroup as input and removes its name, if any.
        """
        G = self._normalize(H)
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
        return self.parent()._G

    def default_name(self):
        return repr(self._C.gens())

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

            sage: G = SymmetricGroup(4)
            sage: B = BurnsideRing(G)
            sage: Z4 = CyclicPermutationGroup(4)
            sage: B[Z4]
            B[((1,2,3,4),)]
            sage: B.set_name(Z4, "Z4")
            sage: B[Z4] #indirect doctest
            B[Z4]
            sage: B.unset_name(Z4)
            sage: B[Z4]
            B[((1,2,3,4),)]
        """
        namefunc = self.parent()._namegetter
        return namefunc(self._C) if namefunc else self.default_name()

    def __le__(self, other):
        r"""
        Return if this element is less or equal to ``other``.

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

class ConjugacyClassesOfSubgroups(Parent):
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
        self._namegetter = None
        Parent.__init__(self, category=FiniteEnumeratedSets())

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
        Return the hash of the representative of the conjugacy class.

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
            sage: C(Z4) #indirect doctest
            ((1,2,3,4),)
        """
        if x.is_subgroup(self._G):
            return self.element_class(self, x)
        raise ValueError(f"unable to convert {x} into {self}: not a subgroup of {self._G}")

    def __iter__(self):
        r"""
        Return iterator over conjugacy classes of subgroups of the group.

        TESTS::

            sage: G = SymmetricGroup(3)
            sage: B = BurnsideRing(G)
            sage: [g for g in B._indices]
            [((),), ((1,2),), ((1,2,3),), 1]
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

            sage: G = SymmetricGroup(4)
            sage: B = BurnsideRing(G)
            sage: Z4 = CyclicPermutationGroup(4)
            sage: TestSuite(B(Z4)).run()
        """
        ConjugacyClassOfSubgroups.__init__(self, parent, C)
    
    def subgroup_of(self):
        return SymmetricGroup(self.grade())

    def default_name(self):
        return f"({self.grade()}, {repr(self._C.gens())})"

    def grade(self):
        return self._C.degree()

    def __hash__(self):
        r"""
        Return the hash of the representative of the conjugacy class.
        """
        return hash((hash(SymmetricGroup(self.grade())), hash(self._C)))

    def __eq__(self, other):
        r"""
        Return if this element is equal to ``other``.

        Two elements compare equal if they are conjugate subgroups in the parent group.
        """
        return (isinstance(other, ConjugacyClassOfSubgroups_SymmetricGroup)
                and self.grade() == other.grade() and _is_conjugate(SymmetricGroup(self._C), self._C, other._C))

class ConjugacyClassesOfSubgroups_SymmetricGroup(ConjugacyClassesOfSubgroups):
    def __init__(self, n):
        ConjugacyClassesOfSubgroups.__init__(self, SymmetricGroup(n))

    Element = ConjugacyClassOfSubgroups_SymmetricGroup

class ConjugacyClassesOfSubgroups_SymmetricGroup_all(UniqueRepresentation, Parent):
    def __init__(self):
        self._namegetter = None
        category = SetsWithGrading().Infinite()
        Parent.__init__(self, category=category)

    def _repr_(self):
        return "Conjugacy classes of subgroups of symmetric groups"

    def subset(self, size):
        return ConjugacyClassesOfSubgroups_SymmetricGroup(size)

    def _element_constructor_(self, H):
        r"""
        H is a subgroup of a symmetric group.
        """
        return self.element_class(self, H)

    def __iter__(self):
        n = 0
        while True:
            yield from self.subset(n)
            n += 1

    def __contains__(self, H):
        r"""
        Returns if H is a subgroup of a symmetric group.
        """
        return H in self.subset(H.degree())

    Element = ConjugacyClassOfSubgroups_SymmetricGroup

class BurnsideRing(CombinatorialFreeModule, SubgroupNameStore):
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
            {48: [Subgroup generated by [(3,5,4,6), (1,2)(4,5,6)] of (Symmetric group of order 6! as a permutation group)],
             720: [Symmetric group of order 6! as a permutation group]}
            sage: b^2
            B[((3,5,4,6), (1,2)(4,5,6))] + B[((5,6), (4,6,5))] + B[((5,6), (3,4), (1,2))]
            sage: B.basis().keys()._cache
            {6: [Subgroup generated by [(5,6), (4,6,5)] of (Symmetric group of order 6! as a permutation group)],
             8: [Subgroup generated by [(5,6), (3,4), (1,2)] of (Symmetric group of order 6! as a permutation group)],
             48: [Subgroup generated by [(3,5,4,6), (1,2)(4,5,6)] of (Symmetric group of order 6! as a permutation group)],
             720: [Symmetric group of order 6! as a permutation group]}

            sage: G = SymmetricGroup(4)
            sage: B = BurnsideRing(G)
            sage: X = Subsets(4, 2)
            sage: b = B.construct_from_action(lambda g, x: X([g(e) for e in x]), X)
            sage: b.tensor(b)
            B[((3,4), (1,2)(3,4))] # B[((3,4), (1,2)(3,4))]
            sage: (b.tensor(b))^2
            4*B[((3,4), (1,2)(3,4))] # B[((3,4), (1,2)(3,4))] + 2*B[((),)] # B[((3,4), (1,2)(3,4))] + 2*B[((3,4), (1,2)(3,4))] # B[((),)] + B[((),)] # B[((),)]

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
        SubgroupNameStore.__init__(self, basis_keys)
        basis_keys._namegetter = self.get_name
        self.set_name(G, "1")

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
            B[((3,4), (1,2)(3,4))]

        Next, we create a group action of `S_4` on itself via conjugation::

            sage: X = G
            sage: a = lambda g, x: g*x*g.inverse()
            sage: B.construct_from_action(a, X)
            B[1] + B[((1,4)(2,3), (1,3)(2,4), (3,4))] + B[((2,4,3),)] + B[((3,4), (1,2)(3,4))] + B[((1,2,3,4),)]

        TESTS::

            sage: G = SymmetricGroup(4)
            sage: B = BurnsideRing(G)
            sage: [list(b.monomial_coefficients().keys())[0]._C.order() for b in B.gens()]
            [1, 2, 2, 3, 4, 4, 4, 6, 8, 12, 24]
            sage: sorted((o, len(l)) for o, l in B._indices._cache.items())
            [(1, 1), (2, 2), (3, 1), (4, 3), (6, 1), (8, 1), (12, 1), (24, 1)]

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

    def product_on_basis(self, g1, g2):
        r"""
        Return the product of the basis elements indexed by ``g1`` and ``g2``.

        For the symmetric group, this is also known as the Hadamard
        or tensor product of group actions.

        EXAMPLES::

            sage: G = SymmetricGroup(3)
            sage: B = BurnsideRing(G)
            sage: matrix([[b * c for b in B.gens()] for c in B.gens()])
            [            6*B[((),)]             3*B[((),)]             2*B[((),)]               B[((),)]]
            [            3*B[((),)] B[((2,3),)] + B[((),)]               B[((),)]            B[((2,3),)]]
            [            2*B[((),)]               B[((),)]        2*B[((1,2,3),)]          B[((1,2,3),)]]
            [              B[((),)]            B[((2,3),)]          B[((1,2,3),)]                   B[1]]

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
        #TODO: Find faster way to multiply
        assert g1.parent() == g2.parent()
        G = g1.subgroup_of()
        dom1 = [frozenset(g) for g in G.cosets(g1._C, side="left")]
        dom2 = [frozenset(g) for g in G.cosets(g2._C, side="left")]
        domain = cartesian_product([dom1, dom2])

        def action(g, pair):
            return (frozenset(g * h for h in pair[0]),
                    frozenset(g * h for h in pair[1]))

        return self.construct_from_action(action, domain)

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

class PolynomialMolecularDecomposition(CombinatorialFreeModule, SubgroupNameStore):
    def __init__(self, base_ring=ZZ):
        basis_keys = ConjugacyClassesOfSubgroups_SymmetricGroup_all()
        category = GradedAlgebrasWithBasis(base_ring)
        CombinatorialFreeModule.__init__(self, base_ring,
                                        basis_keys=basis_keys,
                                        category=category,
                                        prefix="PMD")
        SubgroupNameStore.__init__(self, basis_keys)
        basis_keys._namegetter = self.get_name

    # FIXME: this is currently required, because the implementation of ``basis``
    # in CombinatorialFreeModule overrides that of GradedModulesWithBasis
    basis = FilteredModulesWithBasis.ParentMethods.__dict__['basis']

    def __getitem__(self, x):
        r"""
        Return the basis element indexed by ``x``.
        """
        return self._from_dict({self._indices(x): 1})

    @cached_method
    def one_basis(self):
        r"""
        Returns (0, S0), which indexes the one of this algebra,
        as per :meth:`AlgebrasWithBasis.ParentMethods.one_basis`.
        """
        return self._indices(SymmetricGroup(0))

    # Remember, a basis element here is a molecular species.
    # When two basis elements are multiplied, you get another
    # molecular species, ie a basis element.
    # Any molecular species centered on cardinality n M_n,
    # is equivalent to [S_n/H] where H is some conjugacy class
    # of subgroups of S_n.

    def product_on_basis(self, g1, g2):
        n, m = g1.grade(), g2.grade()
        H, K = g1._C, g2._C

        def construct_element(h, k):
            element = [None for _ in range(n+m)]
            for i in range(n+m):
                if i < n:
                    element[i] = h(i+1)
                else:
                    element[i] = n + k(i-n+1)
            return element
        H_ast_K = [
            construct_element(h, k)
            for h in H
            for k in K
        ]
        # There is no way to create SymmetricGroup(0) using the
        # PermutationGroup constructor as used here, so a special
        # case has to be added.
        if n+m == 0:
            return self._from_dict({self._indices(SymmetricGroup(0)): 1})
        G = PermutationGroup(H_ast_K)
        return self._from_dict({self._indices(G): 1})

    def degree_on_basis(self, x):
        r"""
        x is an instance of ConjugacyClassOfSubgroups_SymmetricGroup.
        """
        return x.grade()

    def _repr_(self):
        return "Polynomial Molecular Decomposition"
