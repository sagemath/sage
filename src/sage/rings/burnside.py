from sage.misc.cachefunc import cached_method
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.richcmp import op_EQ, op_NE
from sage.structure.parent import Parent
from sage.structure.element import Element
from sage.rings.integer_ring import ZZ
from sage.structure.formal_sum import FormalSum
from sage.categories.sets_cat import cartesian_product, Sets
from sage.groups.perm_gps.permgroup import PermutationGroup
from sage.libs.gap.libgap import libgap
from sage.misc.repr import repr_lincomb
from sage.categories.algebras import Algebras
from sage.combinat.free_module import CombinatorialFreeModule

def _is_conjugate(G, H1, H2):
    return libgap.eval('fail') != libgap.RepresentativeAction(G, H1, H2)

class ConjugacyClassOfSubgroups(Element):
    def __init__(self, parent, C):
        Element.__init__(self, parent)
        self._C = C

    def __hash__(self):
        return hash(self._C)

    def _repr_(self):
        return repr(self._C.gens())

    def __le__(self, other):
        return libgap.eval('fail') != libgap.ContainedConjugates(self.parent()._G, other._C, self._C, True)

class ConjugacyClassesOfSubgroups(Parent):
    def __init__(self, G):
        self._G = G
        self._cache = dict() # invariant to subgroups
        Parent.__init__(self, category=Sets().Finite())

    def _group_invariant(self, H):
        return H.order()

    def _element_constructor_(self, x):
        def normalize(H):
            p = self._group_invariant(H)
            if p in self._cache:
                for H0 in self._cache[p]:
                    if _is_conjugate(self._G, H, H0):
                        return H0
                else:
                    self._cache[p].append(H)
            else:
                self._cache[p] = [H]
            return H

        if x.is_subgroup(self._G):
            return self.element_class(self, normalize(x))
        else:
            raise ValueError("x must be a subgroup of " + repr(self._G))

        raise ValueError(f"unable to convert {x} into {self}")

    def __iter__(self):
        return iter(self(H) for H in self._G.conjugacy_classes_subgroups())

    Element = ConjugacyClassOfSubgroups

class BurnsideRingElement(Element):
    def __init__(self, parent, F):
        r"""
        Initialize an element.

        INPUT:

        - ``F`` - a formal sum of representatives of conjugacy #
          classes of subgroups

        EXAMPLES::

            sage: G = SymmetricGroup(4)
            sage: B = BurnsideRing(G)
            sage: B(G)
            1

            sage: X = Subsets(4, 2)
            sage: a = lambda g, x: X([g(e) for e in x])
            sage: B(domain=X, action=a)
            [(3,4), (1,2), (1,2)(3,4)]
        """
        Element.__init__(self, parent)
        self._F = F

    def _repr_(self):
        r"""
        Return a string representation of ``self``.
        """
        return repr_lincomb(([(H.gens() if name is None else name, c)
                              for c, (name, H) in self._F]),
                            repr_monomial=lambda s: s if isinstance(s, str) else repr(list(s)))

    def _richcmp_(self, other, op):
        """

        sage: G = SymmetricGroup(4)
        sage: B = BurnsideRing(G)
        sage: B.zero() == 0
        True
        """
        if op is op_EQ:
            return not bool(self._F - other._F)
        if op is op_NE:
            return bool(self._F - other._F)

    def _acted_upon_(self, scalar, self_on_left):
        r"""
        Scalar multiplication for ``self`` by ``scalar``.

        INPUT:

        - ``scalar`` -- an element of the base ring
        - ``self_on_left`` -- boolean; if ``True``, compute ``self * scalar``

        EXAMPLES::

            sage: G = SymmetricGroup(4)
            sage: B = BurnsideRing(G)
            sage: b = B(G)
            sage: b
            1
            sage: 2*b
            2*1
        """
        B = self.parent()
        F = scalar * self._F
        return B.element_class(B, F)

    def _mul_(self, right):
        r"""
        Return the product of ``self``  and ``right``.

        For the symmetric group, this is also known as the Hadamard
        or tensor product of group actions.

        EXAMPLES::

            sage: G = SymmetricGroup(3)
            sage: B = BurnsideRing(G)
            sage: matrix([[b * c for b in B.gens()] for c in B.gens()])
            [        6*[()]         3*[()]         2*[()]           [()]]
            [        3*[()] [(2,3)] + [()]           [()]        [(2,3)]]
            [        2*[()]           [()]    2*[(1,2,3)]      [(1,2,3)]]
            [          [()]        [(2,3)]      [(1,2,3)]              1]
        """
        #TODO: Find faster way to multiply
        assert right.parent() == self.parent()
        B = self.parent()
        G = B._G

        def mul(H1, H2):
            dom1 = [frozenset(g) for g in G.cosets(H1, side="left")]
            dom2 = [frozenset(g) for g in G.cosets(H2, side="left")]
            domain = cartesian_product([dom1, dom2])

            def action(g, pair):
                return (frozenset(g * h for h in pair[0]),
                        frozenset(g * h for h in pair[1]))

            return B(action=action, domain=domain)._F

        result = FormalSum(0)
        for c1, (_, g1) in self._F:
            for c2, (_, g2) in right._F:
                result += c1 * c2 * mul(g1, g2)

        return B.element_class(B, result)

    def _add_(self, right):
        r"""
        Return the sum of ``self``  and ``right``.
        """
        P = self.parent()
        F = self._F + right._F
        return P.element_class(P, F)

class BurnsideRing(CombinatorialFreeModule):
    def __init__(self, G, base_ring=ZZ):
        """
            sage: G = SymmetricGroup(4)
            sage: B = BurnsideRing(G)
            sage: TestSuite(B).run()
        """
        self._G = G
        self._cache = dict() # invariant to a list of pairs (name, subgroup)
        basis = ConjugacyClassesOfSubgroups(G)
        category = Algebras(base_ring).Commutative()
        CombinatorialFreeModule.__init__(self, base_ring, basis,
                                        element_class=BurnsideRingElement,
                                        category=category)

    def _group_invariant(self, H):
        return H.order()

    def _coerce_map_from_(self, S):
        """
        Return ``True`` if a coercion from ``S`` exists.
        """
        if self.base_ring().has_coerce_map_from(S):
            return True
        return None

    @cached_method
    def _an_element_(self):
        return self.one()

    Element = BurnsideRingElement

    def _element_constructor_(self, x=None, action=None, domain=None, name=None):
        r"""
        Construct an element of the Burnside ring.

        INPUT:

        - ``x`` - data for an element
        - ``action`` - an action on ``domain``
        - ``domain`` - a finite set

        ``x`` can be a subgroup of `G` or a formal sum of such
        subgroups.

        EXAMPLES:

            sage: G = SymmetricGroup(4)
            sage: B = BurnsideRing(G)

        We create a group action of `S_4` on two-element subsets::

            sage: X = Subsets(4, 2)
            sage: a = lambda g, x: X([g(e) for e in x])
            sage: B(domain=X, action=a)
            [(3,4), (1,2), (1,2)(3,4)]

            sage: X = G
            sage: a = lambda g, x: g*x*g.inverse()
            sage: B(domain=X, action=a)
            1 + [(2,4), (1,4)(2,3)] + [(2,3,4)] + [(3,4), (1,2), (1,2)(3,4)] + [(1,3,2,4)]

        TESTS::

            sage: G = SymmetricGroup(4)
            sage: B = BurnsideRing(G)
            sage: [H._F[0][1][1].order() for H in B.gens()]
            [1, 2, 2, 3, 4, 4, 4, 6, 8, 12, 24]
            sage: sorted((o, len(l)) for o, l in B._cache.items())
            [(1, 1), (2, 2), (3, 1), (4, 3), (6, 1), (8, 1), (12, 1), (24, 1)]

            sage: G = SymmetricGroup(4)
            sage: B = BurnsideRing(G)
            sage: B(-3)
            -3*1
        """
        if action is None and domain is None and x in self.base_ring():
            return x * self.one()

        def normalize(H, name=None):
            p = self._group_invariant(H)
            if p in self._cache:
                for name0, H0 in self._cache[p]:
                    if _is_conjugate(self._G, H, H0):
                        assert name is None or name0 == name, "the name should not change"
                        return name0, H0
                else:
                    self._cache[p].append((name, H))
            else:
                self._cache[p] = [(name, H)]
            return name, H

        def find_stabilizer(action, pnt):
            stabilizer = []
            for g in self._G:
                if action(g, pnt) == pnt:
                    stabilizer.append(g)
            H = self._G.subgroup(stabilizer)
            gens = H.gens_small()
            return self._G.subgroup(gens)

        # given a group action
        if action is not None and domain is not None:
            assert name is None
            H = PermutationGroup(self._G.gens(), action=action, domain=domain)
            # decompose H into a sum F of conjugacy classes
            orbit_list = H.orbits()
            # find the stabilizer subgroups
            # TODO: find a way to do this with GAP instead
            stabilizer_list = [find_stabilizer(action, orbit[0]) for orbit in orbit_list]
            # normalize each summand and collect terms
            from collections import Counter
            C = Counter([normalize(stabilizer) for stabilizer in stabilizer_list])
            # create formal sum
            F = FormalSum([(coeff, subgroup) for subgroup, coeff in C.items()])
            return self.element_class(self, F)
        elif action is not None and domain is None:
            raise ValueError("If action is provided then domain must be provided")
        elif action is None and domain is not None:
            raise ValueError("If domain is provided then action must be provided")

        if isinstance(x, list) or isinstance(x, FormalSum):
            assert name is None
            # if x is a list of pairs (coeff, subgroup) or FormalSum
            if not all([subgroup.is_subgroup(self._G) for coeff, subgroup in x]):
                raise ValueError("All groups in list must be subgroups of " + repr(self._G))
            return self.element_class(self, FormalSum([(coeff, normalize(subgroup)) for coeff, subgroup in x]))
        elif x.is_subgroup(self._G):
            # if x is a single subgroup of self._G
            return self.element_class(self, FormalSum([(1, normalize(x, name))]))

        raise ValueError(f"unable to convert {x} into {self}")

    @cached_method
    def one(self):
        return self._element_constructor_(self._G, name="1")

    @cached_method
    def zero(self):
        """
        EXAMPLES::

            sage: G = SymmetricGroup(3)
            sage: B = BurnsideRing(G)
            sage: (0 * B.one()).is_zero()
            True
        """
        return self.element_class(self, FormalSum([]))

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

    @cached_method
    def gens(self):
        r"""
        Return the generators of ``self``.

        These are the conjugacy classes of subgroups of the
        underlying group :meth:`group`.

        EXAMPLES::

            sage: G = SymmetricGroup(3)
            sage: B = BurnsideRing(G)
            sage: B.gens()
            ([()], [(2,3)], [(1,2,3)], 1)
        """
        return tuple(self(H) for H in self._G.conjugacy_classes_subgroups())

    def _repr_(self):
        r"""
        Return a string representation of ``self``.
        """
        return "Burnside ring of " + repr(self._G)
