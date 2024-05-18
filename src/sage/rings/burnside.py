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

    def _normalize(self, H):
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

    def _element_constructor_(self, x):
        if x.is_subgroup(self._G):
            return self.element_class(self, self._normalize(x))
        else:
            raise ValueError("x must be a subgroup of " + repr(self._G))

        raise ValueError(f"unable to convert {x} into {self}")

    def __iter__(self):
        return iter(self(H) for H in self._G.conjugacy_classes_subgroups())

    def _repr_(self):
        r"""
        Return a string representation of ``self``.
        """
        return "Conjugacy classes of subgroups of " + repr(self._G)

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
        self._monomial_coefficients = {parent._indices(H): c for c, H in self._F}

    def _repr_(self):
        r"""
        Return a string representation of ``self``.
        """
        return repr_lincomb(([(f"B[{H.gens()}]" if self.parent()._names[H] is None else f"B[{self.parent()._names[H]}]", c)
                              for c, H in self._F]),
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

    def _add_(self, right):
        r"""
        Return the sum of ``self``  and ``right``.
        """
        P = self.parent()
        F = self._F + right._F
        return P.element_class(P, F)

    def monomial_coefficients(self, copy=True):
        if copy:
            return dict(self._monomial_coefficients)
        return self._monomial_coefficients

class BurnsideRing(CombinatorialFreeModule):
    def __init__(self, G, base_ring=ZZ):
        """
            sage: G = SymmetricGroup(4)
            sage: B = BurnsideRing(G)
            sage: TestSuite(B).run()
        """
        self._G = G
        self._cache = dict() # invariant to subgroups
        self._names = dict() # dictionary mapping subgroups to names
        self.rename_gen(G, "1") # name unit of ring as 1
        basis = ConjugacyClassesOfSubgroups(G)
        category = Algebras(base_ring).Commutative().WithBasis()
        CombinatorialFreeModule.__init__(self, base_ring, basis,
                                        element_class=BurnsideRingElement,
                                        category=category)

    def _group_invariant(self, H):
        return H.order()

    def _normalize(self, H, name=None):
        if name is not None:
            self._names[H] = name
        elif H not in self._names:
            self._names[H] = None
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

    def _element_constructor_(self, x):
        r"""
        Construct an element of the Burnside ring explicitly.

        INPUT:

        - ``x`` - data for an element

        ``x`` can be a subgroup of `G` or a formal sum of such
        subgroups.
        """
        if x in self.base_ring():
            return x * self.one()

        if isinstance(x, list) or isinstance(x, FormalSum):
            # if x is a list of pairs (coeff, subgroup) or FormalSum
            if not all([subgroup.is_subgroup(self._G) for coeff, subgroup in x]):
                raise ValueError("All groups in list must be subgroups of " + repr(self._G))
            return self.element_class(self, FormalSum([(coeff, self._normalize(subgroup)) for coeff, subgroup in x]))
        elif x.is_subgroup(self._G):
            # if x is a single subgroup of self._G
            return self.element_class(self, FormalSum([(1, self._normalize(x))]))

        raise ValueError(f"unable to convert {x} into {self}")

    def _from_dict(self, d, coerce=False, remove_zeros=True):
        r"""
        Construct an element of "self" from an "{index: coefficient}"
        dictionary.

        INPUT:

        * "d" -- a dictionary "{index: coeff}" where each "index" is the
            index of a basis element and each "coeff" belongs to the
            coefficient ring "self.base_ring()"

        * "coerce" -- a boolean (default: "False"), whether to coerce the
            coefficients "coeff" to the coefficient ring

        * "remove_zeros" -- a boolean (default: "True"), if some
            coefficients "coeff" may be zero and should therefore be removed
        """
        assert isinstance(d, dict)
        if coerce:
            R = self.base_ring()
            d = {H._C: R(coeff) for H, coeff in d.items()}
        if remove_zeros:
            d = {H._C: coeff for H, coeff in d.items() if coeff}
        l = [(coeff, H._C) for H, coeff in d.items()]
        return self._element_constructor_(l)

    def construct_from_action(self, action, domain):
        r"""
        Construct an element of the Burnside ring from a group action.

        INPUT:

        - ``action`` - an action on ``domain``
        - ``domain`` - a finite set

        EXAMPLES:

            sage: G = SymmetricGroup(4)
            sage: B = BurnsideRing(G)

        We create a group action of `S_4` on two-element subsets::

            sage: X = Subsets(4, 2)
            sage: a = lambda g, x: X([g(e) for e in x])
            sage: B(domain=X, action=a)
            [(3,4), (1,2), (1,2)(3,4)]

        Next, we create a group action of `S_4` on itself via conjugation::

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
        C = Counter([self._normalize(stabilizer) for stabilizer in stabilizer_list])
        # create formal sum
        F = FormalSum([(coeff, subgroup) for subgroup, coeff in C.items()])
        return self.element_class(self, F)

    def rename_gen(self, subgroup, name):
        r"""
        Rename conjugacy class of ``subgroup`` to ``name``.
        Passing ``None`` to ``name`` will remove any previously assigned name.
        """
        if not isinstance(name, str):
            raise TypeError("name must be a string")
        self._normalize(subgroup, name)

    def monomial(self, H):
        r"""
        Return the basis element indexed by `H`.
        """
        return self.element_class(self, FormalSum([(1, self._normalize(H))]))

    @cached_method
    def one_basis(self):
        r"""
        Returns the generators of the underlying group, which index the one
        of this algebra, as per :meth:`AlgebrasWithBasis.ParentMethods.one_basis`.
        """
        return self._indices(G)

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

    def product_on_basis(self, g1, g2):
        r"""
        Return the product of ``g1``  and ``g2``.

        For the symmetric group, this is also known as the Hadamard
        or tensor product of group actions.

        EXAMPLES::

            sage: G = SymmetricGroup(3)
            sage: B = BurnsideRing(G)
            sage: matrix([[b * c for b in B.gens()] for c in B.gens()])
            [            6*B[((),)]             3*B[((),)]             2*B[((),)]               B[((),)]]
            [            3*B[((),)] B[((1,2),)] + B[((),)]               B[((),)]            B[((1,2),)]]
            [            2*B[((),)]               B[((),)]        2*B[((1,2,3),)]          B[((1,2,3),)]]
            [              B[((),)]            B[((1,2),)]          B[((1,2,3),)]                   B[1]]
        """
        #TODO: Find faster way to multiply
        assert g1.parent() == g2.parent()
        G = g1.parent()._G
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
