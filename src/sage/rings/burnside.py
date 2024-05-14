from sage.structure.element import Element

def _is_conjugate(G, H1, H2):
    return gap("fail") != gap.RepresentativeAction(G, H1, H2)

class BurnsideRingElement(Element):
    def __init__(self, parent, F):
        r"""
        Initialize an element.

        INPUT:

        * ``F`` - a formal sum of representatives of conjugacy #
          classes of subgroups

        EXAMPLES::

            sage: G = SymmetricGroup(4)
            sage: B = BurnsideRing(G)
            sage: B(G)
            Symmetric group of order 4! as a permutation group
             as an element of the Burnside ring
        """
        Element.__init__(self, parent)
        self._F = F

    def _repr_(self):
        r"""
        Return a string representation of ``self``.
        """
        return repr(self._F) + " as an element of the Burnside ring"

    def _acted_upon_(self, scalar, self_on_left):
        r"""
        Scalar multiplication for ``self`` by ``scalar``.

        INPUT:

        - ``scalar`` -- an element of the base ring
        - ``self_on_left`` -- boolean; if ``True``, compute ``self * scalar``

        EXAMPLES:

            sage: B = BurnsideRing(G)
            sage: b = B(G)
            sage: b
            Symmetric group of order 4! as a permutation group
             as an element of the Burnside ring
            sage: 2*b
            2*Symmetric group of order 4!
             as a permutation group as an element of the Burnside ring
        """
        B = self.parent()
        F = scalar * self._F
        return B.element_class(B, F)

    def _mul_(self, right):
        r"""
        Return the product of ``self``  and ``right``.

        For the symmetric group, this is also known as the Hadamard
        or tensor product of group actions.
        """
        #TODO: Find faster way to multiply
        assert right.parent() == self.parent()
        B = self.parent()
        G = B._G
        def mul(H1, H2):
            dom1 = [frozenset(g) for g in G.cosets(H1)]
            dom2 = [frozenset(g) for g in G.cosets(H2)]
            domain = cartesian_product([dom1, dom2])
            def action(g, pair):
                return (frozenset(g*h for h in pair[0]),
                        frozenset(g*h for h in pair[1]))
            return B(action=action, domain=domain)._F
        result = FormalSum(0)
        for c1, g1 in self._F:
            for c2, g2 in right._F:
                result += c1*c2*mul(g1, g2)
        return B.element_class(B, result)

    def _add_(self, right):
        r"""
        Return the sum of ``self``  and ``right``.
        """
        P = self.parent()
        F = self._F + right._F
        return P.element_class(P, F)


class BurnsideRing(Parent):
    def __init__(self, G, base_ring=ZZ):
        self._G = G
        self._cache = dict()
        Parent.__init__(self, base=base_ring)

    def _group_invariant(self, H):
        return H.order()

    Element = BurnsideRingElement

    def _element_constructor_(self, x=None, action=None, domain=None):
        """
        Construct an element of the Burnside ring.

        INPUT:

        - ``x`` - data for an element

        ``x`` can be one of the following:

        * a subgroup of `G`
        * a formal sum of such subgroups `G`

        (* a group action on some set `X`)

        EXAMPLES:

            sage: S4 = SymmetricGroup(4)
            sage: B = BurnsideRing(S4)

        Maybe better like:

            sage: a = lambda g, x: Set([g(e) for e in x])
            sage: B(action=a, domain=Subsets(4, 2))

        We create a group action of `S_4` on two-element subsets::

            sage: H = PermutationGroup(S4.gens(), action=a, domain=Subsets(4, 2))

            # we can only trust the user that this is the homomorphic
            # image of a group action of S4
            sage: B(H)
        """
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
        
        def find_stabilizer(action, pnt):
            stabilizer = []
            for g in self._G:
                if action(g, pnt) == pnt:
                    stabilizer.append(g)
            return self._G.subgroup(stabilizer)

        # given a group action
        if action is not None and domain is not None:
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

        if isinstance(x, list) or isinstance(x,FormalSum):
            # if x is a list of pairs (coeff, subgroup) or FormalSum
            if not all([subgroup.is_subgroup(self._G) for coeff, subgroup in x]):
                raise ValueError("All groups in list must be subgroups of " + repr(self._G))
            return self.element_class(self, FormalSum([(coeff, normalize(subgroup)) for coeff, subgroup in x]))
        elif x.is_subgroup(self._G):
            # if x is a single subgroup of self._G
            return self.element_class(self, FormalSum([(1, normalize(x))]))

        raise ValueError(f"unable to convert {x} into {self}")


    def _repr_(self):
        r"""
        Return a string representation of ``self``.
        """
        return "Burnside ring of " + repr(self._G)