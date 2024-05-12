from sage.structure.element import Element

def _is_conjugate(G, H1, H2):
    return gap("fail") != gap.RepresentativeAction(G, H1, H2)

class BurnsideRingElement(Element):
    def __init__(self, parent, H):
        """
        Initialize an element.

        INPUT:

        * ``H`` - a formal sum of representatives of conjugacy #
          classes of subgroups

        EXAMPLES::

            sage: G = SymmetricGroup(4)
            sage: B = BurnsideRing(G)
            sage: B(G)
        """
        Element.__init__(self, parent)
        self._H = H

    def _mul_(self, right):
        """
        Return the product of ``self``  and ``right``.

        For the symmetric group, this is also known as the Hadamard
        or tensor product of group actions.
        """
        pass

    def _add_(self, right):
        """
        Return the sum of ``self``  and ``right``.
        """
        pass


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

        # given a group action
        if x is None:
            H = PermutationGroup(self._G.gens(), action=a, domain=domain)
            # decompose H into a sum F of conjugacy classes
            # normalize each summand
            return self.element_class(self, F)

        # if x is a single subgroup of self._G
        return self.element_class(self, FormalSum([(1, normalize(x))]))

class PolynomialSpeciesElement(Parent):

    def _mul_(self, right):
        """
        Return the Cauchy product of ``self``  and ``right``.

        EXAMPLES::

            sage: X * X
            X^2
        """
        pass

class PolynomialSpecies(Parent):
    pass

