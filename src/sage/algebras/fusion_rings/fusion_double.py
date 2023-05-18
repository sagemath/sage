"""
The Fusion Ring of the Drinfeld Double of a Finite Group
"""
# ****************************************************************************
#  Copyright (C) 2023 Wenqi Li
#                     Daniel Bump <bump at match.stanford.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from sage.algebras.fusion_rings.generic_fusion_ring import FusionRing
from sage.misc.cachefunc import cached_method
from sage.rings.integer_ring import ZZ
from sage.sets.set import Set

class FusionDouble(FusionRing):
    r"""
    This constructs the Fusion Ring of the modular tensor category of modules over the Drinfeld
    Double of a finite group. Usage is similar to :class:`FusionRing`::.
    Since many of the methods are similar, we assume the reader is familiar with that
    class.

    INPUT:

    - ``G`` -- a finite group
    - ``prefix`` (optional: defaults to 's') a prefix for the names of simple objects.
    - ``inject_varables`` (optional): set TRUE to create variables for the simple objects.

    REFERENCES:

    - [BaKi2001]_ Chapter 3
    - [Mas1995]_
    - [CHW2015]_
    - [Goff1999]_

    EXAMPLES::

        sage: G = DihedralGroup(5)
        sage: H = FusionRing(G, inject_variables=True)
        sage: H.basis()
        Finite family {0: s0, 1: s1, 2: s2, 3: s3, 4: s4, 5: s5, 6: s6, 7: s7, 8: s8, 9: s9, 10: s10, 11: s11, 12: s12, 13: s13, 14: s14, 15: s15}
        sage: for x in H.basis():
        ....:     print ("%s : %s"%(x,x^2))
        ....:
        s0 : s0
        s1 : s0
        s2 : s0 + s1 + s3
        s3 : s0 + s1 + s2
        s4 : s0 + s2 + s3 + s6 + s7 + s8 + s9 + s10 + s11 + s12 + s13 + s14 + s15
        s5 : s0 + s2 + s3 + s6 + s7 + s8 + s9 + s10 + s11 + s12 + s13 + s14 + s15
        s6 : s0 + s1 + s11
        s7 : s0 + s1 + s13
        s8 : s0 + s1 + s15
        s9 : s0 + s1 + s12
        s10 : s0 + s1 + s14
        s11 : s0 + s1 + s6
        s12 : s0 + s1 + s9
        s13 : s0 + s1 + s7
        s14 : s0 + s1 + s10
        s15 : s0 + s1 + s8
        sage: s4*s5
        s1 + s2 + s3 + s6 + s7 + s8 + s9 + s10 + s11 + s12 + s13 + s14 + s15
        sage: s4.ribbon()
        1
        sage: s5.ribbon()
        -1
        sage: s8.ribbon()
        zeta5^3

    If the Fusion Double is multiplicity-free, meaning that the fusion coefficients
    `N_k^{ij}` are bounded by `1`, then the F-matrix may be computed, just as
    for :class:`FusionRing`. There is a caveat here, since even if the fusion
    rules are multiplicity-free, if there are many simple objects, the computation
    may be too large to be practical. At least, this code can compute the F-matrix
    for the Fusion Double of the symmetric group `S_3`, duplicating the result
    of [CHW2015]_ .

    EXAMPLES::

        sage: G1 = SymmetricGroup(3)
        sage: H1 = FusionRing(G1, prefix="u", inject_variables=True)
        sage: F = H1.get_fmatrix()

    The above commands create the F-matrix factory. You can compute the F-matrices
    with the command ::

        sage: H1.find_orthogonal_solution()  # not tested (10-15 minutes)

    Individual F-matrices may be computed thus ::

        sage: F.fmatrix(u3,u3,u3,u4) # not tested

    see :class:`FMatrix` for more information.

    Unfortunately beyond `S_3` the number of simple objects is larger.
    Although the :class:`FusionRing` class and its methods work well
    for groups of moderate size, the FMatrix may not be available.
    For the dihedral group of order 8, there are already 22
    simple objects, and the F-matrix seems out of reach.

    It is an open problem to classify the finite groups whose fusion doubles are
    multiplicity-free. Abelian groups, dihedral groups, dicyclic groups, and all
    groups of order 16 are multiplicity-free.  On the other hand, for groups of order 32,
    some are multiplicity-free and others are not.

    Some groups, as currently implemented in Sage, may be missing methods such as
    centralizers which are needed by this code. To circumvent this, you may
    try to implement them as GAP Permutation groups. Thus ::

        sage: G1 = GL(2,3)
        sage: G2 = G1.as_permutation_group()
        sage: H2 = FusionRing(G2, prefix="b", inject_variables=True)
        sage: b13^2         # long time (43s)
        b0 + b1 + b5 + b6 + b13 + b26 + b30 + b31 + b32 + b33 + b38 + b39
        sage: b13.ribbon()
        zeta3

    In this example, implementing the simple group of order 168 as
    the matrix group ``G1`` will not work with the ``FusionRing``, so we
    recreate it as the permutation group ``G2``. Although the test of
    squaring `b2` takes a long time, the fusion coefficients are cached
    and this FusionRing is not too slow to work with. (Of course the
    F-matrix factory is not available for this group.)
    """
    def __init__(self, G, base_ring=ZZ, prefix="s", cyclotomic_order=None, fusion_labels=None, inject_variables=False):
        """
        EXAMPLES::

            sage: H = FusionRing(DihedralGroup(7))
        """
        self._G = G
        self._elt = {}
        self._chi = {}
        count = 0
        names = {}
        for g in G.conjugacy_classes_representatives():
            for chi in G.centralizer(g).irreducible_characters():
                names[count] = "%s%s"%(prefix, count)
                self._elt[count] = g
                self._chi[count] = chi
                count += 1
        cyclotomic_order = G.exponent()
        super().__init__(names=names, base_ring=ZZ, prefix=prefix, cyclotomic_order=cyclotomic_order, fusion_labels=fusion_labels, inject_variables=inject_variables)

    def _repr_(self):
        """
        EXAMPLES::

            sage: FusionRing(SymmetricGroup(3))
            The Fusion Ring of the Drinfeld Double of Symmetric group of order 3! as a permutation group
        """
        return "The Fusion Ring of the Drinfeld Double of %s"%self._G

    def group(self):
        """
        Returns the underlying group.

        EXAMPLES::

            sage: FusionRing(DiCyclicGroup(4)).group()
            Diyclic group of order 16 as a permutation group
        """
        return self._G

    def IdGroup(self):
        """
        Returns the GAP Small Group identifier. This is a pair ``[n,k]`` where ``n`` is
        the order of the group, and ``k`` is an integer characterizing the
        isomorphism class of the group, available for very many groups.

        EXAMPLES::

            sage: FusionRing(DiCyclicGroup(4)).IdGroup()
            [ 16, 9 ]
        """
        return self._G._libgap_().IdGroup()

    @cached_method
    def s_ij(self, i, j, unitary=False, base_coercion=True):
        r"""
        Return the element of the S-matrix of this fusion ring
        corresponding to the given elements. Without the unitary option
        set true, this is the unnormalized S-matrix entry, denoted `\tilde{s}_{ij}`,
        in [BaKi2001]_ Chapter 3. The normalized S-matrix entries are
        denoted `s_{ij}`.

        INPUT:

        - ``i``, ``j``, -- a pair of basis elements
        - ``unitary`` (optional): set true for the unitary normalized S-matrix.

        EXAMPLES::

            sage: D = FusionRing(SymmetricGroup(3), prefix="t", inject_variables=True)
            sage: [D.s_ij(t2, x) for x in D.basis()]
            [2, 2, 4, 0, 0, -2, -2, -2]
            sage: [D.s_ij(t2, x, unitary=True) for x in D.basis()]
            [1/3, 1/3, 2/3, 0, 0, -1/3, -1/3, -1/3]
        """
        ret = 0
        G = self._G
        [i, j] = [x.support_of_term() for x in [i,j]]
        [a, chi_1] = [self._elt[i], self._chi[i]]
        [b, chi_2] = [self._elt[j], self._chi[j]]
        for g in G:
            if a*g*b*g.inverse() == g*b*g.inverse()*a:
                ret += chi_1(g*b*g.inverse()) * chi_2(g.inverse()*a*g)
        ret *= G.order() / (G.centralizer(a).order() * G.centralizer(b).order())
        if unitary:
            ret /= self.total_q_order(base_coercion=False)
        if (not base_coercion) or (self._basecoer is None):
            return ret
        return self._basecoer(ret)

    @cached_method
    def Nk_ij(self, i, j, k):
        r"""
        Returns the fusion coefficient `N^k_{ij}`, computed using the Verlinde formula:

        .. MATH::
            N^k_{ij} = \sum_l \frac{s(i, \ell)\, s(j, \ell)\, \overline{s(k, \ell)}}{s(I, \ell)},

        EXAMPLES::

            sage: A = FusionRing(AlternatingGroup(4), prefix="aa", inject_variables=True)
            sage: [A.Nk_ij(aa8, aa10, x) for x in A.basis()]
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 1]
        """
        sz = self.one()
        return ZZ(sum(self.s_ij(i, r, unitary=True) * self.s_ij(j, r, unitary=True) * self.s_ijconj(k, r, unitary=True)/self.s_ij(sz, r, unitary=True) for r in self.basis()))

    def char_Nk_ij(self, i, j, k):
        r"""
        Use the character theoretic method to compute the fusion coefficient `N_{ij}^k`.
        This should be functionally equivalent to :meth:`Nk_ij`, and testing shows
        that it is, but it is slower.

        Each simple object, for example `i` corresponds to a conjugacy class `\mathcal{C}_i`
        of the underlying group `G`, and an irreducible character `\chi_i` of the
        centralizer `C(g_i)` of a fixed representative `g_i` of `\mathcal{C}_i`. In addition
        to the fixed representative `g_k` of the class `\mathcal{C}_i`
        and `\mathcal{C}_j`, the formula will make use of variable elements `h_i` and
        `h_j` that are subject to the condition `h_ih_j=g_k`.

        .. MATH::

            \frac{|\mathcal{C}_k|}{|G|}\sum_{\substack{h_i\in\mathcal{C}_i\\ h_j\in\mathcal{C}_j\\ h_ih_j=g_k}}|C(h_i)\cap C(h_j)|\,
            \langle\chi_i^{(h_i)}\chi_j^{(h_j)},\chi_k\rangle_{C(h_i)\cap C(h_j)},

        where `\chi_i^{(h_i)}` is the character `\chi_i` of `C(g_i)` conjugated to a
        character of `C(h_i)`, when `h_i` is a conjugate of the fixed representative `g_i`.
        More exactly, there exists `r_i` such that `r_i g_i r_i^{-1}=h_i`, and then
        `\chi_i^{(h_i)}(x)=\chi_i(r_i^{-1}xr_i)`, and this definition does not
        depend on the choice of `r_i`.

        This formula is due to Christopher Goff [Goff1999]_ when the centralizers are normal,
        and to Wenqi Li in the general case.

        EXAMPLES::

            sage: B = FusionRing(CyclicPermutationGroup(2))
            sage: all(B.char_Nk_ij(x,y,z)==B.Nk_ij(x,y,z) for x in B.basis() for y in B.basis() for z in B.basis())
            True
        """
        G = self._G
        I = G.conjugacy_class(i.g())
        J = G.conjugacy_class(j.g())
        IJ = Set(I_elem * J_elem for I_elem in I for J_elem in J)
        if k.g() not in IJ:
            return 0

        K = G.conjugacy_class(k.g())
        CI = G.centralizer(i.g())
        CJ = G.centralizer(j.g())
        CK = G.centralizer(k.g())

        c = K.cardinality() / G.order()
        summands = [(I_elem, J_elem) for I_elem in I for J_elem in J if I_elem * J_elem == k.g()]
        res = 0
        for p in summands:
            I_elem, J_elem = p
            for g in G:
                if g.inverse() * i.g() * g == I_elem:
                    i_twist = g
                if g.inverse() * j.g() * g == J_elem:
                    j_twist = g
            A = Set(i_twist.inverse() * zi * i_twist for zi in CI)
            B = Set(j_twist.inverse() * zj * j_twist for zj in CJ)
            inner_summands = A.intersection(B).intersection(Set(CK))
            for x in inner_summands:
                res += i.char()(i_twist * x * i_twist.inverse()) * j.char()(j_twist * x * j_twist.inverse()) * k.char()(x).conjugate()
        return c * res

    def one(self):
        """
        The unit element of the ring, which is the first basis element.

        EXAMPLES::

            sage: FusionRing(CyclicPermutationGroup(2), prefix="h").one()
            h0
        """
        return self.basis()[0]

    class Element(FusionRing.Element):
        def g(self):
            r"""
            The data determining a simple object consists of a conjugacy
            class representative `g` and an irreducible character `\chi` of
            the centralizer of `g`.

            Returns the conjugacy class representative of the underlying
            group corresponding to a simple object. See also :meth:`char`.

            EXAMPLES::

                sage: G = QuaternionGroup()
                sage: H = FusionRing(G, prefix="q", inject_variables=True)
                sage: q10.g()
                (1,3)(2,4)(5,7)(6,8)
                sage: q10.char()
                Character of Subgroup generated by [(1,2,3,4)(5,6,7,8), (1,5,3,7)(2,8,4,6)] of
                (Quaternion group of order 8 as a permutation group)
            """
            return self.parent()._elt[self.support_of_term()]

        def char(self):
            r"""
            The data determining a simple object consists of a conjugacy
            class representative `g` and an irreducible character `\chi` of
            the centralizer of `g`.

            Returns the character `chi`. See also :meth:`g`.

            EXAMPLES::

                sage: G = DihedralGroup(5)
                sage: H = FusionRing(G, prefix="d", inject_variables=True)
                sage: d10.g()
                (1,2,3,4,5)
                sage: d10.char()
                Character of Subgroup generated by [(1,2,3,4,5)] of (Dihedral group of order 10 as a permutation group)
            """
            return self.parent()._chi[self.support_of_term()]

        @cached_method
        def ribbon(self, base_coercion=True):
            """
            The twist or ribbon of the simple object.

            EXAMPLES::

                sage: H = FusionRing(CyclicPermutationGroup(3))
                sage: [i.ribbon() for i in H.basis()]
                [1, 1, 1, 1, zeta3, -zeta3 - 1, 1, -zeta3 - 1, zeta3]
            """
            ret = self.char()(self.g()) / self.char()(self.parent()._G.one())
            if (not base_coercion) or (self.parent()._basecoer is None):
                return ret
            return self.parent()._basecoer(ret)
