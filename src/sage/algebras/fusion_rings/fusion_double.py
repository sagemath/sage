"""
The Fusion Ring of the Drinfeld Double of a Finite Group
"""

# ****************************************************************************
#  Copyright (C) 2023 Wenqi Li
#                     Daniel Bump <bump at match.stanford.edu>
#                     Travis Scrimshaw <tcscrims at gmail.com>
#                     Guillermo Aboumrad <gh_willieab>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.categories.algebras_with_basis import AlgebrasWithBasis
from sage.combinat.free_module import CombinatorialFreeModule
from sage.rings.integer_ring import ZZ
from sage.misc.misc import inject_variable
from sage.misc.cachefunc import cached_method
from sage.sets.set import Set
from sage.rings.number_field.number_field import CyclotomicField
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.ideal import Ideal
from sage.matrix.constructor import matrix


class FusionDouble(CombinatorialFreeModule):
    r"""
    The fusion ring corresponding to the Drinfeld double of a finite group.

    This is the fusion ring of the modular tensor category of modules
    over the Drinfeld double of a finite group. Usage is similar
    to :class:`FusionRing`; we refer the reader to that class for more
    information.

    INPUT:

    - ``G`` -- a finite group
    - ``prefix`` -- (default: ``'s'``) a prefix for the names of simple objects
    - ``inject_varables`` -- (optional) set to ``True`` to create variables
      for the simple objects

    REFERENCES:

    - [BaKi2001]_ Chapter 3
    - [Mas1995]_
    - [CHW2015]_
    - [Goff1999]_

    EXAMPLES::

        sage: G = DihedralGroup(5)
        sage: H = FusionDouble(G, inject_variables=True)
        sage: H.basis()
        Finite family {0: s0, 1: s1, 2: s2, 3: s3, 4: s4, 5: s5, 6: s6, 7: s7, 8: s8,
                       9: s9, 10: s10, 11: s11, 12: s12, 13: s13, 14: s14, 15: s15}
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

    If the fusion double is multiplicity-free, meaning that the fusion
    coefficients `N_k^{ij}` are bounded by `1`, then the F-matrix may be
    computed, by solving the pentagon and hexagon relations as described
    in [Bond2007]_ and [Ab2022]_, just as for :class:`FusionRing`.
    There is a caveat here, since even if the fusion rules are multiplicity-free,
    if there are too many F-matrix values to compute, even if many of them are
    zero, in the current implementation singular cannot create enough variables.
    At least, this code can compute the F-matrix for the Fusion Double of the
    symmetric group `S_3`, duplicating the result of [CHW2015]_.

    ::

        sage: G1 = SymmetricGroup(3)
        sage: H1 = FusionDouble(G1, prefix='u', inject_variables=True)
        sage: F = H1.get_fmatrix()

    The above commands create the F-matrix. You can compute all of the
    F-matrices with the command::

        sage: H1.find_orthogonal_solution()  # not tested (10-15 minutes)

    Individual F-matrices may be computed thus::

        sage: F.fmatrix(u3, u3, u3, u4)  # not tested

    See :class:`FMatrix` for more information.

    Unfortunately beyond `S_3` the number of simple objects is seemingly
    impractical. Although the :class:`FusionDouble` class and its methods
    work well for groups of moderate size, the :class:`FMatrix` may not be
    computable. For the dihedral group of order 8, there are already 22
    simple objects, and the F-matrix seems out of reach. The actual limitation
    is that singular will not create a polynomial ring in more than
    `2^{15}-1 = 32767` symbols, and there are more than this many F-matrix
    values to be computed for the dihedral group of order 8, so in the
    current implementation, this FusionRing is out of reach.

    It is an open problem to classify the finite groups whose fusion doubles
    are multiplicity-free. Abelian groups, dihedral groups, dicyclic groups,
    and all groups of order 16 are multiplicity-free.  On the other hand, for
    groups of order 32, some are multiplicity-free and others are not.
    These can all be constructed using :class:`SmallPermutationGroup`.

    EXAMPLES::

        sage: G = SmallPermutationGroup(16,9)
        sage: F = FusionDouble(G, prefix='b', inject_variables=True)
        sage: b13^2 # long time (4s)
        b0 + b2 + b4 + b15 + b16 + b17 + b18 + b24 + b26 + b27
    """
    @staticmethod
    def __classcall_private__(cls, G, prefix='s', inject_variables=False):
        """
        Normalize input to ensure a unique representation.

        EXAMPLES::

            sage: H1 = FusionDouble(DihedralGroup(6), inject_variables=True)
            sage: H2 = FusionDouble(DihedralGroup(6), prefix='s')
            sage: H1 is H2
            True
        """
        F = super().__classcall__(cls, G=G, prefix=prefix)
        if inject_variables:
            F.inject_variables()
        return F

    def __init__(self, G, prefix='s'):
        """
        EXAMPLES::

            sage: H = FusionDouble(DihedralGroup(6))
            sage: TestSuite(H).run()
            sage: H = FusionDouble(DihedralGroup(7))
            sage: TestSuite(H).run()  # long time

            sage: F = FusionDouble(CyclicPermutationGroup(2))
            sage: [F._repr_term(t) for t in F._names]
            ['s0', 's1', 's2', 's3']
            sage: F = FusionDouble(CyclicPermutationGroup(2))
            sage: [F._latex_term(t) for t in F._names]
            ['s_{0}', 's_{1}', 's_{2}', 's_{3}']

            sage: FusionDouble(SymmetricGroup(4)).get_order()
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
        """
        self._G = G
        self._names = {}
        self._elt = {}
        self._chi = {}
        self._unit_index = None  # index of the unit element
        count = ZZ.zero()
        for g in G.conjugacy_classes_representatives():
            for chi in G.centralizer(g).irreducible_characters():
                # NOTE: the trivial char is not necessarily the first one
                self._names[count] = "%s%s" % (prefix, count)
                self._elt[count] = g
                self._chi[count] = chi
                if self._unit_index is None and all(v == 1 for v in chi):
                    self._unit_index = count
                count += ZZ.one()
        self._cyclotomic_order = G.exponent()
        self._basecoer = None
        self._fusion_labels = None
        self._field = None
        cat = AlgebrasWithBasis(ZZ)
        CombinatorialFreeModule.__init__(self, ZZ, list(self._names),
                                         prefix=prefix, bracket=False, category=cat)

    def _repr_(self):
        """
        EXAMPLES::

            sage: FusionDouble(SymmetricGroup(3))
            The Fusion Ring of the Drinfeld Double of Symmetric group of
                order 3! as a permutation group
        """
        return "The Fusion Ring of the Drinfeld Double of %s" % self._G

    def inject_variables(self):
        """
        Create variables for the simple objects in the global name space.

        EXAMPLES::

            sage: F = FusionDouble(DiCyclicGroup(3), prefix='d')
            sage: F.inject_variables()
            sage: d0 + d1 + d5
            d0 + d1 + d5
        """
        for i, name in self._names.items():
            inject_variable(name, self.monomial(i))

    @cached_method
    def _char_cache(self, i, g):
        """
        Return ``self._chi[i](g)``, cached here for speed.

        TESTS::

            sage: D = FusionDouble(SymmetricGroup(4))

            sage: all(D._char_cache(b.support_of_term(), b.g()) == b.char()(b.g())
            ....:     for b in D.basis())
            True
        """
        return self._chi[i](g)

    @cached_method
    def s_ij(self, i, j, unitary=False, base_coercion=True):
        r"""
        Return the element of the `S`-matrix of this fusion ring
        corresponding to the given elements.

        Without the unitary option set true, this is the unnormalized `S`-matrix
        entry, denoted `\tilde{s}_{ij}`, in [BaKi2001]_ Chapter 3. The
        normalized `S`-matrix entries are denoted `s_{ij}`.

        INPUT:

        - ``i``, ``j``, -- a pair of basis elements
        - ``unitary`` -- boolean (default: ``False``); set to ``True`` to
          obtain the unitary `S`-matrix

        EXAMPLES::

            sage: D = FusionDouble(SymmetricGroup(3), prefix='t', inject_variables=True)
            sage: [D.s_ij(t2, x) for x in D.basis()]
            [2, 2, 4, 0, 0, -2, -2, -2]
            sage: [D.s_ij(t2, x, unitary=True) for x in D.basis()]
            [1/3, 1/3, 2/3, 0, 0, -1/3, -1/3, -1/3]
        """
        sum_val = ZZ.zero()
        G = self._G
        [i] = list(i._monomial_coefficients)
        [j] = list(j._monomial_coefficients)
        a = self._elt[i]
        b = self._elt[j]
        for g in G:
            gi = g.inverse()
            conj = g * b * gi
            if a * conj == conj * a:
                sum_val += self._char_cache(i, conj) * self._char_cache(j, gi * a * g)
        if unitary:
            coef = 1 / (G.centralizer(a).order() * G.centralizer(b).order())
        else:
            coef = G.order() / (G.centralizer(a).order() * G.centralizer(b).order())
        ret = coef * sum_val
        if (not base_coercion) or (self._basecoer is None):
            return ret
        return self._basecoer(ret)

    def s_ijconj(self, i, j, unitary=False, base_coercion=True):
        r"""
        Return the conjugate of the element of the `S`-matrix given by
        ``self.s_ij(elt_i, elt_j, base_coercion=base_coercion)``.

        .. SEEALSO::

            :meth:`s_ij`

        EXAMPLES::

            sage: P=FusionDouble(CyclicPermutationGroup(3),prefix='p',inject_variables=True)
            sage: P.s_ij(p1,p3)
            zeta3
            sage: P.s_ijconj(p1,p3)
            -zeta3 - 1
        """
        return self.s_ij(i, j, unitary=unitary, base_coercion=base_coercion).conjugate()

    def s_matrix(self, unitary=False, base_coercion=True):
        r"""
        Return the `S`-matrix of this fusion ring.

        OPTIONAL:

        - ``unitary`` -- boolean (default: ``False``); set to ``True`` to
          obtain the unitary `S`-matrix

        Without the ``unitary`` parameter, this is the matrix denoted
        `\widetilde{s}` in [BaKi2001]_.

        EXAMPLES::

            sage: FusionDouble(SymmetricGroup(3)).s_matrix()
            [ 1  1  2  3  3  2  2  2]
            [ 1  1  2 -3 -3  2  2  2]
            [ 2  2  4  0  0 -2 -2 -2]
            [ 3 -3  0  3 -3  0  0  0]
            [ 3 -3  0 -3  3  0  0  0]
            [ 2  2 -2  0  0  4 -2 -2]
            [ 2  2 -2  0  0 -2 -2  4]
            [ 2  2 -2  0  0 -2  4 -2]
            sage: FusionDouble(SymmetricGroup(3)).s_matrix(unitary=True)
            [ 1/6  1/6  1/3  1/2  1/2  1/3  1/3  1/3]
            [ 1/6  1/6  1/3 -1/2 -1/2  1/3  1/3  1/3]
            [ 1/3  1/3  2/3    0    0 -1/3 -1/3 -1/3]
            [ 1/2 -1/2    0  1/2 -1/2    0    0    0]
            [ 1/2 -1/2    0 -1/2  1/2    0    0    0]
            [ 1/3  1/3 -1/3    0    0  2/3 -1/3 -1/3]
            [ 1/3  1/3 -1/3    0    0 -1/3 -1/3  2/3]
            [ 1/3  1/3 -1/3    0    0 -1/3  2/3 -1/3]
        """
        b = self.basis()
        S = matrix([[self.s_ij(b[x], b[y], unitary=unitary, base_coercion=base_coercion)
                     for x in self.get_order()] for y in self.get_order()])
        return S

    @cached_method
    def N_ijk(self, i, j, k):
        r"""
        The symmetric invariant of three simple objects.

        This is the dimension of

        .. MATH::

           Hom(i \otimes j \otimes k, s_0),

        where `s_0` is the unit element (assuming ``prefix='s'``).
        Method of computation is through the Verlinde formula,
        deducing the values from the known values of the `S`-matrix.

        EXAMPLES::

            sage: A = FusionDouble(AlternatingGroup(4),prefix='a',inject_variables=True)
            sage: [A.N_ijk(a10,a11,x) for x in A.basis()]
            [0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0]

        TESTS::

            sage: F = FusionDouble(SymmetricGroup(4))
            sage: from itertools import product
            sage: B = list(F.basis())
            sage: all(F.N_ijk(i,j,k).parent() is ZZ
            ....:     for i, j, k in product(B[::6], repeat=3))
            True
        """
        sz = self.one()
        return ZZ(sum(self.s_ij(i, r, unitary=True) * self.s_ij(j, r, unitary=True)
                      * self.s_ij(k, r, unitary=True) / self.s_ij(sz, r, unitary=True)
                      for r in self.basis()))

    @cached_method
    def Nk_ij(self, i, j, k, use_characters=False):
        r"""
        Return the fusion coefficient `N^k_{ij}`.

        INPUT:

        - ``i``, ``j``, ``k`` -- basis elements
        - ``use_characters`` -- boolean (default: ``False``); see the algorithm
          description below

        ALGORITHM:

        If ``use_characters=False``, then this is computed using
        the Verlinde formula:

        .. MATH::

            N^k_{ij} = \sum_l \frac{s(i, \ell)\, s(j, \ell)\,
            \overline{s(k, \ell)}}{s(I, \ell)}.

        Otherwise we use a character theoretic method to compute the fusion
        coefficient `N_{ij}^k` as follows. Each simple object, for example
        `i` corresponds to a conjugacy class `\mathcal{C}_i` of the underlying
        group `G`, and an irreducible character `\chi_i` of the centralizer
        `C(g_i)` of a fixed representative `g_i` of `\mathcal{C}_i`. In addition
        to the fixed representative `g_k` of the class `\mathcal{C}_i`
        and `\mathcal{C}_j`, the formula will make use of variable elements
        `h_i` and `h_j` that are subject to the condition `h_i h_j = g_k`.
        See [GoMa2010]_ equation (7).

        .. MATH::

            \frac{|\mathcal{C}_k|}{|G|}
            \sum_{\substack{h_i\in\mathcal{C}_i \\ h_j\in\mathcal{C}_j \\ h_ih_j=g_k}}
            \lvert C(h_i)\cap C(h_j) \rvert \,
            \langle \chi_i^{(h_i)} \chi_j^{(h_j)}, \chi_k \rangle_{C(h_i)\cap C(h_j)},

        where `\chi_i^{(h_i)}` is the character `\chi_i` of `C(g_i)`
        conjugated to a character of `C(h_i)`, when `h_i` is a conjugate
        of the fixed representative `g_i`. More exactly, there exists `r_i`
        such that `r_i g_i r_i^{-1} = h_i`, and then `\chi_i^{(h_i)}(x) =
        \chi_i(r_i^{-1}xr_i)`, and this definition does not depend on the
        choice of `r_i`.

        .. NOTE::

            This should be functionally equivalent, and testing shows
            that it is, but it is slower.

        EXAMPLES::

            sage: A = FusionDouble(AlternatingGroup(4),prefix='aa',inject_variables=True)
            sage: [A.Nk_ij(aa8,aa10,x) for x in A.basis()]
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 1]

            sage: B = FusionDouble(CyclicPermutationGroup(2))
            sage: all(B.Nk_ij(x,y,z,use_characters=True) == B.Nk_ij(x,y,z)
            ....:     for x in B.basis() for y in B.basis() for z in B.basis())
            True
        """
        if not use_characters:
            return self.N_ijk(i, j, self.dual(k))

        G = self._G
        I = G.conjugacy_class(i.g())
        J = G.conjugacy_class(j.g())
        IJ = {I_elem * J_elem for I_elem in I for J_elem in J}
        if k.g() not in IJ:
            return ZZ.zero()

        K = G.conjugacy_class(k.g())
        CI = G.centralizer(i.g())
        CJ = G.centralizer(j.g())
        CK = G.centralizer(k.g())

        c = K.cardinality() / G.order()
        summands = [(I_elem, J_elem) for I_elem in I for J_elem in J if I_elem * J_elem == k.g()]
        res = ZZ.zero()
        ichar = i.char()
        jchar = j.char()
        kchar = k.char()
        for p in summands:
            I_elem, J_elem = p
            for g in G:
                ginv = g.inverse()
                if ginv * i.g() * g == I_elem:
                    i_twist = g
                if ginv * j.g() * g == J_elem:
                    j_twist = g
            A = Set(i_twist.inverse() * zi * i_twist for zi in CI)
            B = Set(j_twist.inverse() * zj * j_twist for zj in CJ)
            inner_summands = A.intersection(B).intersection(Set(CK))
            i_twist_inv = i_twist.inverse()
            j_twist_inv = j_twist.inverse()
            res += sum(ichar(i_twist * x * i_twist_inv) * jchar(j_twist * x * j_twist_inv) * kchar(x).conjugate()
                       for x in inner_summands)
        return c * res

    @cached_method
    def field(self):
        """
        Return a cyclotomic field large enough to contain the values
        of R-matrices and twists that can arise for this fusion ring.

        EXAMPLES::

            sage: FusionDouble(SymmetricGroup(3)).field()
            Cyclotomic Field of order 24 and degree 8
        """
        return CyclotomicField(4 * self._cyclotomic_order)

    def fvars_field(self):
        r"""
        Return a field containing the ``CyclotomicField`` computed by
        :meth:`field` as well as all the F-symbols of the associated
        ``FMatrix`` factory object.

        This method is only available if ``self`` is multiplicity-free.

        EXAMPLES::

            sage: FusionDouble(SymmetricGroup(3)).fvars_field()
            Cyclotomic Field of order 24 and degree 8
        """
        if self.is_multiplicity_free(verbose=False):
            return self.get_fmatrix().field()
        raise NotImplementedError("method is only available for multiplicity free fusion rings")

    def root_of_unity(self, r, base_coercion=True):
        r"""
        Return `e^{i\pi r}` as an element of ``self.field()`` if possible.

        INPUT:

        - ``r`` -- a rational number

        EXAMPLES::

            sage: H = FusionDouble(DihedralGroup(6))
            sage: H.field()
            Cyclotomic Field of order 24 and degree 8
            sage: for n in [1..7]:
            ....:     try:
            ....:         print (n,H.root_of_unity(2/n))
            ....:     except ValueError as err:
            ....:         print (n,err)
            ....:
            1 1
            2 -1
            3 zeta24^4 - 1
            4 zeta24^6
            5 not a root of unity in the field
            6 zeta24^4
            7 not a root of unity in the field
        """
        n = 2 * r * self._cyclotomic_order
        if n not in ZZ:
            raise ValueError("not a root of unity in the field")
        ret = self.field().gen() ** n
        if (not base_coercion) or (self._basecoer is None):
            return ret
        return self._basecoer(ret)

    @cached_method
    def r_matrix(self, i, j, k, base_coercion=True):
        r"""
        Return the R-matrix entry corresponding to the subobject ``k``
        in the tensor product of ``i`` with ``j``. This method is only
        correct if the fusion coefficient ``N_{ij}^k\leq 1``. See the
        :class:`FusionRing` method for more information, including
        the reason for this caveat, and the algorithm.

        EXAMPLES::

            sage: C = FusionDouble(SymmetricGroup(3),prefix='c',inject_variables=True)
            sage: c4*c5
            c3 + c4
            sage: [C.r_matrix(c4,c5,k) for k in [c3,c4]]
            [-zeta24^6, 1]
            sage: c6^2
            c0 + c1 + c6
            sage: [C.r_matrix(c6,c6,k) for k in [c0,c1,c6]]
            [zeta3, -zeta3, -zeta3 - 1]
        """
        if self.Nk_ij(i, j, k) == 0:
            return self.field().zero() if (not base_coercion) or (self._basecoer is None) else self.fvars_field().zero()

        if i != j:
            ret = self.root_of_unity((k.twist() - i.twist() - j.twist()) / 2)
        else:
            i0 = self.one()
            B = self.basis()
            ret = sum(y.ribbon()**2 / (i.ribbon() * x.ribbon()**2)
                      * self.s_ij(i0, y) * self.s_ij(i, z) * self.s_ijconj(x, z)
                      * self.s_ijconj(k, x) * self.s_ijconj(y, z) / self.s_ij(i0, z)
                      for x in B for y in B for z in B) / (self.total_q_order()**4)

        if (not base_coercion) or (self._basecoer is None):
            return ret
        return self._basecoer(ret)

    def global_q_dimension(self, base_coercion=True):
        r"""
        Return the global quantum dimension, which is the sum of the squares of the
        quantum dimensions of the simple objects.
        For the Drinfeld double, it is the square of the order of the underlying quantum group.

        EXAMPLES::

            sage: G = SymmetricGroup(4)
            sage: H = FusionDouble(G)
            sage: H.global_q_dimension()
            576
            sage: sum(x.q_dimension()^2 for x in H.basis())
            576
        """
        ret = self._G.order() ** 2
        if (not base_coercion) or (self._basecoer is None):
            return ret
        return self._basecoer(ret)

    def total_q_order(self, base_coercion=True):
        r"""
        Return the positive square root of :meth:`self.global_q_dimension()
        <global_q_dimension>` as an element of :meth:`self.field() <field>`.

        For the Drinfeld double of a finite group `G`, this equals the
        cardinality of `G`. This is also equal to `\sum d_i^2 \theta_i^{\pm 1}`,
        where `i` runs through the simple objects, `d_i` is the quantum
        dimension, and `\theta_i` is the twist. This sum with `\theta_i` is
        denoted `p_-` in [BaKi2001]_ Chapter 3.

        EXAMPLES::

            sage: FusionDouble(DihedralGroup(7)).total_q_order()
            14
        """
        ret = self._G.order()
        if (not base_coercion) or (self._basecoer is None):
            return ret
        return self._basecoer(ret)

    D_minus = D_plus = total_q_order

    def is_multiplicity_free(self, verbose=False):
        """
        Return ``True`` if all fusion coefficients are at most 1.

        EXAMPLES::

            sage: FusionDouble(SymmetricGroup(3)).is_multiplicity_free()
            True
            sage: FusionDouble(SymmetricGroup(4)).is_multiplicity_free()
            False

            sage: FusionDouble(SymmetricGroup(3)).is_multiplicity_free(True)
            Checking multiplicity freeness
            True
            sage: FusionDouble(SymmetricGroup(4)).is_multiplicity_free(True)
            Checking multiplicity freeness
            N(s2,s13,s13) = 2
            False
        """
        if verbose:
            print("Checking multiplicity freeness")
            from itertools import product
            for (i, j, k) in product(self.basis(), repeat=3):
                if self.N_ijk(i, j, k) > 1:
                    print("N(%s,%s,%s) = %s" % (i, j, k, self.N_ijk(i, j, k)))
                    return False
            return True

        return all(self.N_ijk(i,j,k) <= 1 for i in self.basis()
                   for j in self.basis() for k in self.basis())

    @cached_method
    def one_basis(self):
        r"""
        The unit element of the ring, which is the first basis element.

        EXAMPLES::

            sage: FusionDouble(CyclicPermutationGroup(2), prefix='h').one()
            h1
        """
        return self._unit_index

    @cached_method
    def dual(self, i):
        r"""
        Return the dual object ``i^\ast`` to ``i``.

        The dual is also available as an element method of ``i``.

        EXAMPLES::

            sage: K = FusionDouble(CyclicPermutationGroup(3),prefix='k')
            sage: [(x,K.dual(x)) for x in K.basis()]
            [(k0, k0),
            (k1, k2),
            (k2, k1),
            (k3, k6),
            (k4, k8),
            (k5, k7),
            (k6, k3),
            (k7, k5),
            (k8, k4)]
            sage: all(K.dual(x)==x.dual() for x in K.basis())
            True
        """
        sz = self.one()
        for j in self.basis():
            if self.N_ijk(i, j, sz) > 0:
                return j

    def product_on_basis(self, a, b):
        """
        Return the product of two basis elements corresponding to keys `a` and `b`.

        INPUT:

        - ``a``, ``b`` -- keys for the dictionary ``self._names`` representing
          simple objects

        EXAMPLES::

            sage: Q=FusionDouble(SymmetricGroup(3),prefix='q',inject_variables=True)
            sage: q3*q4
            q1 + q2 + q5 + q6 + q7
            sage: Q._names
            {0: 'q0', 1: 'q1', 2: 'q2', 3: 'q3', 4: 'q4', 5: 'q5', 6: 'q6', 7: 'q7'}
            sage: Q.product_on_basis(3,4)
            q1 + q2 + q5 + q6 + q7
        """
        d = {k.support_of_term(): val for k in self.basis()
             if (val := self.N_ijk(self.monomial(a),self.monomial(b),self.dual(k)))}
        return self._from_dict(d, remove_zeros=False)

    def group(self):
        """
        Return the underlying group.

        EXAMPLES::

            sage: FusionDouble(DiCyclicGroup(4)).group()
            Dicyclic group of order 16 as a permutation group
        """
        return self._G

    def get_fmatrix(self, *args, **kwargs):
        r"""
        Construct an :class:`FMatrix` factory to solve the pentagon and
        hexagon relations and organize the resulting F-symbols.

        EXAMPLES::

            sage: f = FusionDouble(SymmetricGroup(3)).get_fmatrix(); f
            F-Matrix factory for The Fusion Ring of the Drinfeld Double of
                Symmetric group of order 3! as a permutation group
        """
        if not hasattr(self, 'fmats') or kwargs.get('new', False):
            kwargs.pop('new', None)
            from sage.algebras.fusion_rings.f_matrix import FMatrix
            self.fmats = FMatrix(self, *args, **kwargs)
        return self.fmats

    class Element(CombinatorialFreeModule.Element):
        def is_simple_object(self):
            r"""
            Determine whether ``self`` is a simple object (basis element) of the fusion ring.

            EXAMPLES::

                sage: H = FusionDouble(CyclicPermutationGroup(2), prefix='g', inject_variables=True)
                sage: [x.is_simple_object() for x in [g0, g1, g0+g1]]
                [True, True, False]
            """
            return len(self._monomial_coefficients) == 1

        def g(self):
            r"""
            The data determining a simple object consists of a conjugacy
            class representative `g` and an irreducible character `\chi` of
            the centralizer of `g`.

            Returns the conjugacy class representative of the underlying
            group corresponding to a simple object. See also :meth:`char`.

            EXAMPLES::

                sage: G = QuaternionGroup()
                sage: H = FusionDouble(G, prefix='e', inject_variables=True)
                sage: e10.g()
                (1,3)(2,4)(5,7)(6,8)
                sage: e10.char()
                Character of Subgroup generated by [(1,2,3,4)(5,6,7,8), (1,5,3,7)(2,8,4,6)]
                    of (Quaternion group of order 8 as a permutation group)
            """
            return self.parent()._elt[self.support_of_term()]

        def char(self):
            r"""
            Return the character `\chi` corresponding to ``self``.

            The data determining a simple object consists of a conjugacy
            class representative `g` and an irreducible character `\chi` of
            the centralizer of `g`.

            .. SEEALSO:: :meth:`g`

            EXAMPLES::

                sage: G = DihedralGroup(5)
                sage: H = FusionDouble(G, prefix='f', inject_variables=True)
                sage: f10.g()
                (1,2,3,4,5)
                sage: f10.char()
                Character of Subgroup generated by [(1,2,3,4,5)] of
                    (Dihedral group of order 10 as a permutation group)
            """
            return self.parent()._chi[self.support_of_term()]

        def ribbon(self, base_coercion=True):
            """
            The twist or ribbon of the simple object.

            EXAMPLES::

                sage: H = FusionDouble(CyclicPermutationGroup(3))
                sage: [i.ribbon() for i in H.basis()]
                [1, 1, 1, 1, zeta3, -zeta3 - 1, 1, -zeta3 - 1, zeta3]
            """
            i = self.support_of_term()
            P = self.parent()
            ret = P._char_cache(i, self.g()) / P._char_cache(i, P._G.one())
            if (not base_coercion) or (P._basecoer is None):
                return ret
            return self.parent()._basecoer(ret)

        def twist(self, reduced=True):
            r"""
            Return a rational number `h` such that `\theta = e^{i \pi h}`
            is the twist of ``self``.

            The quantity `e^{i \pi h}` is also available using :meth:`ribbon`.

            This method is only available for simple objects.

            EXAMPLES::

                sage: Q=FusionDouble(CyclicPermutationGroup(3))
                sage: [x.twist() for x in Q.basis()]
                [0, 0, 0, 0, 2/3, 4/3, 0, 4/3, 2/3]
                sage: [x.ribbon() for x in Q.basis()]
                [1, 1, 1, 1, zeta3, -zeta3 - 1, 1, -zeta3 - 1, zeta3]

            TESTS::

                sage: H = FusionDouble(AlternatingGroup(4))
                sage: sum(H.basis()).twist()
                Traceback (most recent call last):
                ...
                ValueError: quantum twist is only available for simple objects
            """
            if not self.is_simple_object():
                raise ValueError("quantum twist is only available for simple objects")
            P = self.parent()
            zeta = P.field().gen()
            rib = self.ribbon()
            norm = 2 * P._cyclotomic_order
            for k in range(4*P._cyclotomic_order):
                if zeta ** k == rib:
                    return k / norm

        def dual(self):
            """
            Return the dual of ``self``.

            This method is only available for simple objects.

            EXAMPLES::

                sage: G = CyclicPermutationGroup(4)
                sage: H = FusionDouble(G, prefix='j')
                sage: [x for x in H.basis() if x == x.dual()]
                [j0, j1, j8, j9]

            TESTS::

                sage: H = FusionDouble(AlternatingGroup(4))
                sage: sum(H.basis()).dual()
                Traceback (most recent call last):
                ...
                ValueError: dual is only available for simple objects
            """
            if not self.is_simple_object():
                raise ValueError("dual is only available for simple objects")
            return self.parent().dual(self)

        @cached_method
        def q_dimension(self, base_coercion=True):
            """
            Return the q-dimension of ``self``.

            This method is only available for simple objects.

            EXAMPLES::

                sage: G = AlternatingGroup(4)
                sage: H = FusionDouble(G)
                sage: [x.q_dimension() for x in H.basis()]
                [1, 1, 1, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4]
                sage: sum(x.q_dimension()^2 for x in H.basis()) == G.order()^2
                True

            TESTS::

                sage: H = FusionDouble(AlternatingGroup(4))
                sage: sum(H.basis()).q_dimension()
                Traceback (most recent call last):
                ...
                ValueError: quantum dimension is only available for simple objects
            """
            if not self.is_simple_object():
                raise ValueError("quantum dimension is only available for simple objects")
            return self.parent().s_ij(self,self.parent().one())
