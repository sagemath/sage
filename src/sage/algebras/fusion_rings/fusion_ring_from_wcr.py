"""
Fusion Rings defined by Weyl Character Rings
"""
# ****************************************************************************
#  Copyright (C) 2019 Daniel Bump <bump at match.stanford.edu>
#                     Guillermo Aboumrad <gh_willieab>
#                     Travis Scrimshaw <tcscrims at gmail.com>
#                     Nicolas Thiery <nthiery at users.sf.net>
#                2022 Guillermo Aboumrad <gh_willieab>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from sage.algebras.fusion_rings.generic_fusion_ring import FusionRing
from sage.combinat.q_analogues import q_int
from sage.combinat.root_system.weyl_characters import WeylCharacterRing
from sage.misc.cachefunc import cached_method
from sage.rings.integer_ring import ZZ

class FusionRingFromWCR(FusionRing):
    r"""
    Return the Fusion Ring (Verlinde Algebra) of level ``k``.

    INPUT:

    - ``ct`` -- the Cartan type of a simple (finite-dimensional) Lie algebra
    - ``k`` -- a nonnegative integer
    - ``conjugate`` -- (default ``False``) set ``True`` to obtain
      the complex conjugate ring
    - ``cyclotomic_order`` -- (default computed depending on ``ct`` and ``k``)
    - ``fusion_labels`` --  (default None) either a tuple of strings to use as labels of the
      basis of simple objects, or a string from which the labels will be
      constructed
    - ``inject_variables`` -- (default ``False``): use with ``fusion_labels``.
      If ``inject_variables`` is ``True``, the fusion labels will be variables
      that can be accessed from the command line

    The cyclotomic order is an integer `N` such that all computations
    will return elements of the cyclotomic field of `N`-th roots of unity.
    Normally you will never need to change this but consider changing it
    if :meth:`root_of_unity` raises a ``ValueError``.

    This algebra has a basis (sometimes called *primary fields* but here
    called *simple objects*) indexed by the weights of level `\leq k`.
    These arise as the fusion algebras of Wess-Zumino-Witten (WZW) conformal
    field theories, or as Grothendieck groups of tilting modules for quantum
    groups at roots of unity. The :class:`FusionRing` class is implemented as
    a variant of the :class:`WeylCharacterRing`.

    REFERENCES:

    - [BaKi2001]_ Chapter 3
    - [DFMS1996]_ Chapter 16
    - [EGNO2015]_ Chapter 8
    - [Feingold2004]_
    - [Fuchs1994]_
    - [Row2006]_
    - [Walton1990]_
    - [Wan2010]_

    EXAMPLES::

        sage: A22 = FusionRing(("A2", 2))
        sage: [f1, f2] = A22.fundamental_weights()
        sage: M = [A22(x) for x in [0*f1, 2*f1, 2*f2, f1+f2, f2, f1]]
        sage: [M[3] * x for x in M]
        [A22(1,1),
         A22(0,1),
         A22(1,0),
         A22(0,0) + A22(1,1),
         A22(0,1) + A22(2,0),
         A22(1,0) + A22(0,2)]

    You may assign your own labels to the basis elements. In the next
    example, we create the `SO(5)` fusion ring of level `2`, check the
    weights of the basis elements, then assign new labels to them while
    injecting them into the global namespace::

        sage: B22 = FusionRing(("B2", 2))
        sage: b = [B22(x) for x in B22.get_order()]; b
        [B22(0,0), B22(1,0), B22(0,1), B22(2,0), B22(1,1), B22(0,2)]
        sage: [x.weight() for x in b]
        [(0, 0), (1, 0), (1/2, 1/2), (2, 0), (3/2, 1/2), (1, 1)]
        sage: B22.fusion_labels(['I0', 'Y1', 'X', 'Z', 'Xp', 'Y2'], inject_variables=True)
        sage: b = [B22(x) for x in B22.get_order()]; b
        [I0, Y1, X, Z, Xp, Y2]
        sage: [(x, x.weight()) for x in b]
        [(I0, (0, 0)),
         (Y1, (1, 0)),
         (X, (1/2, 1/2)),
         (Z, (2, 0)),
         (Xp, (3/2, 1/2)),
         (Y2, (1, 1))]
        sage: X * Y1
        X + Xp
        sage: Z * Z
        I0

    A fixed order of the basis keys is available with :meth:`get_order`.
    This is the order used by methods such as :meth:`s_matrix`. You may
    use :meth:`CombinatorialFreeModule.set_order` to reorder the basis::

        sage: B22.set_order([x.weight() for x in [I0, Y1, Y2, X, Xp, Z]])
        sage: [B22(x) for x in B22.get_order()]
        [I0, Y1, Y2, X, Xp, Z]

    To reset the labels, you may run :meth:`fusion_labels` with no parameter::

        sage: B22.fusion_labels()
        sage: [B22(x) for x in B22.get_order()]
        [B22(0,0), B22(1,0), B22(0,2), B22(0,1), B22(1,1), B22(2,0)]

    To reset the order to the default, simply set it to the list of basis
    element keys::

        sage: B22.set_order(B22.basis().keys().list())
        sage: [B22(x) for x in B22.get_order()]
        [B22(0,0), B22(1,0), B22(0,1), B22(2,0), B22(1,1), B22(0,2)]

    The fusion ring has a number of methods that reflect its role
    as the Grothendieck ring of a *modular tensor category* (MTC). These
    include twist methods :meth:`Element.twist` and :meth:`Element.ribbon`
    for its elements related to the ribbon structure, and the
    S-matrix :meth:`s_ij`.

    There are two natural normalizations of the S-matrix. Both
    are explained in Chapter 3 of [BaKi2001]_. The one that is computed
    by the method :meth:`s_matrix`, or whose individual entries
    are computed by :meth:`s_ij` is denoted `\tilde{s}` in
    [BaKi2001]_. It is not unitary.

    The unitary S-matrix is `s=D^{-1/2}\tilde{s}` where

    .. MATH::

        D = \sum_V d_i(V)^2.

    The sum is over all simple objects `V` with
    `d_i(V)` the *quantum dimension*. We will call quantity `D`
    the *global quantum dimension* and `\sqrt{D}` the
    *total quantum order*. They are  computed by :meth:`global_q_dimension`
    and :meth:`total_q_order`. The unitary S-matrix `s` may be obtained
    using :meth:`s_matrix` with the option ``unitary=True``.

    Let us check the Verlinde formula, which is [DFMS1996]_ (16.3). This
    famous identity states that

    .. MATH::

        N^k_{ij} = \sum_l \frac{s(i, \ell)\, s(j, \ell)\, \overline{s(k, \ell)}}{s(I, \ell)},

    where `N^k_{ij}` are the fusion coefficients, i.e. the structure
    constants of the fusion ring, and ``I`` is the unit object.
    The S-matrix has the property that if `i*` denotes the dual
    object of `i`, implemented in Sage as ``i.dual()``, then

    .. MATH::

        s(i*, j) = s(i, j*) = \overline{s(i, j)}.

    This is equation (16.5) in [DFMS1996]_. Thus with `N_{ijk}=N^{k*}_{ij}`
    the Verlinde formula is equivalent to

    .. MATH::

        N_{ijk} = \sum_l \frac{s(i, \ell)\, s(j, \ell)\, s(k, \ell)}{s(I, \ell)},

    In this formula `s` is the normalized unitary S-matrix
    denoted `s` in [BaKi2001]_. We may define a function that
    corresponds to the right-hand side, except using
    `\tilde{s}` instead of `s`::

        sage: def V(i, j, k):
        ....:     R = i.parent()
        ....:     return sum(R.s_ij(i, l) * R.s_ij(j, l) * R.s_ij(k, l) / R.s_ij(R.one(), l)
        ....:                for l in R.basis())

    This does not produce ``self.N_ijk(i, j, k)`` exactly, because of the
    missing normalization factor. The following code to check the
    Verlinde formula takes this into account::

       sage: def test_verlinde(R):
       ....:     b0 = R.one()
       ....:     c = R.global_q_dimension()
       ....:     return all(V(i, j, k) == c * R.N_ijk(i, j, k) for i in R.basis()
       ....:                for j in R.basis() for k in R.basis())

    Every fusion ring should pass this test::

        sage: test_verlinde(FusionRing(("A2", 1)))
        True
        sage: test_verlinde(FusionRing(("B4", 2))) # long time (.56s)
        True

    As an exercise, the reader may verify the examples in
    Section 5.3 of [RoStWa2009]_. Here we check the example
    of the Ising modular tensor category, which is related
    to the BPZ minimal model `M(4, 3)` or to an `E_8` coset
    model. See [DFMS1996]_ Sections 7.4.2 and 18.4.1.
    [RoStWa2009]_ Example 5.3.4 tells us how to
    construct it as the conjugate of the `E_8` level 2
    :class:`FusionRing`::

        sage: I = FusionRing(("E8", 2), conjugate=True)
        sage: I.fusion_labels(["i0", "p", "s"], inject_variables=True)
        sage: b = I.basis().list(); b
        [i0, p, s]
        sage: Matrix([[x*y for x in b] for y in b]) # long time (.93s)
        [    i0      p      s]
        [     p     i0      s]
        [     s      s i0 + p]
        sage: [x.twist() for x in b]
        [0, 1, 1/8]
        sage: [x.ribbon() for x in b]
        [1, -1, zeta128^8]
        sage: [I.r_matrix(i, j, k) for (i, j, k) in [(s, s, i0), (p, p, i0), (p, s, s), (s, p, s), (s, s, p)]]
        [-zeta128^56, -1, -zeta128^32, -zeta128^32, zeta128^24]
        sage: I.r_matrix(s, s, i0) == I.root_of_unity(-1/8)
        True
        sage: I.global_q_dimension()
        4
        sage: I.total_q_order()
        2
        sage: [x.q_dimension()^2 for x in b]
        [1, 1, 2]
        sage: I.s_matrix()
        [                       1                        1 -zeta128^48 + zeta128^16]
        [                       1                        1  zeta128^48 - zeta128^16]
        [-zeta128^48 + zeta128^16  zeta128^48 - zeta128^16                        0]
        sage: I.s_matrix().apply_map(lambda x:x^2)
        [1 1 2]
        [1 1 2]
        [2 2 0]

    The term *modular tensor category* refers to the fact that associated
    with the category there is a projective representation of the modular
    group `SL(2, \ZZ)`. We recall that this group is generated by

    .. MATH::

        S = \begin{pmatrix} & -1\\1\end{pmatrix}, \qquad
        T = \begin{pmatrix} 1 & 1\\ &1 \end{pmatrix}

    subject to the relations `(ST)^3 = S^2`, `S^2T = TS^2`, and `S^4 = I`.
    Let `s` be the normalized S-matrix, and
    `t` the diagonal matrix whose entries are the twists of the simple
    objects. Let `s` the unitary S-matrix and `t` the matrix of twists,
    and `C` the conjugation matrix :meth:`conj_matrix`. Let

    .. MATH::

        D_+ = \sum_i d_i^2 \theta_i, \qquad D_- = d_i^2 \theta_i^{-1},

    where `d_i` and `\theta_i` are the quantum dimensions and twists of the
    simple objects. Let `c` be the Virasoro central charge, a rational number
    that is computed in :meth:`virasoro_central_charge`. It is known that

    .. MATH::

        \sqrt{\frac{D_+}{D_-}} = e^{i\pi c/4}.

    It is proved in [BaKi2001]_ Equation (3.1.17) that

    .. MATH::

        (st)^3 = e^{i\pi c/4} s^2, \qquad
        s^2 = C, \qquad C^2 = 1, \qquad Ct = tC.

    Therefore `S \mapsto s, T \mapsto t` is a projective representation
    of `SL(2, \ZZ)`. Let us confirm these identities for the Fibonacci MTC
    ``FusionRing("G2", 1)``::

        sage: R = FusionRing(("G2", 1))
        sage: S = R.s_matrix(unitary=True)
        sage: T = R.twists_matrix()
        sage: C = R.conj_matrix()
        sage: c = R.virasoro_central_charge(); c
        14/5
        sage: (S*T)^3 == R.root_of_unity(c/4) * S^2
        True
        sage: S^2 == C
        True
        sage: C*T == T*C
        True
    """
    def __init__(self, ct, k, conjugate=False, base_ring=ZZ, prefix=None, cyclotomic_order=None, fusion_labels=None, inject_variables=False):
        self._WCR = WeylCharacterRing(ct, base_ring, prefix,
                                        k=k, conjugate=conjugate,
                                        cyclotomic_order=cyclotomic_order,
                                        style="coroots")
        names = dict(self._WCR.basis())
        super().__init__(names=names, base_ring=base_ring,
                        prefix=prefix, conjugate=conjugate,
                        cyclotomic_order=cyclotomic_order,
                        fusion_labels=fusion_labels,
                        inject_variables=inject_variables)
        self._cyclotomic_order = self._WCR._cyclotomic_order

    def _repr_(self):
        """
        EXAMPLES::

            sage: FusionRing(("A1", 3))
            The Fusion Ring of Type A1 and level 3 with Integer Ring coefficients
        """
        return self._WCR._repr_()

    def _test_verlinde(self, **options):
        """
        Check the Verlinde formula for this :class:`FusionRing` instance.

        EXAMPLES::

            sage: G22 = FusionRing(("G2", 2))
            sage: G22._test_verlinde()
        """
        tester = self._tester(**options)
        c = self.global_q_dimension()
        i0 = self.one()
        from sage.misc.misc import some_tuples
        B = self.basis()
        for x, y, z in some_tuples(B, 3, tester._max_runs):
            v = sum(self.s_ij(x, w) * self.s_ij(y, w) * self.s_ij(z, w) / self.s_ij(i0, w) for w in B)
            tester.assertEqual(v, c * self.N_ijk(x, y, z))

    def _test_total_q_order(self, **options):
        r"""
        Check that the total quantum order is real and positive.

        The total quantum order is the positive square root
        of the global quantum dimension. This indirectly test the
        Virasoro central charge.

        EXAMPLES::

            sage: G22 = FusionRing(("G2", 2))
            sage: G22._test_total_q_order()
        """
        tester = self._tester(**options)
        tqo = self.total_q_order(base_coercion=False)
        tester.assertTrue(tqo.is_real_positive())
        tester.assertEqual(tqo**2, self.global_q_dimension(base_coercion=False))

    def cartan_type(self):
        """
        Return the Cartan type of ``self``.

        EXAMPLES::

            sage: G22 = FusionRing(("G2", 2))
            sage: G22.cartan_type()
            ['G', 2]
        """
        return self._WCR._cartan_type

    def fundamental_weights(self):
        """
        Return the fundamental weights of the :class:`WeylCharacterRing`
        associated to ``self``.

        EXAMPLES::

            sage: G22 = FusionRing(("G2", 2))
            sage: G22.fundamental_weights()
            Finite family {1: (1, 0, -1), 2: (2, -1, -1)}
        """
        return self._WCR.fundamental_weights()

    def fusion_level(self):
        r"""
        Return the level `k` of ``self``.

        EXAMPLES::

            sage: B22 = FusionRing(('B2', 2))
            sage: B22.fusion_level()
            2
        """
        return self._WCR._k

    def fusion_l(self):
        r"""
        Return the product `\ell = m_g(k + h^\vee)`, where `m_g` denotes the
        square of the ratio of the lengths of long to short roots of
        the underlying Lie algebra, `k` denotes the level of the FusionRing,
        and `h^\vee` denotes the dual Coxeter number of the underlying Lie
        algebra.

        This value is used to define the associated root `2\ell`-th
        of unity `q = e^{i\pi/\ell}`.

        EXAMPLES::

            sage: B22 = FusionRing(('B2', 2))
            sage: B22.fusion_l()
            10
            sage: D52 = FusionRing(('D5', 2))
            sage: D52.fusion_l()
            10
        """
        return self._WCR._l

    @cached_method
    def Nk_ij(self, elt_i, elt_j, elt_k):
        r"""
        Return the fusion coefficient `N^k_{ij}`.

        These are the structure coefficients of the fusion ring, so

        .. MATH::

            i * j = \sum_{k} N_{ij}^k k.

        EXAMPLES::

            sage: A22 = FusionRing(("A2", 2))
            sage: b = A22.basis().list()
            sage: all(x*y == sum(A22.Nk_ij(x, y, k)*k for k in b) for x in b for y in b)
            True
        """
        mc = (self._WCR(elt_i.weight()) * self._WCR(elt_j.weight()))._monomial_coefficients
        return mc.get(elt_k.weight(), 0)

    def one(self):
        return self.basis()[self._WCR.one_basis()]

    @cached_method
    def s_ij(self, elt_i, elt_j, unitary=False, base_coercion=True):
        r"""
        Return the element of the S-matrix of this fusion ring corresponding to
        the given elements.

        This is computed using the formula

        .. MATH::

            s_{i, j} = \frac{1}{\theta_i\theta_j} \sum_k N_{ik}^j d_k \theta_k,

        where `\theta_k` is the twist and `d_k` is the quantum
        dimension. See [Row2006]_ Equation (2.2) or [EGNO2015]_
        Proposition 8.13.8.

        INPUT:

        - ``elt_i``, ``elt_j`` -- elements of the fusion basis

        EXAMPLES::

            sage: G21 = FusionRing(("G2", 1))
            sage: b = G21.basis()
            sage: [G21.s_ij(x, y) for x in b for y in b]
            [1, -zeta60^14 + zeta60^6 + zeta60^4, -zeta60^14 + zeta60^6 + zeta60^4, -1]
        """
        ijtwist = elt_i.twist() + elt_j.twist()
        ret = sum(k.q_dimension(base_coercion=False) * self.Nk_ij(elt_i, k, elt_j)
                   * self.root_of_unity(k.twist() - ijtwist, base_coercion=False)
                   for k in self.basis())
        if unitary:
            ret /= self.total_q_order(base_coercion=False)
        if (not base_coercion) or (self._basecoer is None):
            return ret
        return self._basecoer(ret)

    def virasoro_central_charge(self):
        r"""
        Return the Virasoro central charge of the WZW conformal
        field theory associated with the Fusion Ring.

        If `\mathfrak{g}` is the corresponding semisimple Lie algebra, this is

        .. MATH::

            \frac{k\dim\mathfrak{g}}{k+h^\vee},

        where `k` is the level and `h^\vee` is the dual Coxeter number.
        See [DFMS1996]_ Equation (15.61).

        Let `d_i` and `\theta_i` be the quantum dimensions and
        twists of the simple objects. By Proposition 2.3 in [RoStWa2009]_,
        there exists a rational number `c` such that
        `D_+ / \sqrt{D} = e^{i\pi c/4}`, where `D_+ = \sum d_i^2 \theta_i`
        is computed in :meth:`D_plus` and `D = \sum d_i^2 > 0` is computed
        by :meth:`global_q_dimension`. Squaring this identity and
        remembering that `D_+ D_- = D` gives

        .. MATH::

            D_+ / D_- = e^{i\pi c/2}.

        EXAMPLES::

            sage: R = FusionRing(("A1", 2))
            sage: c = R.virasoro_central_charge(); c
            3/2
            sage: Dp = R.D_plus(); Dp
            2*zeta32^6
            sage: Dm = R.D_minus(); Dm
            -2*zeta32^10
            sage: Dp / Dm == R.root_of_unity(c/2)
            True
        """
        dim_g = len(self._WCR.space().roots()) + self._WCR.cartan_type().rank()
        return self._conj * self.fusion_level() * dim_g / (self.fusion_level() + self._WCR._h_check)

    class Element(FusionRing.Element):

        @cached_method
        def q_dimension(self, base_coercion=True):
            r"""
            Return the quantum dimension as an element of the cyclotomic
            field of the `2\ell`-th roots of unity, where `l = m (k+h^\vee)`
            with `m=1, 2, 3` depending on whether type is simply, doubly or
            triply laced, `k` is the level and `h^\vee` is the dual
            Coxeter number.
            """
            if not self.is_simple_object():
                raise ValueError("quantum dimension is only available for simple objects of a FusionRing")
            P = self.parent()._WCR
            lam = self.weight()
            space = P.space()
            rho = space.rho()
            powers = {}
            for alpha in space.positive_roots():
                val = alpha.inner_product(lam + rho)
                if val in powers:
                    powers[val] += 1
                else:
                    powers[val] = 1
                val = alpha.inner_product(rho)
                if val in powers:
                    powers[val] -= 1
                else:
                    powers[val] = -1
            R = ZZ['q']
            q = R.gen()
            expr = R.fraction_field().one()
            for val in powers:
                exp = powers[val]
                if exp > 0:
                    expr *= q_int(P._nf * val, q)**exp
                elif exp < 0:
                    expr /= q_int(P._nf * val, q)**(-exp)
            expr = R(expr)
            expr = expr.substitute(q=q**4) / (q**(2*expr.degree()))
            zet = self.parent().field().gen() ** (self.parent()._cyclotomic_order/P._l)
            ret = expr.substitute(q=zet)

            if (not base_coercion) or (self.parent()._basecoer is None):
                return ret
            return self.parent()._basecoer(ret)

        @cached_method
        def ribbon(self, base_coercion=True):
            r"""
            Return the twist or ribbon element of ``self``.

            If `h` is the rational number modulo 2 produced by
            ``self.twist()``, this  method produces `e^{i\pi h}`.

            .. SEEALSO::

                An additive version of this is available as :meth:`twist`.

            EXAMPLES::

                sage: F = FusionRing(("A1", 3))
                sage: [x.twist() for x in F.basis()]
                [0, 3/10, 4/5, 3/2]
                sage: [x.ribbon(base_coercion=False) for x in F.basis()]
                [1, zeta40^6, zeta40^12 - zeta40^8 + zeta40^4 - 1, -zeta40^10]
                sage: [F.root_of_unity(x, base_coercion=False) for x in [0, 3/10, 4/5, 3/2]]
                [1, zeta40^6, zeta40^12 - zeta40^8 + zeta40^4 - 1, -zeta40^10]
            """
            ret = self.parent().root_of_unity(self.twist(), base_coercion=False)
            if (not base_coercion) or (self.parent()._basecoer is None):
                return ret
            return self.parent()._basecoer(ret)

        def twist(self, reduced=True):
            r"""
            Return a rational number `h` such that `\theta = e^{i \pi h}`
            is the twist of ``self``. The quantity `e^{i \pi h}` is
            also available using :meth:`ribbon`.

            This method is only available for simple objects. If
            `\lambda` is the weight of the object, then
            `h = \langle \lambda, \lambda+2\rho \rangle`, where
            `\rho` is half the sum of the positive roots.
            As in [Row2006]_, this requires normalizing
            the invariant bilinear form so that
            `\langle \alpha, \alpha \rangle = 2` for short roots.

            INPUT:

            - ``reduced`` -- (default: ``True``) boolean; if ``True``
              then return the twist reduced modulo 2
            """
            if not self.is_simple_object():
                raise ValueError("Quantum twist is only available for simple objects of a FusionRing")
            P = self.parent()._WCR
            rho = P.space().rho()
            # We copy self.weight() to skip the test (which was already done
            # by self.is_simple_object()).
            lam = next(iter(self._monomial_coefficients))
            inner = lam.inner_product(lam + 2*rho)
            twist = P._conj * P._nf * inner / self.parent().fusion_l()
            # Reduce modulo 2
            if reduced:
                f = twist.floor()
                twist -= f
                return twist + (f % 2)
            return twist

        def weight(self):
            r"""
            Return the parametrizing dominant weight in the level `k` alcove.

            This method is only available for basis elements.

            EXAMPLES::

                sage: A21 = FusionRing(("A2", 1))
                sage: [x.weight() for x in A21.basis().list()]
                [(0, 0, 0), (2/3, -1/3, -1/3), (1/3, 1/3, -2/3)]
            """
            if len(self._monomial_coefficients) != 1:
                raise ValueError("Fusion weight is valid for basis elements only")
            return next(iter(self._monomial_coefficients))
