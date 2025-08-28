r"""
Ariki-Koike Algebra Representations

AUTHORS:

- Travis Scrimshaw (2023-12-28): Initial version
"""

#*****************************************************************************
#       Copyright (C) 2023 Travis Scrimshaw <tcscrims at gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.misc_c import prod
from sage.misc.latex import latex
from sage.categories.modules import Modules
from sage.rings.integer_ring import ZZ
from sage.combinat.free_module import CombinatorialFreeModule
from sage.combinat.partition_tuple import PartitionTuples
from sage.combinat.permutation import Permutations
from sage.combinat.tableau_tuple import StandardTableauTuples


class SpechtModule(CombinatorialFreeModule):
    r"""
    Specht module of the Ariki-Koike algebra.

    Let `H_{r,n}(q, u)` be the Ariki-Koike algebra with parameters `q`
    and `u = (u_1, \ldots, u_r)` (note our indexing convention for
    the `u` parameters differs from
    :mod:`sage.algebras.hecke_algebras.ariki_koike_algebra`) over a
    commutative ring `R`. Let `\lambda` be a partition tuple of level
    `r` and size `n`. The *Specht module* of shape `\lambda` is the (right)
    `H_{r,n}(q,u)`-representation `S^{\lambda}` given as free `R`-module
    with basis given by the standard tableau (tuples) of shape `\lambda`.

    We will now describe the right action of the Ariki-Koike algebra,
    but we first need to set some notation and definitions.
    Let `t` be a standard tableau tuple of level `r` and size `n`.
    Define the *residue* of `i` in `t` to be `r_t(i) = q^{c-r} u_k`,
    where `i` is in cell `(r, c)` of the `k`-th tableau.

    The action of `L_i` is given by `t \cdot L_i = r_T(i) t`. For `T_i`,
    we need to consider the following cases. If `i, i+1` are in the same
    row (resp. column), then `t \cdot T_i = q t` (resp. `t \cdot T_i = -t`).
    Otherwise if we swap `i, i+1`, the resulting tableau tuple `s` is again
    standard and the action is given by

    .. MATH::

        t \cdot T_i = \frac{(q - 1) r_t(i)}{r_s(i) - r_t(i)} t
        + \frac{q r_t(i) - r_s(i)}{r_s(i) - r_t(i)} s.

    Note that `r_s(i) = r_t(i+1)`.

    Over a field of characteristic `0`, the set of Specht modules for all
    partition tuples of level `r` and size `n` form the complete set
    of irreducible modules for `H_{r,n}(q, u)` [AK1994]_. (The condition
    on the base ring can be weakened; see Theorem 3.2 of [Mathas2002]_.)

    EXAMPLES:

    We construct the Specht module `S^{(2,1,21)}` for `H_{3,6}(q, u)` with
    generic parameters `q, u` over `\GF(3)` and perform some basic
    computations. We change the tableaux to use the compact representation
    to condense the output::

        sage: TableauTuples.options.display = 'compact'

        sage: R = PolynomialRing(GF(3), 'u', 3)
        sage: u = R.gens()
        sage: q = R['q'].gen()
        sage: H = algebras.ArikiKoike(3, 6, q, u, use_fraction_field=True)
        sage: LT = H.LT()
        sage: T0, T1, T2, T3, T4, T5 = LT.T()
        sage: S = H.specht_module([[2], [1], [2,1]])
        sage: S.dimension()
        120
        sage: elt = S.an_element(); elt
        S[1,2|3|4,5/6] - S[1,3|2|4,5/6] + S[1,3|4|2,5/6]
        sage: elt * LT.L(3)
        u1*S[1,2|3|4,5/6] + (-u0*q)*S[1,3|2|4,5/6] + u0*q*S[1,3|4|2,5/6]
        sage: elt * T2
        (((-u0-u1)*q-u1)/(-u0*q+u1))*S[1,2|3|4,5/6]
         + (((-u0+u2)*q)/(u0*q-u2))*S[1,2|4|3,5/6]
         + ((-u0*q^2-u0*q-u1)/(-u0*q+u1))*S[1,3|2|4,5/6]
         + ((u0*q^2-u0*q)/(u0*q-u2))*S[1,3|4|2,5/6]
        sage: (elt * T3) * T2 == elt * (T3 * T2)
        True
        sage: elt * T2 * T3 * T2 == elt * T3 * T2 * T3
        True
        sage: elt * T0 * T1 * T0 * T1 == elt * T1 * T0 * T1 * T0
        True
        sage: elt * T2 * T5 == elt * T5 * T2
        True

        sage: TableauTuples.options._reset()

    REFERENCES:

    - [AK1994]_
    - [DJM1998]_
    - [DR2001]_
    - [Mathas2002]_
    - [Mathas2004]_
    """
    @staticmethod
    def __classcall_private__(cls, AK, la):
        """
        Normalize input to ensure a unique representation.

        EXAMPLES::

            sage: AK = algebras.ArikiKoike(3, 6)
            sage: S1 = AK.specht_module([[3], [1], [1,1]])
            sage: S2 = AK.specht_module(PartitionTuple([[3], [1], [1,1]]))
            sage: S1 is S2
            True
        """
        la = PartitionTuples(AK._r, AK._n)(la)
        return super().__classcall__(cls, AK, la)

    def __init__(self, AK, la):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: AK = algebras.ArikiKoike(3, 6, use_fraction_field=True)
            sage: S = AK.specht_module([[3], [1], [1,1]])
            sage: TestSuite(S).run()  # long time
            sage: Sp = AK.specht_module([[], [2,1,1], [2]])
            sage: TestSuite(Sp).run()  # long time
        """
        self._shape = la
        self._AK = AK
        self._q = AK.q()
        self._u = AK.u()
        self._Pn = Permutations(la.size())
        indices = StandardTableauTuples(la)
        R = AK.base_ring()
        cat = Modules(R).FiniteDimensional().WithBasis()
        CombinatorialFreeModule.__init__(self, R, indices, category=cat, prefix='S')

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: AK = algebras.ArikiKoike(3, 8)
            sage: AK.specht_module([[3], [], [2,2,1]])
            Specht module of shape ([3], [], [2, 2, 1]) for
             Ariki-Koike algebra of rank 3 and order 8 with q=q and u=(u0, u1, u2)
              over ... over Integer Ring
        """
        return "Specht module of shape {} for {}".format(self._shape, self._AK)

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: AK = algebras.ArikiKoike(3, 8)
            sage: S = AK.specht_module([[3], [], [2,2,1]])
            sage: latex(S)
            S^{{\def\lr#1{\multicolumn{1}{|@{\hspace{.6ex}}c@{\hspace{.6ex}}|}{\raisebox{-.3ex}{$#1$}}}
            \raisebox{-.6ex}{$\begin{array}[b]{*{3}c}\cline{1-3}
            \lr{\phantom{x}}&\lr{\phantom{x}}&\lr{\phantom{x}}\\\cline{1-3}
            \end{array}$},\emptyset,\raisebox{-.6ex}{$\begin{array}[b]{*{2}c}\cline{1-2}
            \lr{\phantom{x}}&\lr{\phantom{x}}\\\cline{1-2}
            \lr{\phantom{x}}&\lr{\phantom{x}}\\\cline{1-2}
            \lr{\phantom{x}}\\\cline{1-1}
            \end{array}$}
            }}_{\mathcal{H}_{3,8}(q)}
        """
        return "S^{{{}}}_{{{}}}".format(latex(self._shape), latex(self._AK))

    def _test_representation(self, **options):
        r"""
        Test that the relations of the Ariki-Koike algebra are satisfied.

        EXAMPLES::

            sage: q = ZZ['q'].fraction_field().gen()
            sage: AK = algebras.ArikiKoike(2, 4, q, [q^2+1, q-3], q.parent())
            sage: S = AK.specht_module([[2,1], [1]])
            sage: S._test_representation(elements=S.basis())
        """
        tester = self._tester(**options)
        n = self._shape.size()
        q = self._q
        from sage.misc.misc import some_tuples
        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

        # Build the polynomial for testing the T0 action
        z = PolynomialRing(self.base_ring(), 'DUMMY').gen()
        T0_poly = -prod(z - val for val in self._u)

        def apply_T0_power(b, exp):
            for i in range(exp):
                b = b.T(0)
            return b

        AKelts = self._AK.some_elements()
        for b in tester.some_elements():
            t0 = self.linear_combination((apply_T0_power(b, exp), c)
                                         for exp, c in enumerate(T0_poly))
            tester.assertEqual(t0, self.zero())

            tester.assertEqual(b.T([0, 1, 0, 1]), b.T([1, 0, 1, 0]))
            tester.assertEqual(b.T(1).T(1), (q-1)*b.T(1) + q*b)
            for i in range(2, n):
                tester.assertEqual(b.T(i).T(i), (q-1)*b.T(i) + q*b)
                tester.assertEqual(b.T(i).T(0), b.T(0).T(i))
                if i < n - 1:
                    tester.assertEqual(b.T([i, i+1, i]), b.T([i+1, i, i+1]))
                    for j in range(i+2, n):
                        tester.assertEqual(b.T([i, j]), b.T([j, i]))

            for (x, y) in some_tuples(AKelts, 2, tester._max_runs):
                tester.assertEqual(b*(x*y), (b*x)*y)

    def _L_on_basis(self, i, t):
        """
        Return the action of `L_i` on the basis element indexed by
        the standard tableau tuple ``t``.

        EXAMPLES::

            sage: AK = algebras.ArikiKoike(3, 10)
            sage: S = AK.specht_module([[2,1], [], [3,2,2]])
            sage: P = S.basis().keys()
            sage: t = P([[[2,4],[8]], [], [[1,3,7],[5,6],[9,10]]])
            sage: S._L_on_basis(1, t)
            u2*S[([[2, 4], [8]], [], [[1, 3, 7], [5, 6], [9, 10]])]
            sage: S._L_on_basis(4, t)
            u0*q*S[([[2, 4], [8]], [], [[1, 3, 7], [5, 6], [9, 10]])]
            sage: S._L_on_basis(6, t)
            u2*S[([[2, 4], [8]], [], [[1, 3, 7], [5, 6], [9, 10]])]
            sage: S._L_on_basis(8, t)
            (u0*q^-1)*S[([[2, 4], [8]], [], [[1, 3, 7], [5, 6], [9, 10]])]
            sage: S._L_on_basis(9, t)
            (u2*q^-2)*S[([[2, 4], [8]], [], [[1, 3, 7], [5, 6], [9, 10]])]
        """
        c = t.cells_containing(i)[0]
        if len(c) == 2:  # it is of level 1 and a regular tableau
            c = (0,) + c
        res = self._q**(c[2]-c[1]) * self._u[c[0]]
        R = self.base_ring()
        return self.element_class(self, {t: R(res)})

    def _T_on_basis(self, i, t):
        r"""
        Return the action of `T_i` on the basis element indexed by
        the standard tableau tuple ``t``.

        EXAMPLES::

            sage: AK = algebras.ArikiKoike(3, 10, use_fraction_field=True)
            sage: S = AK.specht_module([[2,1], [], [3,2,2]])
            sage: P = S.basis().keys()
            sage: t = P([[[2,4],[8]], [], [[1,5,7],[3,6],[9,10]]])
            sage: S._T_on_basis(0, t) == S._L_on_basis(1, t)
            True
            sage: S._T_on_basis(1, t)
            ((u2*q-u0)/(u0-u2))*S[([[1, 4], [8]], [], [[2, 5, 7], [3, 6], [9, 10]])]
             + ((u0*q-u0)/(u0-u2))*S[([[2, 4], [8]], [], [[1, 5, 7], [3, 6], [9, 10]])]
            sage: S._T_on_basis(2, t)
            ((u2*q-u2)/(-u0*q+u2))*S[([[2, 4], [8]], [], [[1, 5, 7], [3, 6], [9, 10]])]
             + ((u0*q^2-u2)/(-u0*q+u2))*S[([[3, 4], [8]], [], [[1, 5, 7], [2, 6], [9, 10]])]
            sage: S._T_on_basis(5, t)
            -S[([[2, 4], [8]], [], [[1, 5, 7], [3, 6], [9, 10]])]
            sage: S._T_on_basis(7, t)
            ((u2*q^4-u0)/(-u2*q^3+u0))*S[([[2, 4], [7]], [], [[1, 5, 8], [3, 6], [9, 10]])]
             + ((u0*q-u0)/(-u2*q^3+u0))*S[([[2, 4], [8]], [], [[1, 5, 7], [3, 6], [9, 10]])]
            sage: S._T_on_basis(9, t)
            q*S[([[2, 4], [8]], [], [[1, 5, 7], [3, 6], [9, 10]])]
        """
        R = self.base_ring()
        if i == 0:
            return self._L_on_basis(1, t)

        ct = t.cells_containing(i)[0]
        cs = t.cells_containing(i+1)[0]
        if len(ct) == 2:  # it is of level 1 and a regular tableau
            ct = (0,) + ct
            cs = (0,) + cs

        if ct[0] == cs[0] and ct[2] == cs[2]:  # same column
            return self.element_class(self, {t: -R.one()})

        if ct[0] == cs[0] and ct[1] == cs[1]:  # same row
            return self.element_class(self, {t: self._q})

        # result is standard
        s = t.symmetric_group_action_on_entries(self._Pn.simple_reflection(i))
        assert s.parent() is t.parent()

        def res(cell):
            return self._q**(cell[2]-cell[1]) * self._u[cell[0]]

        # Note that the residue of i in t is given by the cell c
        #   and of i in s corresponds to cell cp because the
        #   corresponding action of the permutation on t.
        one = self.base_ring().one()
        denom = res(cs) - res(ct)
        coefft = (self._q - one) * res(cs) / denom
        coeffs = (self._q * res(ct) - res(cs)) / denom
        return self.element_class(self, {t: R(coefft), s: R(coeffs)})

    def ariki_koike_algebra(self):
        r"""
        Return the Ariki-Koike algebra that ``self`` is a representation of.

        EXAMPLES::

            sage: AK = algebras.ArikiKoike(3, 6)
            sage: S = AK.specht_module([[2], [], [3,1]])
            sage: S.ariki_koike_algebra() is AK
            True
        """
        return self._AK

    class Element(CombinatorialFreeModule.Element):
        def _acted_upon_(self, scalar, self_on_left):
            r"""
            Return the action of ``scalar`` on ``self``.

            EXAMPLES::

                sage: TableauTuples.options.display = 'compact'
                sage: AK = algebras.ArikiKoike(4, 6, use_fraction_field=True)
                sage: q = AK.q()
                sage: LT = AK.LT()
                sage: T = AK.T()
                sage: S = AK.specht_module([[], [2], [1], [2,1]])
                sage: elt = S.an_element()
                sage: 5 * elt
                5*S[-|1,2|3|4,5/6] + 10*S[-|1,3|2|4,5/6] + 5*S[-|1,3|4|2,5/6]
                 + 15*S[-|2,3|1|4,5/6]
                sage: elt * (q - 2)
                (q-2)*S[-|1,2|3|4,5/6] + (2*q-4)*S[-|1,3|2|4,5/6]
                 + (q-2)*S[-|1,3|4|2,5/6] + (3*q-6)*S[-|2,3|1|4,5/6]
                sage: elt * LT.an_element() == elt * T(LT.an_element())
                True
                sage: T.an_element() * elt
                Traceback (most recent call last):
                ...
                TypeError: unsupported operand parent(s) for *: 'Ariki-Koike algebra ... over Integer Ring'
                sage: TableauTuples.options._reset()

            TESTS::

                sage: AK = algebras.ArikiKoike(2, 4, use_fraction_field=True)
                sage: LT = AK.LT()
                sage: T = AK.T()
                sage: S = AK.specht_module([[1], [2,1]])
                sage: B = list(LT.basis())[::55]
                sage: all(b * x == b * T(x) for b in S.basis() for x in B)  # long time
                True
            """
            ret = super()._acted_upon_(scalar, self_on_left)
            if ret is not None:
                return ret
            if not self_on_left:  # only a right action
                return None
            P = self.parent()
            if scalar not in P._AK:
                return None
            scalar = P._AK(scalar)
            if scalar.parent() is P._AK.LT():
                return P.linear_combination((self.L(sum(([i]*val for i, val in enumerate(m[0], start=1)), [])).T(m[1].reduced_word()), c)
                                            for m, c in scalar)
            elif scalar.parent() is P._AK.T():
                AKT = P._AK.T()
                return P.linear_combination((self.T(AKT._basis_to_word(m)), c)
                                            for m, c in scalar)
            return self * P._AK.LT()(scalar)

        def L(self, i):
            r"""
            Return the (right) action of `L_i` on ``self``.

            INPUT:

            - ``i`` -- an integer or a list of integers

            EXAMPLES::

                sage: TableauTuples.options.display = 'compact'  # compact tableau printing
                sage: AK = algebras.ArikiKoike(3, 6, use_fraction_field=True)
                sage: S = AK.specht_module([[2], [], [3,1]])
                sage: elt = S.an_element(); elt
                S[1,2|-|3,4,5/6] + 2*S[1,3|-|2,4,5/6]
                 + 3*S[2,3|-|1,4,5/6] + S[2,4|-|1,3,5/6]
                sage: elt.L(1)
                u0*S[1,2|-|3,4,5/6] + 2*u0*S[1,3|-|2,4,5/6]
                 + 3*u2*S[2,3|-|1,4,5/6] + u2*S[2,4|-|1,3,5/6]
                sage: elt.L(2)
                u0*q*S[1,2|-|3,4,5/6] + 2*u2*S[1,3|-|2,4,5/6]
                 + 3*u0*S[2,3|-|1,4,5/6] + u0*S[2,4|-|1,3,5/6]
                sage: elt.L(6)
                u2/q*S[1,2|-|3,4,5/6] + 2*u2/q*S[1,3|-|2,4,5/6]
                 + 3*u2/q*S[2,3|-|1,4,5/6] + u2/q*S[2,4|-|1,3,5/6]
                sage: elt.L([3,3,3])
                u2^3*S[1,2|-|3,4,5/6] + 2*u0^3*q^3*S[1,3|-|2,4,5/6]
                + 3*u0^3*q^3*S[2,3|-|1,4,5/6] + u2^3*q^3*S[2,4|-|1,3,5/6]
                sage: LT = AK.LT()
                sage: elt.L([3,3,3]) == elt * (LT.L(3)^3)
                True
                sage: TableauTuples.options._reset()  # reset
            """
            if not self:  # action on 0 is 0
                return self
            if i not in ZZ:
                ret = self
                for val in i:
                    ret = ret.L(val)
                return ret
            P = self.parent()
            return P.linear_combination((P._L_on_basis(i, t), c) for t, c in self)

        def T(self, i):
            r"""
            Return the (right) action of `T_i` on ``self``.

            INPUT:

            - ``i`` -- an integer or a list of integers

            EXAMPLES::

                sage: TableauTuples.options.display = 'compact'  # compact tableau printing
                sage: AK = algebras.ArikiKoike(3, 10, use_fraction_field=True)
                sage: q = AK.q()
                sage: S = AK.specht_module([[2,1], [], [3,2,2]])
                sage: P = S.basis().keys()
                sage: t = P([[[2,4],[8]], [], [[1,5,7],[3,6],[9,10]]])
                sage: b = S.basis()[t]
                sage: b.T(2)
                ((u2*q-u2)/(-u0*q+u2))*S[2,4/8|-|1,5,7/3,6/9,10]
                 + ((u0*q^2-u2)/(-u0*q+u2))*S[3,4/8|-|1,5,7/2,6/9,10]
                sage: b.T(6)
                ((-q)/(q+1))*S[2,4/8|-|1,5,6/3,7/9,10]
                 + (q^2/(q+1))*S[2,4/8|-|1,5,7/3,6/9,10]
                sage: b.T([2,1,2]) == b.T([1,2,1])
                True
                sage: b.T(9)
                q*S[2,4/8|-|1,5,7/3,6/9,10]
                sage: all(b.T([i,i]) == (q-1)*b.T(i) + q*b for i in range(1,10))
                True
                sage: b.T(0)
                u2*S[2,4/8|-|1,5,7/3,6/9,10]
                sage: b.T([0,1,0,1]) == b.T([1,0,1,0])
                True
                sage: TableauTuples.options._reset()  # reset
            """
            if not self:  # action on 0 is 0
                return self
            if i not in ZZ:
                ret = self
                for val in i:
                    ret = ret.T(val)
                return ret
            P = self.parent()
            return P.linear_combination((P._T_on_basis(i, t), c) for t, c in self)
