r"""
Quantum Oscillator Algebras

AUTHORS:

- Travis Scrimshaw (2023-12): initial version
"""

#*****************************************************************************
#  Copyright (C) 2023 Travis Scrimshaw <tcscrims at gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.misc.misc_c import prod
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.integer_ring import ZZ
from sage.categories.algebras import Algebras
from sage.combinat.free_module import CombinatorialFreeModule
from sage.categories.cartesian_product import cartesian_product
from sage.sets.family import Family
from sage.sets.non_negative_integers import NonNegativeIntegers


class QuantumOscillatorAlgebra(CombinatorialFreeModule):
    r"""
    The quantum oscillator algebra.

    Let `R` be a commutative algebra and `q \in R` be a unit.
    The *quantum oscillator algebra*, or `q`-oscillator algebra,
    is the unital associative `R`-algebra with generators `a^+`,
    `a^-` and `k^{\pm 1}` satisfying the following relations:

    .. MATH::

        k a^{\pm} = q^{\pm 1} a^{\pm} k, \qquad
        a^- a^+ = 1 - q^2 k^2, \qquad
        a^+ a^- = 1 - k^2.

    INPUT:

    - ``q`` -- (optional) the parameter `q`
    - ``R`` -- (default: `\QQ(q)`) the base ring that contains ``q``

    EXAMPLES:

    We construct the algebra and perform some basic computations::

        sage: O = algebras.QuantumOscillator()
        sage: ap, am, k, ki = O.algebra_generators()
        sage: q = O.q()
        sage: k^-3 * ap * ki * am^2 * k - q^3 * ap * k^3
        q^5*a-*k^-3 - q^3*a-*k^-1 - q^3*a+*k^3

    We construct representations of the type `A_1` quantum coordinate ring
    using the quantum oscillator algebra and verify the quantum determinant::

        sage: pi = matrix([[am, k], [-q*k, ap]]); pi
        [  a-    k]
        [-q*k   a+]
        sage: pi[0,0] * pi[1,1] - q * pi[0,1] * pi[1,0]
        1

    Next, we use this to build representations for type `A_2`::

        sage: def quantum_det(M):
        ....:     n = M.nrows()
        ....:     return sum((-q)**sigma.length()
        ....:                * prod(M[i,sigma[i]-1] for i in range(n))
        ....:                for sigma in Permutations(n))
        sage: def build_repr(wd, gens):
        ....:     n = gens[0].nrows()
        ....:     ret = gens[wd[0]-1]
        ....:     for ind in wd[1:]:
        ....:         g = gens[ind-1]
        ....:         temp = [[None]*n for _ in range(n)]
        ....:         for i in range(n):
        ....:             for j in range(n):
        ....:                 temp[i][j] = sum(tensor([ret[i,k], g[k,j]]) for k in range(n))
        ....:         ret = matrix(temp)
        ....:     return ret
        sage: pi1 = matrix.block_diagonal(pi, matrix.identity(1)); pi1
        [  a-    k|   0]
        [-q*k   a+|   0]
        [---------+----]
        [   0    0|   1]
        sage: pi2 = matrix.block_diagonal(matrix.identity(1), pi); pi2
        [   1|   0    0]
        [----+---------]
        [   0|  a-    k]
        [   0|-q*k   a+]
        sage: quantum_det(pi1) == 1
        True
        sage: quantum_det(pi2) == 1
        True
        sage: pi12 = build_repr([1,2], [pi1, pi2]); pi12
        [  a- # 1   k # a-    k # k]
        [-q*k # 1  a+ # a-   a+ # k]
        [       0 -q*1 # k   1 # a+]
        sage: quantum_det(pi12)
        1 # 1
        sage: pi121 = build_repr([1,2,1], [pi1, pi2]); pi121
        [   a- # 1 # a- - q*k # a- # k      a- # 1 # k + k # a- # a+                     k # k # 1]
        [-q*k # 1 # a- - q*a+ # a- # k   -q*k # 1 # k + a+ # a- # a+                    a+ # k # 1]
        [                q^2*1 # k # k                 -q*1 # k # a+                    1 # a+ # 1]
        sage: quantum_det(pi121)
        1 # 1 # 1
        sage: pi212 = build_repr([2,1,2], [pi1, pi2]); pi212
        [                   1 # a- # 1                    1 # k # a-                     1 # k # k]
        [                -q*a- # k # 1    a- # a+ # a- - q*k # 1 # k      a- # a+ # k + k # 1 # a+]
        [                q^2*k # k # 1 -q*k # a+ # a- - q*a+ # 1 # k   -q*k # a+ # k + a+ # 1 # a+]
        sage: quantum_det(pi212)
        1 # 1 # 1

    REFERENCES:

    - [Kuniba2022]_ Section 3.2
    """
    @staticmethod
    def __classcall_private__(cls, q=None, R=None):
        r"""
        Standardize input to ensure a unique representation.

        TESTS::

            sage: O1 = algebras.QuantumOscillator()
            sage: q = PolynomialRing(ZZ, 'q').fraction_field().gen()
            sage: O2 = algebras.QuantumOscillator(q=q)
            sage: O3 = algebras.QuantumOscillator(q, q.parent())
            sage: O1 is O2 and O2 is O3
            True
        """
        if q is None:
            q = PolynomialRing(ZZ, 'q').fraction_field().gen()
        if R is None:
            R = q.parent()
        q = R(q)

        return super().__classcall__(cls, q, R)

    def __init__(self, q, R):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: O = algebras.QuantumOscillator()
            sage: TestSuite(O).run()
        """
        self._q = q
        self._k_poly = PolynomialRing(R, 'k')
        indices = cartesian_product([ZZ, ZZ])

        cat = Algebras(R).WithBasis()
        CombinatorialFreeModule.__init__(self, R, indices, category=cat)
        self._assign_names(('ap', 'am', 'k', 'ki'))

    def _repr_(self) -> str:
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: algebras.QuantumOscillator()
            Quantum oscillator algebra with q=q over
             Fraction Field of Univariate Polynomial Ring in q over Integer Ring
        """
        return "Quantum oscillator algebra with q={} over {}".format(
            self._q, self.base_ring())

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: O = algebras.QuantumOscillator()
            sage: latex(O)
            \operatorname{Osc}_{q}
        """
        return "\\operatorname{Osc}_{%s}" % self._q

    def q(self):
        r"""
        Return the `q` of ``self``.

        EXAMPLES::

            sage: O = algebras.QuantumOscillator()
            sage: O.q()
            q
            sage: O = algebras.QuantumOscillator(q=QQ(-5))
            sage: O.q()
            -5
        """
        return self._q

    @cached_method
    def algebra_generators(self):
        r"""
        Return the algebra generators of ``self``.

        EXAMPLES::

            sage: O = algebras.QuantumOscillator()
            sage: O.algebra_generators()
            Finite family {'am': a-, 'ap': a+, 'k': k, 'ki': k^-1}
        """
        d = {'ap': self.monomial((ZZ.one(), ZZ.zero())),
             'am': self.monomial((-ZZ.one(), ZZ.zero())),
             'k': self.monomial((ZZ.zero(), ZZ.one())),
             'ki': self.monomial((ZZ.zero(), -ZZ.one()))}
        return Family(d)

    @cached_method
    def gens(self) -> tuple:
        r"""
        Return the generators of ``self``.

        EXAMPLES::

            sage: O = algebras.QuantumOscillator()
            sage: O.gens()
            (a+, a-, k, k^-1)
        """
        return tuple(self.algebra_generators())

    @cached_method
    def one_basis(self) -> tuple:
        r"""
        Return the index of the basis element of `1`.

        EXAMPLES::

            sage: O = algebras.QuantumOscillator()
            sage: O.one_basis()
            (0, 0)
        """
        return (ZZ.zero(), ZZ.zero())

    def some_elements(self) -> tuple:
        r"""
        Return some elements of ``self``.

        EXAMPLES::

            sage: O = algebras.QuantumOscillator()
            sage: O.some_elements()
            (a+, a-, k, k^-1, 1, a+^3, a-^4, k^2, k^-5, a+*k,
             a-^4*k^-3, 1 + 3*k + 2*a+ + a+*k)
        """
        ap, am, k, ki = self.gens()
        return (ap, am, k, ki, self.one(),
                ap**3, am**4, k**2, ki**5, ap*k, am**4*ki**3,
                self.an_element())

    def fock_space_representation(self):
        r"""
        Return the Fock space representation of ``self``.

        .. SEEALSO::

            :class:`~sage.algebras.quantum_oscillator.FockSpaceRepresentation`

        EXAMPLES::

            sage: O = algebras.QuantumOscillator()
            sage: O.fock_space_representation()
            Fock space representation of Quantum oscillator algebra with q=q
             over Fraction Field of Univariate Polynomial Ring in q over Integer Ring
        """
        return FockSpaceRepresentation(self)

    def _repr_term(self, m) -> str:
        r"""
        Return a string representation of the basis element indexed by ``m``.

        EXAMPLES::

            sage: O = algebras.QuantumOscillator()
            sage: O._repr_term((1, 3))
            'a+*k^3'
            sage: O._repr_term((-1, 1))
            'a-*k'
            sage: O._repr_term((5, 0))
            'a+^5'
            sage: O._repr_term((-4, -2))
            'a-^4*k^-2'
            sage: O._repr_term((0, -4))
            'k^-4'
            sage: O._repr_term((0, 0))
            '1'

            sage: O(5)
            5
        """
        a, k = m

        astr = ''
        if a == 1:
            astr = 'a+'
        elif a > 1:
            astr = 'a+^{}'.format(a)
        elif a == -1:
            astr = 'a-'
        elif a < -1:
            astr = 'a-^{}'.format(-a)

        kstr = ''
        if k == 1:
            kstr = 'k'
        elif k != 0:
            kstr = 'k^{}'.format(k)

        if astr:
            if kstr:
                return astr + '*' + kstr
            return astr
        if kstr:
            return kstr
        return '1'

    def _latex_term(self, m):
        r"""
        Return a latex representation for the basis element indexed by ``m``.

        EXAMPLES::

            sage: O = algebras.QuantumOscillator()
            sage: O._latex_term((1, 3))
            'a^+ k^{3}'
            sage: O._latex_term((-1, 1))
            'a^- k'
            sage: O._latex_term((5, 0))
            '(a^+)^{5}'
            sage: O._latex_term((-4, -2))
            '(a^-)^{4} k^{-2}'
            sage: O._latex_term((0, -4))
            'k^{-4}'
            sage: O._latex_term((0, 0))
            '1'

            sage: latex(O(5))
            5
        """
        a, k = m

        astr = ''
        if a == 1:
            astr = 'a^+'
        elif a > 1:
            astr = '(a^+)^{{{}}}'.format(a)
        elif a == -1:
            astr = 'a^-'
        elif a < -1:
            astr = '(a^-)^{{{}}}'.format(-a)

        kstr = ''
        if k == 1:
            kstr = 'k'
        elif k != 0:
            kstr = 'k^{{{}}}'.format(k)

        if astr:
            if kstr:
                return astr + ' ' + kstr
            return astr
        if kstr:
            return kstr
        return '1'

    @cached_method
    def product_on_basis(self, ml, mr):
        r"""
        Return the product of the basis elements indexed by ``ml`` and ``mr``.

        EXAMPLES::

            sage: O = algebras.QuantumOscillator()
            sage: ap, am, k, ki = O.algebra_generators()
            sage: O.product_on_basis((-2, 3), (-4, 5))
            1/q^12*a-^6*k^8
            sage: O.product_on_basis((2, 3), (4, -5))
            q^12*a+^6*k^-2
            sage: O.product_on_basis((2, 3), (0, -3))
            a+^2
            sage: k^5 * ki^10
            k^-5
            sage: k^10 * ki^5
            k^5
            sage: ap^3 * k^5
            a+^3*k^5
            sage: am^3 * k^5
            a-^3*k^5
            sage: k^5 * ap^3
            q^15*a+^3*k^5
            sage: k^5 * am^3
            1/q^15*a-^3*k^5
            sage: ki^5 * ap^3
            1/q^15*a+^3*k^-5
            sage: ki^5 * am^3
            q^15*a-^3*k^-5
            sage: ap * am
            1 - k^2
            sage: am * ap
            1 - q^2*k^2

            sage: (ap + am + k + ki)^2
            a-^2 + (q+1)*a-*k^-1 + ((q+1)/q)*a-*k + k^-2 + 4 - q^2*k^2
             + ((q+1)/q)*a+*k^-1 + (q+1)*a+*k + a+^2

            sage: (ap)^3 * (am)^5
            a-^2 + ((-q^4-q^2-1)/q^8)*a-^2*k^2 + ((q^4+q^2+1)/q^14)*a-^2*k^4 - 1/q^18*a-^2*k^6
            sage: (ap)^5 * (am)^3
            a+^2 + ((-q^4-q^2-1)/q^4)*a+^2*k^2 + ((q^4+q^2+1)/q^6)*a+^2*k^4 - 1/q^6*a+^2*k^6
            sage: (am)^3 * (ap)^5
            a+^2 + (-q^10-q^8-q^6)*a+^2*k^2 + (q^18+q^16+q^14)*a+^2*k^4 - q^24*a+^2*k^6
            sage: (am)^5 * (ap)^3
            a-^2 + (-q^6-q^4-q^2)*a-^2*k^2 + (q^10+q^8+q^6)*a-^2*k^4 - q^12*a-^2*k^6
        """
        q = self._q
        k = self._k_poly.gen()
        al, kl = ml
        ar, kr = mr
        coeff = q ** (kl * ar)
        if (al <= 0 and ar <= 0) or (al >= 0 and ar >= 0):
            return self.element_class(self, {(al + ar, kl + kr): coeff})
        # now al and ar have different signs
        if al < 0:  # a^- * a^+ case
            kp = self._k_poly.prod(1 - q**(2*(ar-i)) * k**2 for i in range(min(-al,ar)))
        else:  # a^+ * a^- case
            kp = self._k_poly.prod(1 - q**(2*(ar+i)) * k**2 for i in range(1,min(al,-ar)+1))
        a = al + ar
        return self.element_class(self, {(a, kl+kr+i): c * coeff for i, c in enumerate(kp) if c})

    class Element(CombinatorialFreeModule.Element):
        def __invert__(self):
            r"""
            Return the inverse if ``self`` is a basis element.

            EXAMPLES::

                sage: O = algebras.QuantumOscillator()
                sage: ap, am, k, ki = O.algebra_generators()
                sage: k.inverse()
                k^-1
                sage: ~k^5
                k^-5
                sage: ~ki^2
                k^2
                sage: O.zero().inverse()
                Traceback (most recent call last):
                ...
                ZeroDivisionError
                sage: ~ap
                Traceback (most recent call last):
                ...
                NotImplementedError: only implemented for monomials in k
                sage: ~(k + ki)
                Traceback (most recent call last):
                ...
                NotImplementedError: only implemented for monomials in k
            """
            if not self:
                raise ZeroDivisionError
            if len(self) != 1 or self.leading_support()[0] != 0:
                raise NotImplementedError("only implemented for monomials in k")

            ((a, k), coeff), = list(self._monomial_coefficients.items())
            O = self.parent()
            return O.element_class(O, {(a, -k): coeff.inverse_of_unit()})


class FockSpaceRepresentation(CombinatorialFreeModule):
    r"""
    The unique Fock space representation of the
    :class:`~sage.algebras.quantum_oscillator.QuantumOscillatorAlgebra`.
    """
    def __init__(self, oscillator_algebra):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: O = algebras.QuantumOscillator()
            sage: F = O.fock_space_representation()
            sage: TestSuite(F).run()
        """
        self._O = oscillator_algebra
        ind = NonNegativeIntegers()
        CombinatorialFreeModule.__init__(self, oscillator_algebra.base_ring(), ind, prefix='', bracket=['|', '>'],
                                         latex_bracket=[r'\lvert', r'\rangle'])

    def _test_representation(self, **options):
        r"""
        Test that ``self`` is a representation of the quantum
        oscillator algebra.

        EXAMPLES::

            sage: O = algebras.QuantumOscillator(q=GF(7)(3))
            sage: F = O.fock_space_representation()
            sage: F._test_representation()
        """
        tester = self._tester(**options)
        S = self._O.some_elements()
        num_trials = 0
        from itertools import product
        for a, b in product(S, repeat=2):
            for elt in tester.some_elements():
                num_trials += 1
                if num_trials > tester._max_runs:
                    return
                tester.assertEqual((a*b)*elt, a*(b*elt))

    def _repr_(self) -> str:
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: O = algebras.QuantumOscillator(q=GF(5)(2))
            sage: O.fock_space_representation()
            Fock space representation of Quantum oscillator algebra
             with q=2 over Finite Field of size 5
        """
        return "Fock space representation of {}".format(self._O)

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: O = algebras.QuantumOscillator()
            sage: F = O.fock_space_representation()
            sage: latex(F)
            \mathfrak{F}_{q}
        """
        return r"\mathfrak{{F}}_{{{}}}".format(self._O._q)

    def vacuum(self):
        r"""
        Return the vacuum element `|0\rangle` of ``self``.

        EXAMPLES::

            sage: O = algebras.QuantumOscillator()
            sage: F = O.fock_space_representation()
            sage: F.vacuum()
            |0>
        """
        return self.basis()[0]

    def some_elements(self):
        r"""
        Return some elements of ``self``.

        EXAMPLES::

            sage: O = algebras.QuantumOscillator()
            sage: F = O.fock_space_representation()
            sage: F.some_elements()
            (|0>, |1>, |52>, |0> + 2*|1> + 3*|2> + |42>)
        """
        B = self.basis()
        return (B[0], B[1], B[52], self.an_element())

    class Element(CombinatorialFreeModule.Element):
        def _acted_upon_(self, scalar, self_on_left=True):
            r"""
            Return the action of ``scalar`` on ``self``.

            EXAMPLES::

                sage: O = algebras.QuantumOscillator()
                sage: ap, am, k, ki = O.gens()
                sage: F = O.fock_space_representation()
                sage: B = F.basis()
                sage: [ap * B[i] for i in range(3)]
                [|1>, |2>, |3>]
                sage: [am * B[i] for i in range(3)]
                [0, (-q^2+1)*|0>, (-q^4+1)*|1>]
                sage: [k * B[i] for i in range(3)]
                [|0>, q*|1>, q^2*|2>]
                sage: [ki * B[i] for i in range(3)]
                [|0>, 1/q*|1>, 1/q^2*|2>]
                sage: (am)^3 * B[5]
                (-q^24+q^18+q^16+q^14-q^10-q^8-q^6+1)*|2>
                sage: (7*k^3 + am) * (B[0] + B[1] + B[2])
                (-q^2+8)*|0> + (-q^4+7*q^3+1)*|1> + 7*q^6*|2>
                sage: 5 * (B[2] + B[3])
                5*|2> + 5*|3>
            """
            # Check for scalars first
            ret = super()._acted_upon_(scalar, self_on_left)
            if ret is not None:
                return ret
            P = self.parent()
            if self_on_left or scalar not in P._O:  # needs to be a left Osc-action
                return None
            scalar = P._O(scalar)
            q = P._O._q

            ret = []
            for om, oc in scalar:
                a, k = om
                for fm, fc in self:
                    if fm < -a:  # the result will be 0
                        continue
                    c = q ** (fm*k)
                    if a < 0:
                        c *= prod(1 - q**(2*(fm-i)) for i in range(-a))
                    if c:
                        ret.append((fm+a, oc * fc * c))
            return P.sum_of_terms(ret)
