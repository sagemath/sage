# sage.doctest: needs sage.modules
r"""
Quantum Clifford Algebras

AUTHORS:

- Travis Scrimshaw (2021-05): initial version
"""

#*****************************************************************************
#  Copyright (C) 2021 Travis Scrimshaw <tcscrims at gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.integer_ring import ZZ
from sage.categories.algebras import Algebras
from sage.categories.fields import Fields
from sage.combinat.free_module import CombinatorialFreeModule
from sage.categories.cartesian_product import cartesian_product
from sage.sets.family import Family
from sage.rings.fraction_field import FractionField
from sage.sets.finite_enumerated_set import FiniteEnumeratedSet
from itertools import product
from sage.combinat.subset import powerset


class QuantumCliffordAlgebra(CombinatorialFreeModule):
    r"""
    The quantum Clifford algebra.

    The *quantum Clifford algebra*, or `q`-Clifford algebra,
    of rank `n` and twist `k` is the unital associative algebra
    `\mathrm{Cl}_{q}(n, k)` over a field `F` with generators
    `\psi_a, \psi_a^*, \omega_a` for `a = 1, \ldots, n` that
    satisfy the following relations:

    .. MATH::

        \begin{aligned}
        \omega_a \omega_b & = \omega_b \omega_a,
        & \omega_a^{4k} & = (1 + q^{-2k}) \omega_a^{2k} - q^{-2k},
        \\ \omega_a \psi_b & = q^{\delta_{ab}} \psi_b \omega_a,
        & \omega_a \psi^*_b & = \psi^*_b \omega_a,
        \\ \psi_a \psi_b & + \psi_b \psi_a = 0,
        & \psi^*_a \psi^*_b & + \psi^*_b \psi^*_a = 0,
        \\ \psi_a \psi^*_a & + q^k \psi^*_a \psi_a = \omega_a^{-k},
        & \psi^*_a \psi_a & + q^{-k} \psi^*_a \psi_a = \omega_a^k,
        \\ \psi_a \psi^*_b & + \psi_b^* \psi_a = 0
        & & \text{if } a \neq b.
        \end{aligned}

    When `k = 2`, we recover the original definition given by Hayashi in
    [Hayashi1990]_. The `k = 1` version was used in [Kwon2014]_.

    INPUT:

    - ``n`` -- positive integer; the rank
    - ``k`` -- positive integer (default: 1); the twist
    - ``q`` -- (optional) the parameter `q`
    - ``F`` -- (default: `\QQ(q)`) the base field that contains ``q``

    EXAMPLES:

    We construct the rank 3 and twist 1 `q`-Clifford algebra::

        sage: Cl = algebras.QuantumClifford(3)
        sage: Cl
        Quantum Clifford algebra of rank 3 and twist 1 with q=q over
         Fraction Field of Univariate Polynomial Ring in q over Integer Ring
        sage: q = Cl.q()

    Some sample computations::

        sage: p0, p1, p2, d0, d1, d2, w0, w1, w2 = Cl.gens()
        sage: p0 * p1
        psi0*psi1
        sage: p1 * p0
        -psi0*psi1
        sage: p0 * w0 * p1 * d0 * w2
        (1/(q^3-q))*psi1*w2 + (-q/(q^2-1))*psi1*w0^2*w2
        sage: w0^4
        -1/q^2 + ((q^2+1)/q^2)*w0^2

    We construct the homomorphism from `U_q(\mathfrak{sl}_3)` to
    `\mathrm{Cl}(3, 1)` given in (3.17) of [Hayashi1990]_::

        sage: e1 = p0*d1; e2 = p1*d2
        sage: f1 = p1*d0; f2 = p2*d1
        sage: k1 = w0*~w1; k2 = w1*~w2
        sage: k1i = w1*~w0; k2i = w2*~w1
        sage: (e1, e2, f1, f2, k1, k2, k1i, k2i)
        (psi0*psid1, psi1*psid2,
         -psid0*psi1, -psid1*psi2,
         (q^2+1)*w0*w1 - q^2*w0*w1^3, (q^2+1)*w1*w2 - q^2*w1*w2^3,
         (q^2+1)*w0*w1 - q^2*w0^3*w1, (q^2+1)*w1*w2 - q^2*w1^3*w2)

    We check that `k_i` and `k_i^{-1}` are inverses::

        sage: k1 * k1i
        1
        sage: k2 * k2i
        1

    The relations between `e_i`, `f_i`, and `k_i`::

        sage: k1 * f1 == q^-2 * f1 * k1
        True
        sage: k2 * f1 == q^1 * f1 * k2
        True
        sage: k2 * e1 == q^-1 * e1 * k2
        True
        sage: k1 * e1 == q^2 * e1 * k1
        True
        sage: e1 * f1 - f1 * e1 == (k1 - k1i)/(q-q^-1)
        True
        sage: e2 * f1 - f1 * e2
        0

    The `q`-Serre relations::

        sage: e1 * e1 * e2 - (q^1 + q^-1) * e1 * e2 * e1 + e2 * e1 * e1
        0
        sage: f1 * f1 * f2 - (q^1 + q^-1) * f1 * f2 * f1 + f2 * f1 * f1
        0

    This also can be constructed at the special point when `q^{2k} = 1`,
    but the basis used is different::

        sage: Cl = algebras.QuantumClifford(1, 1, -1)
        sage: Cl.inject_variables()
        Defining psi0, psid0, w0
        sage: psi0 * psid0
        psi0*psid0
        sage: psid0 * psi0
        -w0 + psi0*psid0
        sage: w0^2
        1
    """
    @staticmethod
    def __classcall_private__(cls, n, k=1, q=None, F=None):
        r"""
        Standardize input to ensure a unique representation.

        TESTS::

            sage: Cl1 = algebras.QuantumClifford(3)
            sage: q = PolynomialRing(ZZ, 'q').fraction_field().gen()
            sage: Cl2 = algebras.QuantumClifford(3, q=q)
            sage: Cl3 = algebras.QuantumClifford(3, 1, q, q.parent())
            sage: Cl1 is Cl2 and Cl2 is Cl3
            True
        """
        if q is None:
            q = PolynomialRing(ZZ, 'q').fraction_field().gen()
        if F is None:
            F = q.parent()
        q = F(q)
        if F not in Fields():
            F = FractionField(F)
            q = F(q)

        if bool(q**(2*k) == 1):
            return QuantumCliffordAlgebraRootUnity(n, k, q, F)
        return QuantumCliffordAlgebraGeneric(n, k, q, F)

    def __init__(self, n, k, q, F, psi, indices):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: Cl = algebras.QuantumClifford(1, 2)
            sage: TestSuite(Cl).run(elements=Cl.basis())
        """
        self._n = n
        self._k = k
        self._q = q
        self._w_poly = PolynomialRing(F, n, 'w')
        self._psi = psi
        indices = FiniteEnumeratedSet(indices)

        cat = Algebras(F).FiniteDimensional().Semisimple().WithBasis()
        CombinatorialFreeModule.__init__(self, F, indices, category=cat)
        self._assign_names(self.algebra_generators().keys())

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: algebras.QuantumClifford(3)
            Quantum Clifford algebra of rank 3 and twist 1 with q=q over
             Fraction Field of Univariate Polynomial Ring in q over Integer Ring
        """
        return "Quantum Clifford algebra of rank {} and twist {} with q={} over {}".format(
            self._n, self._k, self._q, self.base_ring())

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: Cl = algebras.QuantumClifford(3)
            sage: latex(Cl)
            \operatorname{Cl}_{q}(3, 1)
        """
        return "\\operatorname{Cl}_{%s}(%s, %s)" % (self._q, self._n, self._k)

    def q(self):
        r"""
        Return the `q` of ``self``.

        EXAMPLES::

            sage: Cl = algebras.QuantumClifford(3)
            sage: Cl.q()
            q

            sage: Cl = algebras.QuantumClifford(3, q=QQ(-5))
            sage: Cl.q()
            -5
        """
        return self._q

    def twist(self):
        r"""
        Return the twist `k` of ``self``.

        EXAMPLES::

            sage: Cl = algebras.QuantumClifford(3, 2)
            sage: Cl.twist()
            2
        """
        return self._k

    def rank(self):
        r"""
        Return the rank `k` of ``self``.

        EXAMPLES::

            sage: Cl = algebras.QuantumClifford(3, 2)
            sage: Cl.rank()
            3
        """
        return self._n

    def dimension(self):
        r"""
        Return the dimension of ``self``.

        EXAMPLES::

            sage: Cl = algebras.QuantumClifford(3)
            sage: Cl.dimension()
            512

            sage: Cl = algebras.QuantumClifford(4, 2)  # long time
            sage: Cl.dimension()  # long time
            65536
        """
        return ZZ(8*self._k) ** self._n

    @cached_method
    def algebra_generators(self):
        r"""
        Return the algebra generators of ``self``.

        EXAMPLES::

            sage: Cl = algebras.QuantumClifford(3)
            sage: Cl.algebra_generators()
            Finite family {'psi0': psi0, 'psi1': psi1, 'psi2': psi2,
                           'psid0': psid0, 'psid1': psid1, 'psid2': psid2,
                           'w0': w0, 'w1': w1, 'w2': w2}
        """
        one = (0,) * self._n  # one in the corresponding free abelian group
        zero = [0] * self._n
        d = {}
        for i in range(self._n):
            r = list(zero)  # Make a copy
            r[i] = 1
            d['psi%s' % i] = self.monomial((self._psi(r), one))
            r[i] = -1
            d['psid%s' % i] = self.monomial((self._psi(r), one))
        zero = self._psi(zero)
        for i in range(self._n):
            temp = list(zero)  # Make a copy
            temp[i] = 1
            d['w%s' % i] = self.monomial((zero, tuple(temp)))
        return Family(sorted(d), lambda i: d[i])

    @cached_method
    def gens(self):
        r"""
        Return the generators of ``self``.

        EXAMPLES::

            sage: Cl = algebras.QuantumClifford(3)
            sage: Cl.gens()
            (psi0, psi1, psi2, psid0, psid1, psid2, w0, w1, w2)
        """
        return tuple(self.algebra_generators())

    @cached_method
    def one_basis(self):
        r"""
        Return the index of the basis element of `1`.

        EXAMPLES::

            sage: Cl = algebras.QuantumClifford(3)
            sage: Cl.one_basis()
            ((0, 0, 0), (0, 0, 0))
        """
        return (self._psi([0]*self._n), (0,)*self._n)


class QuantumCliffordAlgebraGeneric(QuantumCliffordAlgebra):
    r"""
    The quantum Clifford algebra when `q^{2k} \neq 1`.

    The *quantum Clifford algebra*, or `q`-Clifford algebra,
    of rank `n` and twist `k` is the unital associative algebra
    `\mathrm{Cl}_{q}(n, k)` over a field `F` with generators
    `\psi_a, \psi_a^*, \omega_a` for `a = 1, \ldots, n` that
    satisfy the following relations:

    .. MATH::

        \begin{aligned}
        \omega_a \omega_b & = \omega_b \omega_a,
        & \omega_a^{4k} & = (1 + q^{-2k}) \omega_a^{2k} - q^{-2k},
        \\ \omega_a \psi_b & = q^{\delta_{ab}} \psi_b \omega_a,
        & \omega_a \psi^*_b & = \psi^*_b \omega_a,
        \\ \psi_a \psi_b & + \psi_b \psi_a = 0,
        & \psi^*_a \psi^*_b & + \psi^*_b \psi^*_a = 0,
        \\ \psi_a \psi^*_a & = \frac{q^k \omega_a^{3k} - q^{-k} \omega_a^k}{q^k - q^{-k}},
        & \psi^*_a \psi_a & = \frac{q^{2k} (\omega_a - \omega_a^{3k})}{q^k - q^{-k}},
        \\ \psi_a \psi^*_b & + \psi_b^* \psi_a = 0
        & & \text{if } a \neq b,
        \end{aligned}

    where `q \in F` such that `q^{2k} \neq 1`.

    When `k = 2`, we recover the original definition given by Hayashi in
    [Hayashi1990]_. The `k = 1` version was used in [Kwon2014]_.
    """
    def __init__(self, n, k, q, F):
        r"""
        Initialize ``self``.

        TESTS::

            sage: Cl = algebras.QuantumClifford(1,3)
            sage: TestSuite(Cl).run(elements=Cl.basis())  # long time

            sage: Cl = algebras.QuantumClifford(3)
            sage: elts = Cl.some_elements() + list(Cl.algebra_generators())
            sage: TestSuite(Cl).run(elements=elts)  # long time

            sage: Cl = algebras.QuantumClifford(2, 4)
            sage: elts = Cl.some_elements() + list(Cl.algebra_generators())
            sage: TestSuite(Cl).run(elements=elts)  # long time
        """
        psi = cartesian_product([(-1,0,1)]*n)
        indices = [(tuple(p), tuple(w))
                   for p in psi
                   for w in product(*[list(range((4-2*abs(p[i]))*k)) for i in range(n)])]
        super().__init__(n, k, q, F, psi, indices)

    def _repr_term(self, m):
        r"""
        Return a string representation of the basis element indexed by ``m``.

        EXAMPLES::

            sage: Cl = algebras.QuantumClifford(3, 3)
            sage: Cl._repr_term( ((1, 0, -1), (0, 2, 5)) )
            'psi0*psid2*w1^2*w2^5'
            sage: Cl._repr_term( ((1, 0, -1), (0, 0, 0)) )
            'psi0*psid2'
            sage: Cl._repr_term( ((0, 0, 0), (0, 2, 5)) )
            'w1^2*w2^5'
            sage: Cl._repr_term( ((0, 0, 0), (0, 0, 0)) )
            '1'

            sage: Cl(5)
            5
        """
        p, v = m
        rp = '*'.join('psi%s' % i if p[i] > 0 else 'psid%s' % i
                      for i in range(self._n) if p[i] != 0)
        gen_str = lambda e: '' if e == 1 else '^%s' % e
        rv = '*'.join('w%s' % i + gen_str(v[i]) for i in range(self._n) if v[i] != 0)
        if rp:
            if rv:
                return rp + '*' + rv
            return rp
        if rv:
            return rv
        return '1'

    def _latex_term(self, m):
        r"""
        Return a latex representation for the basis element indexed by ``m``.

        EXAMPLES::

            sage: Cl = algebras.QuantumClifford(3, 3)
            sage: Cl._latex_term( ((1, 0, -1), (0, -2, 5)) )
            '\\psi_{0}\\psi^{\\dagger}_{2}\\omega_{1}^{-2}\\omega_{2}^{5}'
            sage: Cl._latex_term( ((1, 0, -1), (0, 0, 0)) )
            '\\psi_{0}\\psi^{\\dagger}_{2}'
            sage: Cl._latex_term( ((0, 0, 0), (0, -2, 5)) )
            '\\omega_{1}^{-2}\\omega_{2}^{5}'
            sage: Cl._latex_term( ((0, 0, 0), (0, 0, 0)) )
            '1'

            sage: latex(Cl(5))
            5
        """
        p, v = m
        rp = ''.join('\\psi_{%s}' % i if p[i] > 0 else '\\psi^{\\dagger}_{%s}' % i
                     for i in range(self._n) if p[i] != 0)
        gen_str = lambda e: '' if e == 1 else '^{%s}' % e
        rv = ''.join('\\omega_{%s}' % i + gen_str(v[i])
                     for i in range(self._n) if v[i] != 0)
        if not rp and not rv:
            return '1'
        return rp + rv

    @cached_method
    def product_on_basis(self, m1, m2):
        r"""
        Return the product of the basis elements indexed by ``m1`` and ``m2``.

        EXAMPLES::

            sage: Cl = algebras.QuantumClifford(3)
            sage: Cl.inject_variables()
            Defining psi0, psi1, psi2, psid0, psid1, psid2, w0, w1, w2
            sage: psi0^2  # indirect doctest
            0
            sage: psid0^2
            0
            sage: w0 * psi0
            q*psi0*w0
            sage: w0 * psid0
            1/q*psid0*w0
            sage: w2 * w0
            w0*w2
            sage: w0^4
            -1/q^2 + ((q^2+1)/q^2)*w0^2
        """
        p1, w1 = m1
        p2, w2 = m2

        # Check for \psi_i^2 == 0 and for the dagger version
        if any(p1[i] != 0 and p1[i] == p2[i] for i in range(self._n)):
            return self.zero()

        # \psi_i is represented by a 1 in p1[i] and p2[i]
        # \psi_i^{\dagger} is represented by a -1 in p1[i] and p2[i]

        k = self._k
        q_power = 0
        sign = 1
        pairings = []
        supported = []
        p = [0] * self._n
        # Move w1 * p2 to q^q_power * p2 * w1
        # Also find pairs \psi_i \psi_i^{\dagger} (or vice versa)
        for i in range(self._n):
            if p2[i] != 0:
                supported.append(i)
                q_power += w1[i] * p2[i]
                if p1[i] != 0:
                    # We make pairings 1-based because we cannot distinguish 0 and -0
                    pairings.append((i+1) * p1[i])
            # we know p1[i] != p2[i] if nonzero, so their sum is -1, 0, 1
            p[i] = p1[i] + p2[i]

        supported.append(self._n-1) # To get between the last support and the end
        # Get the sign of moving \psi_i and \psi_i^{\dagger} into position
        for i in reversed(range(1, len(supported))):
            if i % 2 != 0:
                for j in reversed(range(supported[i-1]+1, supported[i]+1)):
                    if p1[j] != 0:
                        sign = (-1)**i * sign

        # We move the pairs \psi_i \psi_i^{\dagger} (or the reverse) to the
        #   end of the \psi part. This does not change the sign because they
        #   move in pairs.
        # We take the products.
        vp = self._w_poly.gens()
        poly = self._w_poly.one()
        q = self._q
        for i in pairings:
            if i < 0:
                i = -i - 1 # Go back to 0-based
                vpik = -q**(2*k) * vp[i]**(3*k) + (1 + q**(2*k)) * vp[i]**k
                poly *= -(vp[i]**k - vpik) / (q**k - q**(-k))
            else:
                i -= 1 # Go back to 0-based
                vpik = -q**(2*k) * vp[i]**(3*k) + (1 + q**(2*k)) * vp[i]**k
                poly *= (q**k * vp[i]**k - q**(-k) * vpik) / (q**k - q**(-k))

        v = list(w1)
        for i in range(self._n):
            v[i] += w2[i]

        # For all \psi_i v_i^k, convert this to either k = 0, 1
        #   and same for \psi_i^{\dagger}
        for i in range(self._n):
            if p[i] > 0 and v[i] != 0:
                q_power -= 2 * k * (v[i] // (2*k))
                v[i] = v[i] % (2*k)
            if p[i] < 0 and v[i] != 0:
                v[i] = v[i] % (2*k)

        poly *= self._w_poly.monomial(*v)
        poly = poly.reduce([vp[i]**(4*k) - (1 + q**(-2*k)) * vp[i]**(2*k) + q**(-2*k)
                            for i in range(self._n)])
        pdict = poly.monomial_coefficients()
        ret = {(self._psi(p), tuple(e)): pdict[e] * q**q_power * sign
               for e in pdict}

        return self._from_dict(ret)

    class Element(CombinatorialFreeModule.Element):
        def inverse(self):
            r"""
            Return the inverse if ``self`` is a basis element.

            EXAMPLES::

                sage: Cl = algebras.QuantumClifford(2)
                sage: Cl.inject_variables()
                Defining psi0, psi1, psid0, psid1, w0, w1
                sage: w0^-1
                (q^2+1)*w0 - q^2*w0^3
                sage: w0^-1 * w0
                1
                sage: w0^-2
                (q^2+1) - q^2*w0^2
                sage: w0^-2 * w0^2
                1
                sage: w0^-2 * w0 == w0^-1
                True
                sage: w = w0 * w1
                sage: w^-1
                (q^4+2*q^2+1)*w0*w1 + (-q^4-q^2)*w0*w1^3
                 + (-q^4-q^2)*w0^3*w1 + q^4*w0^3*w1^3
                sage: w^-1 * w
                1
                sage: w * w^-1
                1

                sage: (2*w0)^-1
                ((q^2+1)/2)*w0 - q^2/2*w0^3

                sage: (w0 + w1)^-1
                Traceback (most recent call last):
                ...
                ValueError: cannot invert self (= w1 + w0)
                sage: (psi0 * w0)^-1
                Traceback (most recent call last):
                ...
                ValueError: cannot invert self (= psi0*w0)

                sage: Cl = algebras.QuantumClifford(1, 2)
                sage: Cl.inject_variables()
                Defining psi0, psid0, w0
                sage: (psi0 + psid0).inverse()
                psid0*w0^2 + q^2*psi0*w0^2

                sage: Cl = algebras.QuantumClifford(2, 2)
                sage: Cl.inject_variables()
                Defining psi0, psi1, psid0, psid1, w0, w1
                sage: w0^-1
                (q^4+1)*w0^3 - q^4*w0^7
                sage: w0 * w0^-1
                1
            """
            if not self:
                raise ZeroDivisionError
            if len(self) != 1:
                return super().__invert__()
            Cl = self.parent()
            ((p, w), coeff), = list(self._monomial_coefficients.items())
            if any(p[i] != 0 for i in range(Cl._n)):
                return super().__invert__()
            poly = Cl._w_poly.monomial(*w)
            wp = Cl._w_poly.gens()
            q = Cl._q
            k = Cl._k
            poly = poly.subs({wi: -q**(2*k) * wi**(4*k-1) + (1 + q**(2*k)) * wi**(2*k-1)
                              for wi in wp})
            poly = poly.reduce([wi**(4*k) - (1 + q**(-2*k)) * wi**(2*k) + q**(-2*k)
                                for wi in wp])
            pdict = poly.monomial_coefficients()
            coeff = coeff.inverse_of_unit()
            ret = {(p, tuple(e)): coeff * c for e, c in pdict.items()}
            return Cl.element_class(Cl, ret)

        __invert__ = inverse


class QuantumCliffordAlgebraRootUnity(QuantumCliffordAlgebra):
    r"""
    The quantum Clifford algebra when `q^{2k} = 1`.

    The *quantum Clifford algebra*, or `q`-Clifford algebra,
    of rank `n` and twist `k` is the unital associative algebra
    `\mathrm{Cl}_{q}(n, k)` over a field `F` with generators
    `\psi_a, \psi_a^*, \omega_a` for `a = 1, \ldots, n` that
    satisfy the following relations:

    .. MATH::

        \begin{aligned}
        \omega_a \omega_b & = \omega_b \omega_a,
        & \omega_a^{2k} & = 1,
        \\ \omega_a \psi_b & = q^{\delta_{ab}} \psi_b \omega_a,
        & \omega_a \psi^*_b & = \psi^*_b \omega_a,
        \\ \psi_a \psi_b & + \psi_b \psi_a = 0,
        & \psi^*_a \psi^*_b & + \psi^*_b \psi^*_a = 0,
        \\ \psi_a \psi^*_a & + q^k \psi^*_a \psi_a = \omega_a^k
        & \psi_a \psi^*_b & + \psi_b^* \psi_a = 0 \quad (a \neq b),
        \end{aligned}

    where `q \in F` such that `q^{2k} = 1`. This has further relations of

    .. MATH::

        \begin{aligned}
        \psi^*_a \psi_a \psi^*_a & = \psi^*_a \omega_a^k,
        \\
        \psi_a \psi^*_a \psi_a & = q^k \psi_a \omega_a^k,
        \\
        (\psi_a \psi^*_a)^2 & = \psi_a \psi^*_a \omega_a^k.
        \end{aligned}
    """
    def __init__(self, n, k, q, F):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: Cl = algebras.QuantumClifford(1,2,-1)
            sage: TestSuite(Cl).run(elements=Cl.basis())

            sage: z = CyclotomicField(3).gen()
            sage: Cl = algebras.QuantumClifford(1,3,z)
            sage: TestSuite(Cl).run(elements=Cl.basis())

            sage: Cl = algebras.QuantumClifford(3,1,-1)
            sage: elts = Cl.some_elements() + list(Cl.algebra_generators())
            sage: TestSuite(Cl).run(elements=elts)  # long time

            sage: Cl = algebras.QuantumClifford(2,4,-1)
            sage: elts = Cl.some_elements() + list(Cl.algebra_generators())
            sage: TestSuite(Cl).run(elements=elts)  # long time
        """
        psi = cartesian_product([(-1,0,1,2)]*n)
        indices = [(tuple(p), tuple(w))
                   for p in psi
                   for w in product(list(range(2*k)), repeat=n)]
        super().__init__(n, k, q, F, psi, indices)

    def _repr_term(self, m):
        r"""
        Return a string representation of the basis element indexed by ``m``.

        EXAMPLES::

            sage: Cl = algebras.QuantumClifford(3, 3, -1)
            sage: Cl._repr_term( ((1, 0, -1), (0, 2, 5)) )
            'psi0*psid2*w1^2*w2^5'
            sage: Cl._repr_term( ((1, 0, 2), (0, 0, 0)) )
            'psi0*psi2*psid2'
            sage: Cl._repr_term( ((0, 0, 0), (0, 2, 5)) )
            'w1^2*w2^5'
            sage: Cl._repr_term( ((0, 0, 0), (0, 0, 0)) )
            '1'

            sage: Cl(5)
            5
        """
        p, v = m

        def ppr(i):
            val = p[i]
            if val == -1:
                return 'psid%s' % i
            elif val == 1:
                return 'psi%s' % i
            elif val == 2:
                return 'psi%s*psid%s' % (i,i)

        rp = '*'.join(ppr(i) for i in range(self._n) if p[i] != 0)
        gen_str = lambda e: '' if e == 1 else '^%s' % e
        rv = '*'.join('w%s' % i + gen_str(v[i]) for i in range(self._n) if v[i] != 0)
        if rp:
            if rv:
                return rp + '*' + rv
            return rp
        if rv:
            return rv
        return '1'

    def _latex_term(self, m):
        r"""
        Return a latex representation for the basis element indexed by ``m``.

        EXAMPLES::

            sage: Cl = algebras.QuantumClifford(3, 3, -1)
            sage: Cl._latex_term( ((1, 0, -1), (0, 2, 5)) )
            '\\psi_{0}\\psi^{\\dagger}_{2}\\omega_{1}^{2}\\omega_{2}^{5}'
            sage: Cl._latex_term( ((1, 0, 2), (0, 0, 0)) )
            '\\psi_{0}\\psi_{2}\\psi^{\\dagger}_{2}'
            sage: Cl._latex_term( ((0, 0, 0), (0, 2, 5)) )
            '\\omega_{1}^{2}\\omega_{2}^{5}'
            sage: Cl._latex_term( ((0, 0, 0), (0, 0, 0)) )
            '1'

            sage: latex(Cl(5))
            5
        """
        p, v = m

        def ppr(i):
            val = p[i]
            if val == -1:
                return '\\psi^{\\dagger}_{%s}' % i
            elif val == 1:
                return '\\psi_{%s}' % i
            elif val == 2:
                return '\\psi_{%s}\\psi^{\\dagger}_{%s}' % (i, i)

        rp = ''.join(ppr(i) for i in range(self._n) if p[i] != 0)
        gen_str = lambda e: '' if e == 1 else '^{%s}' % e
        rv = ''.join('\\omega_{%s}' % i + gen_str(v[i])
                     for i in range(self._n) if v[i] != 0)
        if not rp and not rv:
            return '1'
        return rp + rv

    @cached_method
    def product_on_basis(self, m1, m2):
        r"""
        Return the product of the basis elements indexed by ``m1`` and ``m2``.

        EXAMPLES::

            sage: z = CyclotomicField(3).gen()
            sage: Cl = algebras.QuantumClifford(3, 3, z)
            sage: Cl.inject_variables()
            Defining psi0, psi1, psi2, psid0, psid1, psid2, w0, w1, w2
            sage: psi0^2  # indirect doctest
            0
            sage: psid0^2
            0
            sage: w0 * psi0
            -(-zeta3)*psi0*w0
            sage: w0 * psid0
            -(zeta3+1)*psid0*w0
            sage: psi0 * psid0
            psi0*psid0
            sage: psid0 * psi0
            w0^3 - psi0*psid0
            sage: w2 * w0
            w0*w2
            sage: w0^6
            1
            sage: psi0 * psi1
            psi0*psi1
            sage: psi1 * psi0
            -psi0*psi1
            sage: psi1 * (psi0 * psi2)
            -psi0*psi1*psi2

            sage: z = CyclotomicField(6).gen()
            sage: Cl = algebras.QuantumClifford(3, 3, z)
            sage: Cl.inject_variables()
            Defining psi0, psi1, psi2, psid0, psid1, psid2, w0, w1, w2

            sage: psid1 * (psi1 * psid1)
            psid1*w1^3
            sage: (psi1* psid1) * (psi1 * psid1)
            psi1*psid1*w1^3
            sage: (psi1 * psid1) * psi1
            -psi1*w1^3
        """
        p1, w1 = m1
        p2, w2 = m2
        k = self._k
        tk = 2 * k

        # \psi_i is represented by a 1 in p1[i] and p2[i]
        # \psi_i^{\dagger} is represented by a -1 in p1[i] and p2[i]
        # \psi_i \psi_i^{\dagger} is a 2 in p1[i] and p2[i]

        # Check for \psi_i^2 == 0 and for the dagger version
        if any((p1[i] % 2 != 0 and p1[i] == p2[i])
               or (p1[i] == 2 and p2[i] == -1) or (p2[i] == 2 and p1[i] == 1)
               for i in range(self._n)):
            return self.zero()

        # Reduce any v_i^{2k} = 1
        v = [(w1[i] + w2[i]) % tk for i in range(self._n)]

        q_power = 0
        sign = 1
        pairings = []
        p = [0] * self._n
        # Move w1 * p2 to q^q_power * p2 * w1
        # Also find pairs \psi_i \psi_i^{\dagger} (or vice versa) and
        #   count the sign from moving \psi_i and \psi_i^{\dagger} into position
        num_cross = 0
        total_cross = 0
        for i in range(self._n):
            num_cross += p2[i]
            if p2[i] == 2:
                # By the above check, we cannot have p1[i] == 1
                if p1[i] != 0:
                    v[i] = (v[i] + k) % tk
                    p[i] = p1[i]
                else:
                    p[i] = p2[i]
            elif p2[i] != 0: # == +1, -1
                q_power += w1[i] * p2[i]
                # By the above check, we cannot have p1[i] == p2[i]
                if p1[i] == -1:
                    pairings.append(i)
                    total_cross -= 1 # correction
                    p[i] = None
                elif p1[i] == 1:
                    total_cross -= 1 # correction
                    p[i] = 2
                elif p1[i] == 2:
                    q_power += k
                    v[i] = (v[i] + k) % tk
                    p[i] = p2[i]
                else:
                    p[i] = p2[i]
            else:
                p[i] = p1[i]  # since p2[i] == 0

            if abs(p1[i]) == 1:
                total_cross += num_cross

        # total_cross does not need to actually be the total number of crossings; just correct mod 2
        if total_cross % 2:
            sign = -sign

        # Replace \psi_i^{\dagger} \psi_i = q^k ( w^k - \psi_i \psi_i^{\dagger} )
        def key(X):
            e = list(v)  # Make a copy
            for i in pairings:
                if i in X:
                    p[i] = 2
                else:
                    p[i] = 0
                    e[i] = (e[i] + k) % tk
            return (self._psi(p), tuple(e))

        q = self._q
        ret = {key(X): (-1)**len(X) * sign * q**(q_power+k*(len(pairings) % 2))
               for X in powerset(pairings)}

        return self._from_dict(ret)

    class Element(QuantumCliffordAlgebra.Element):
        def inverse(self):
            r"""
            Return the inverse if ``self`` is a basis element.

            EXAMPLES::

                sage: Cl = algebras.QuantumClifford(3, 3, -1)
                sage: Cl.inject_variables()
                Defining psi0, psi1, psi2, psid0, psid1, psid2, w0, w1, w2
                sage: w0^-1
                w0^5
                sage: w0^-1 * w0
                1
                sage: w0^-2
                w0^4
                sage: w0^-2 * w0^2
                1
                sage: w0^-2 * w0 == w0^-1
                True
                sage: w = w0 * w1^3
                sage: w^-1
                w0^5*w1^3
                sage: w^-1 * w
                1
                sage: w * w^-1
                1

                sage: (2*w0)^-1
                1/2*w0^5

                sage: Cl = algebras.QuantumClifford(3, 1, -1)
                sage: Cl.inject_variables()
                Defining psi0, psi1, psi2, psid0, psid1, psid2, w0, w1, w2

                sage: (w0 + w1)^-1
                Traceback (most recent call last):
                ...
                ValueError: cannot invert self (= w1 + w0)
                sage: (psi0 * w0)^-1
                Traceback (most recent call last):
                ...
                ValueError: cannot invert self (= psi0*w0)

                sage: z = CyclotomicField(6).gen()
                sage: Cl = algebras.QuantumClifford(1, 3, z)
                sage: Cl.inject_variables()
                Defining psi0, psid0, w0
                sage: (psi0 + psid0).inverse()
                psid0*w0^3 - psi0*w0^3

                sage: Cl = algebras.QuantumClifford(2, 2, -1)
                sage: Cl.inject_variables()
                Defining psi0, psi1, psid0, psid1, w0, w1
                sage: w0^-1
                w0^3
                sage: w0 * w0^-1
                1
            """
            if not self:
                raise ZeroDivisionError
            if len(self) != 1:
                return super().__invert__()
            Cl = self.parent()
            ((p, w), coeff), = list(self._monomial_coefficients.items())
            if any(p[i] != 0 for i in range(Cl._n)):
                return super().__invert__()
            tk = 2 * Cl._k
            w = tuple([tk - val if val else 0 for val in w])
            return Cl.element_class(Cl, {(p, w): coeff.inverse_of_unit()})

        __invert__ = inverse
