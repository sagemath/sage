# -*- coding: utf-8 -*-
r"""
`p`-adic `L`-functions of elliptic curves

To an elliptic curve `E` over the rational numbers and a prime `p`, one
can associate a `p`-adic L-function; at least if `E` does not have additive
reduction at `p`. This function is defined by interpolation of L-values of `E`
at twists. Through the main conjecture of Iwasawa theory it should also be
equal to a characteristic series of a certain Selmer group.

If `E` is ordinary, then it is an element of the Iwasawa algebra
`\Lambda(\ZZ_p^\times) = \ZZ_p[\Delta][\![T]\!]`, where `\Delta` is the group
of `(p-1)`-st roots of unity in `\ZZ_p^\times`, and `T = [\gamma] - 1` where
`\gamma = 1 + p` is a generator of `1 + p\ZZ_p`. (There is a slightly different
description for `p = 2`.)

One can decompose this algebra as the direct product of the subalgebras
corresponding to the characters of `\Delta`, which are simply the powers
`\tau^\eta` (`0 \le \eta \le p-2`) of the Teichmueller character `\tau: \Delta
\to \ZZ_p^\times`. Projecting the L-function into these components gives `p-1`
power series in `T`, each with coefficients in `\ZZ_p`.

If `E` is supersingular, the series will have coefficients in a quadratic
extension of `\QQ_p`, and the coefficients will be unbounded. In this case we
have only implemented the series for `\eta = 0`. We have also implemented the
`p`-adic L-series as formulated by Perrin-Riou [BP1993]_, which has coefficients in
the Dieudonné module `D_pE = H^1_{dR}(E/\QQ_p)` of `E`. There is a different
description by Pollack [Pol2003]_ which is not available here.

According to the `p`-adic version of the Birch and Swinnerton-Dyer conjecture
[MTT1986]_, the order of vanishing of the `L`-function at the trivial character
(i.e. of the series for `\eta = 0` at `T = 0`) is just the rank of `E(\QQ)`, or
this rank plus one if the reduction at `p` is split multiplicative.

See [SW2013]_ for more details.

AUTHORS:

- William Stein (2007-01-01): first version

- Chris Wuthrich (22/05/2007): changed minor issues and added supersingular things

- Chris Wuthrich (11/2008): added quadratic_twists

- David Loeffler (01/2011): added nontrivial Teichmueller components

"""

######################################################################
#       Copyright (C) 2007 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  https://www.gnu.org/licenses/
######################################################################

from sage.rings.integer_ring import   ZZ
from sage.rings.rational_field import QQ
from sage.rings.padics.factory import Qp
from sage.rings.infinity import infinity
from sage.rings.all import LaurentSeriesRing, PowerSeriesRing, PolynomialRing, Integers

from sage.rings.integer import Integer
from sage.arith.all import valuation, binomial, kronecker_symbol, gcd, prime_divisors

from sage.structure.sage_object import SageObject
from sage.structure.richcmp import richcmp_method, richcmp

from sage.misc.all import denominator
from sage.misc.verbose import verbose, get_verbose
import sage.arith.all as arith

from sage.modules.free_module_element import vector
import sage.matrix.all as matrix
import sage.schemes.hyperelliptic_curves.monsky_washnitzer
from sage.functions.log import log
from sage.functions.other import floor
from sage.misc.cachefunc import cached_method

@richcmp_method
class pAdicLseries(SageObject):
    r"""
    The `p`-adic L-series of an elliptic curve.

    EXAMPLES:

    An ordinary example::

        sage: e = EllipticCurve('389a')
        sage: L = e.padic_lseries(5)
        sage: L.series(0)
        Traceback (most recent call last):
        ...
        ValueError: n (=0) must be a positive integer
        sage: L.series(1)
        O(T^1)
        sage: L.series(2)
        O(5^4) + O(5)*T + (4 + O(5))*T^2 + (2 + O(5))*T^3 + (3 + O(5))*T^4 + O(T^5)
        sage: L.series(3, prec=10)
        O(5^5) + O(5^2)*T + (4 + 4*5 + O(5^2))*T^2 + (2 + 4*5 + O(5^2))*T^3 + (3 + O(5^2))*T^4 + (1 + O(5))*T^5 + O(5)*T^6 + (4 + O(5))*T^7 + (2 + O(5))*T^8 + O(5)*T^9 + O(T^10)
        sage: L.series(2,quadratic_twist=-3)
        2 + 4*5 + 4*5^2 + O(5^4) + O(5)*T + (1 + O(5))*T^2 + (4 + O(5))*T^3 + O(5)*T^4 + O(T^5)

    A prime p such that E[p] is reducible::

        sage: L = EllipticCurve('11a').padic_lseries(5)
        sage: L.series(1)
        5 + O(5^2) + O(T)
        sage: L.series(2)
        5 + 4*5^2 + O(5^3) + O(5^0)*T + O(5^0)*T^2 + O(5^0)*T^3 + O(5^0)*T^4 + O(T^5)
        sage: L.series(3)
        5 + 4*5^2 + 4*5^3 + O(5^4) + O(5)*T + O(5)*T^2 + O(5)*T^3 + O(5)*T^4 + O(T^5)

    An example showing the calculation of nontrivial Teichmueller twists::

        sage: E = EllipticCurve('11a1')
        sage: lp = E.padic_lseries(7)
        sage: lp.series(4,eta=1)
        3 + 7^3 + 6*7^4 + 3*7^5 + O(7^6) + (2*7 + 7^2 + O(7^3))*T + (1 + 5*7^2 + O(7^3))*T^2 + (4 + 4*7 + 4*7^2 + O(7^3))*T^3 + (4 + 3*7 + 7^2 + O(7^3))*T^4 + O(T^5)
        sage: lp.series(4,eta=2)
        5 + 6*7 + 4*7^2 + 2*7^3 + 3*7^4 + 2*7^5 + O(7^6) + (6 + 4*7 + 7^2 + O(7^3))*T + (3 + 2*7^2 + O(7^3))*T^2 + (1 + 4*7 + 7^2 + O(7^3))*T^3 + (6 + 6*7 + 6*7^2 + O(7^3))*T^4 + O(T^5)
        sage: lp.series(4,eta=3)
        O(7^6) + (5 + 4*7 + 2*7^2 + O(7^3))*T + (6 + 5*7 + 2*7^2 + O(7^3))*T^2 + (5*7 + O(7^3))*T^3 + (7 + 4*7^2 + O(7^3))*T^4 + O(T^5)

    (Note that the last series vanishes at `T = 0`, which is consistent with ::

        sage: E.quadratic_twist(-7).rank()
        1

    This proves that `E` has rank 1 over `\QQ(\zeta_7)`.)

    TESTS:

    The load-dumps test::

        sage: lp = EllipticCurve('11a').padic_lseries(5)
        sage: lp == loads(dumps(lp))
        True
    """
    def __init__(self, E, p, implementation = 'eclib', normalize='L_ratio'):
        r"""
        INPUT:

        -  ``E`` -- an elliptic curve
        -  ``p`` -- a prime of good reduction
        -  ``implementation`` -- string (default:'eclib'); either 'eclib' to use
           John Cremona's ``eclib`` for the computation of modular
           symbols, 'num' to use numerical modular symbols
           or 'sage' to use Sage's own implementation
        -  ``normalize`` -- ``'L_ratio'`` (default), ``'period'`` or ``'none'``;
           this is describes the way the modular symbols
           are normalized. See ``modular_symbol`` of
           an elliptic curve over Q for more details.

        EXAMPLES::

            sage: E = EllipticCurve('11a1')
            sage: Lp = E.padic_lseries(3)
            sage: Lp.series(2,prec=3)
            2 + 3 + 3^2 + 2*3^3 + O(3^4) + (1 + O(3))*T + (1 + O(3))*T^2 + O(T^3)
        """
        self._E = E
        self._p = ZZ(p)
        self._normalize = normalize
        if implementation not in ['eclib', 'sage', 'num']:
            raise ValueError("Implementation should be one of 'eclib', 'num' or 'sage'")
        self._implementation = implementation
        if not self._p.is_prime():
            raise ValueError("p (=%s) must be a prime" % p)
        if E.conductor() % (self._p)**2 == 0:
            raise NotImplementedError("p (=%s) must be a prime of semi-stable reduction" % p)

        try:
            E.label()
        except LookupError:
            if implementation != 'num':
                print("Warning : Curve outside Cremona's table. Computations of modular symbol space might take very long !")

        self._modular_symbol = E.modular_symbol(sign=+1,
                                                implementation=implementation,
                                                normalize=normalize)

    def __add_negative_space(self):
        r"""
        A helper function not designed for direct use.

        This function add the attribute ``_negative_modular_symbol`` to the class. This may take time
        and will only be needed when twisting with negative fundamental discriminants.

        EXAMPLES::

            sage: E = EllipticCurve('11a1')
            sage: lp = E.padic_lseries(5)
            sage: lp.modular_symbol(1/7,sign=-1)  #indirect doctest
            -1/2
        """
        self._negative_modular_symbol = self._E.modular_symbol(sign=-1, implementation="sage", normalize=self._normalize)

    def __richcmp__(self, other, op):
        r"""
        Compare ``self`` and ``other``.

        TESTS::

            sage: lp1 = EllipticCurve('11a1').padic_lseries(5)
            sage: lp2 = EllipticCurve('11a1').padic_lseries(7)
            sage: lp3 = EllipticCurve('11a2').padic_lseries(5)
            sage: lp1 == lp1
            True
            sage: lp1 == lp2
            False
            sage: lp1 == lp3
            False
        """
        if type(self) != type(other):
            return NotImplemented
        return richcmp((self._E, self._p), (other._E, other._p), op)

    def elliptic_curve(self):
        r"""
        Return the elliptic curve to which this `p`-adic L-series is associated.

        EXAMPLES::

            sage: L = EllipticCurve('11a').padic_lseries(5)
            sage: L.elliptic_curve()
            Elliptic Curve defined by y^2 + y = x^3 - x^2 - 10*x - 20 over Rational Field
        """
        return self._E

    def prime(self):
        r"""
        Return the prime `p` as in 'p-adic L-function'.

        EXAMPLES::

            sage: L = EllipticCurve('11a').padic_lseries(5)
            sage: L.prime()
            5
        """
        return self._p

    def _repr_(self):
        r"""
        Return print representation.

        EXAMPLES::

            sage: e = EllipticCurve('37a')
            sage: e.padic_lseries(3)._repr_()
            '3-adic L-series of Elliptic Curve defined by y^2 + y = x^3 - x over Rational Field'
            sage: e.padic_lseries(3,normalize='none')
            3-adic L-series of Elliptic Curve defined by y^2 + y = x^3 - x over Rational Field (not normalized)
            sage: L = e.padic_lseries(3,normalize='none')
            sage: L.rename('(factor)*L_3(T)')
            sage: L
            (factor)*L_3(T)
        """
        s = "%s-adic L-series of %s" % (self._p, self._E)
        if not self._normalize == 'L_ratio':
            s += ' (not normalized)'
        return s

    def modular_symbol(self, r, sign=+1, quadratic_twist=+1):
        r"""
        Return the modular symbol evaluated at `r`.

        This is used to compute this `p`-adic L-series.

        Note that the normalization is not correct at this
        stage: use ``_quotient_of periods_to_twist`` to correct.

        Note also that this function does not check if the condition
        on the quadratic_twist=D is satisfied. So the result will only
        be correct if for each prime `\ell` dividing `D`, we have
        `ord_{\ell}(N)<= ord_{\ell}(D)`, where `N` is the conductor of the curve.

        INPUT:

        -  ``r`` -- a cusp given as either a rational number or oo

        -  ``sign`` -- +1 (default) or -1 (only implemented without twists)

        -  ``quadratic_twist`` -- a fundamental discriminant of a quadratic field or +1 (default)

        EXAMPLES::

            sage: E = EllipticCurve('11a1')
            sage: lp = E.padic_lseries(5)
            sage: [lp.modular_symbol(r) for r in [0,1/5,oo,1/11]]
            [1/5, 6/5, 0, 0]
            sage: [lp.modular_symbol(r,sign=-1) for r in [0,1/3,oo,1/7]]
            [0, 1/2, 0, -1/2]
            sage: [lp.modular_symbol(r,quadratic_twist=-20) for r in [0,1/5,oo,1/11]]
            [1, 1, 0, 1/2]

            sage: E = EllipticCurve('20a1')
            sage: Et = E.quadratic_twist(-4)
            sage: lpt = Et.padic_lseries(5)
            sage: eta = lpt._quotient_of_periods_to_twist(-4)
            sage: lpt.modular_symbol(0) == lp.modular_symbol(0,quadratic_twist=-4) / eta
            True
        """
        if quadratic_twist == +1:
            if sign == +1:
                return self._modular_symbol(r)
            elif sign == -1:
                try:
                    m = self._negative_modular_symbol
                except (KeyError, AttributeError):
                    if not hasattr(self, '_modular_symbol_negative'):
                        self.__add_negative_space()
                        m = self._negative_modular_symbol
                return m(r)
        else:
            D = quadratic_twist
            if sign == -1:
                raise NotImplementedError("Quadratic twists for negative modular symbols are not yet implemented.")
            if D > 0:
                m = self._modular_symbol
                return sum([kronecker_symbol(D, u) * m(r + ZZ(u) / D)
                            for u in range(1, D)])

            else:
                try:
                    m = self._negative_modular_symbol
                except (KeyError, AttributeError):
                    if not hasattr(self, '_modular_symbol_negative'):
                        self.__add_negative_space()
                        m = self._negative_modular_symbol
                return -sum([kronecker_symbol(D, u) * m(r + ZZ(u) / D)
                             for u in range(1, -D)])

    def measure(self, a, n, prec, quadratic_twist=+1, sign = +1):
        r"""
        Return the measure on `\ZZ_p^{\times}` defined by

           `\mu_{E,\alpha}^+ ( a + p^n \ZZ_p  ) =
           \frac{1}{\alpha^n} \left [\frac{a}{p^n}\right]^{+} -
           \frac{1}{\alpha^{n+1}} \left[\frac{a}{p^{n-1}}\right]^{+}`

        where `[\cdot]^{+}` is the modular symbol. This is used to define
        this `p`-adic L-function (at least when the reduction is good).

        The optional argument ``sign`` allows the minus symbol `[\cdot]^{-}` to
        be substituted for the plus symbol.

        The optional argument ``quadratic_twist`` replaces `E` by the twist in
        the above formula, but the twisted modular symbol is computed using a
        sum over modular symbols of `E` rather than finding the modular symbols
        for the twist. Quadratic twists are only implemented if the sign is
        `+1`.

        Note that the normalization is not correct at this
        stage: use  ``_quotient_of periods`` and ``_quotient_of periods_to_twist``
        to correct.

        Note also that this function does not check if the condition
        on the ``quadratic_twist=D`` is satisfied. So the result will only
        be correct if for each prime `\ell` dividing `D`, we have
        `ord_{\ell}(N)<= ord_{\ell}(D)`, where `N` is the conductor of the curve.

        INPUT:

        -  ``a`` -- an integer

        -  ``n`` -- a non-negative integer

        -  ``prec`` -- an integer

        -  ``quadratic_twist`` (default = 1) -- a fundamental discriminant of a quadratic field,
           should be coprime to the conductor of `E`

        - ``sign`` (default = 1) -- an integer, which should be `\pm 1`.

        EXAMPLES::

            sage: E = EllipticCurve('37a')
            sage: L = E.padic_lseries(5)
            sage: L.measure(1,2, prec=9)
            2 + 3*5 + 4*5^3 + 2*5^4 + 3*5^5 + 3*5^6 + 4*5^7 + 4*5^8 + O(5^9)
            sage: L.measure(1,2, quadratic_twist=8,prec=15)
            O(5^15)
            sage: L.measure(1,2, quadratic_twist=-4,prec=15)
            4 + 4*5 + 4*5^2 + 3*5^3 + 2*5^4 + 5^5 + 3*5^6 + 5^8 + 2*5^9 + 3*5^12 + 2*5^13 + 4*5^14 + O(5^15)

            sage: E = EllipticCurve('11a1')
            sage: a = E.quadratic_twist(-3).padic_lseries(5).measure(1,2,prec=15)
            sage: b = E.padic_lseries(5).measure(1,2, quadratic_twist=-3,prec=15)
            sage: a == b * E.padic_lseries(5)._quotient_of_periods_to_twist(-3)
            True
        """
        s = ZZ(sign)
        if s not in [1, -1]:
            raise ValueError("Sign must be +- 1")
        if quadratic_twist != 1 and s != 1:
            raise NotImplementedError("Quadratic twists not implemented for sign -1")

        if quadratic_twist < 0:
            s = ZZ(-1)

        try:
            p, alpha, z, w, f = self.__measure_data[(n, prec, s)]
        except (KeyError, AttributeError):
            if not hasattr(self, '__measure_data'):
                self.__measure_data = {}
            p = self._p
            alpha = self.alpha(prec=prec)
            z = 1/(alpha**n)
            w = p**(n-1)
            if s == +1:
                f = self._modular_symbol
            else:
                try:
                    f = self._negative_modular_symbol
                except (KeyError, AttributeError):
                    if not hasattr(self, '_modular_symbol_negative'):
                        self.__add_negative_space()
                        f = self._negative_modular_symbol
            self.__measure_data[(n, prec, s)] = (p, alpha, z, w, f)

        if quadratic_twist == 1:
            if self._E.conductor() % p == 0:
                return z * f(a/(p*w))
            return z * ( f(a/(p*w)) - f(a/w) / alpha)
        else:
            D = quadratic_twist
            if self.is_ordinary():
                chip = kronecker_symbol(D,p)
            else:
                chip = 1 # alpha is +- sqrt(-p) anyway
            if self._E.conductor() % p == 0:
                mu = chip**n * z * sum([kronecker_symbol(D,u) * f(a/(p*w)+ZZ(u)/D) for u in range(1,D.abs())])
            else:
                mu = chip**n * z * sum([kronecker_symbol(D,u) *( f(a/(p*w)+ZZ(u)/D) - chip /alpha * f(a/w+ZZ(u)/D) ) for u in range(1,D.abs())])
            return s*mu

    def alpha(self, prec=20):
        r"""
        Return a `p`-adic root `\alpha` of the polynomial `x^2 - a_p x
        + p` with `ord_p(\alpha) < 1`.  In the ordinary case this is
        just the unit root.

        INPUT:

        - ``prec`` -- positive integer, the `p`-adic precision of the root.

        EXAMPLES:

        Consider the elliptic curve 37a::

            sage: E = EllipticCurve('37a')

        An ordinary prime::

            sage: L = E.padic_lseries(5)
            sage: alpha = L.alpha(10); alpha
            3 + 2*5 + 4*5^2 + 2*5^3 + 5^4 + 4*5^5 + 2*5^7 + 5^8 + 5^9 + O(5^10)
            sage: alpha^2 - E.ap(5)*alpha + 5
            O(5^10)

        A supersingular prime::

            sage: L = E.padic_lseries(3)
            sage: alpha = L.alpha(10); alpha
            alpha + O(alpha^21)
            sage: alpha^2 - E.ap(3)*alpha + 3
            O(alpha^22)

        A reducible prime::

            sage: L = EllipticCurve('11a').padic_lseries(5)
            sage: L.alpha(5)
            1 + 4*5 + 3*5^2 + 2*5^3 + 4*5^4 + O(5^5)
        """
        try:
            return self._alpha[prec]
        except AttributeError:
            self._alpha = {}
        except KeyError:
            pass
        E = self._E
        p = self._p
        a_p = E.ap(p)
        K = Qp(p, prec, print_mode='series')

        if E.conductor() % p == 0:
            self._alpha[prec] = K(a_p)
            return K(a_p)

        R = ZZ['x']
        f = R([p, -a_p, 1])
        if E.is_ordinary(p):
            G = f.factor_padic(p, prec + 5)
            for pr, e in G:
                a = -pr[0]
                if a.valuation() < 1:
                    self._alpha[prec] = K(a)
                    return K(a)
            raise RuntimeError("bug in p-adic L-function alpha")
        else: # supersingular case
            f = f.change_ring(K)
            A = K.extension(f, names="alpha")
            a = A.gen()
            self._alpha[prec] = a
            return a

    def order_of_vanishing(self):
        r"""
        Return the order of vanishing of this `p`-adic L-series.

        The output of this function is provably correct, due to a
        theorem of Kato [Kat2004]_.

        .. NOTE:: currently `p` must be a prime of good ordinary reduction.

        REFERENCES:

        - [MTT1986]_

        - [Kat2004]_

        EXAMPLES::

            sage: L = EllipticCurve('11a').padic_lseries(3)
            sage: L.order_of_vanishing()
            0
            sage: L = EllipticCurve('11a').padic_lseries(5)
            sage: L.order_of_vanishing()
            0
            sage: L = EllipticCurve('37a').padic_lseries(5)
            sage: L.order_of_vanishing()
            1
            sage: L = EllipticCurve('43a').padic_lseries(3)
            sage: L.order_of_vanishing()
            1
            sage: L = EllipticCurve('37b').padic_lseries(3)
            sage: L.order_of_vanishing()
            0
            sage: L = EllipticCurve('389a').padic_lseries(3)
            sage: L.order_of_vanishing()
            2
            sage: L = EllipticCurve('389a').padic_lseries(5)
            sage: L.order_of_vanishing()
            2
            sage: L = EllipticCurve('5077a').padic_lseries(5, implementation = 'eclib')
            sage: L.order_of_vanishing()
            3
        """
        try:
            return self.__ord
        except AttributeError:
            pass

        if not self.is_ordinary():
            raise NotImplementedError
        E = self.elliptic_curve()
        if not E.is_good(self.prime()):
            raise ValueError("prime must be of good reduction")
        r = E.rank()
        n = 1
        while True:
            f = self.series(n)
            v = f.valuation()
            if v < n and v < r:
                raise RuntimeError("while computing p-adic order of vanishing, got a contradiction: the curve is %s, the curve has rank %s, but the p-adic L-series vanishes to order <= %s" % (E, r, v))
            if v == r:
                self.__ord = v
                return v
            n += 1

    def teichmuller(self, prec):
        r"""
        Return Teichmuller lifts to the given precision.

        INPUT:

        - ``prec`` - a positive integer.

        OUTPUT:

        - a list of `p`-adic numbers, the cached Teichmuller lifts

        EXAMPLES::

            sage: L = EllipticCurve('11a').padic_lseries(7)
            sage: L.teichmuller(1)
            [0, 1, 2, 3, 4, 5, 6]
            sage: L.teichmuller(2)
            [0, 1, 30, 31, 18, 19, 48]
        """
        p = self._p
        K = Qp(p, prec, print_mode='series')
        return [Integer(0)] + \
               [a.residue(prec).lift() for a in K.teichmuller_system()]

    def _e_bounds(self, n, prec):
        r"""
        A helper function not designed for direct use.

        It computes the valuations of the coefficients of `\omega_n = (1+T)^{p^n}-1`.

        EXAMPLES::

            sage: E = EllipticCurve('11a1')
            sage: Lp = E.padic_lseries(2)
            sage: Lp._e_bounds(1,10)
            [+Infinity, 1, 0, 0, 0, 0, 0, 0, 0, 0]
            sage: Lp._e_bounds(2,10)
            [+Infinity, 2, 1, 1, 0, 0, 0, 0, 0, 0]
            sage: Lp._e_bounds(3,10)
            [+Infinity, 3, 2, 2, 1, 1, 1, 1, 0, 0]
            sage: Lp._e_bounds(4,10)
            [+Infinity, 4, 3, 3, 2, 2, 2, 2, 1, 1]
        """
        # trac 10280: replace with new corrected code, note that the sequence has to be decreasing.
        pn = self._p**n
        enj = infinity
        res = [enj]
        for j in range(1,prec):
            bino = valuation(binomial(pn,j),self._p)
            if bino < enj:
                enj = bino
            res.append(enj)
        return res

    def _get_series_from_cache(self, n, prec, D, eta):
        r"""
        A helper function not designed for direct use.

        It picks up the series in the cache if it has been previously computed.

        EXAMPLES::

            sage: E = EllipticCurve('11a1')
            sage: Lp = E.padic_lseries(5)
            sage: Lp._pAdicLseries__series = {}  # clear cached series
            sage: Lp._get_series_from_cache(3,5,1,0)
            sage: Lp.series(3,prec=5)
            5 + 4*5^2 + 4*5^3 + O(5^4) + O(5)*T + O(5)*T^2 + O(5)*T^3 + O(5)*T^4 + O(T^5)
            sage: Lp._get_series_from_cache(3,5,1,0)
            5 + 4*5^2 + 4*5^3 + O(5^4) + O(5)*T + O(5)*T^2 + O(5)*T^3 + O(5)*T^4 + O(T^5)
        """
        try:
            return self.__series[(n,prec,D,eta)]
        except AttributeError:
            self.__series = {}
        except KeyError:
            for _n, _prec, _D, _eta in self.__series:
                if _n == n and _D == D and _eta == eta and _prec >= prec:
                    return self.__series[(_n,_prec,_D,_eta)].add_bigoh(prec)
        return None

    def _set_series_in_cache(self, n, prec, D, eta, f):
        r"""
        A helper function not designed for direct use.

        It picks up the series in the cache if it has been previously computed.

        EXAMPLES::

            sage: E = EllipticCurve('11a1')
            sage: Lp = E.padic_lseries(5)
            sage: Lp.series(3,prec=5)
            5 + 4*5^2 + 4*5^3 + O(5^4) + O(5)*T + O(5)*T^2 + O(5)*T^3 + O(5)*T^4 + O(T^5)
            sage: Lp._set_series_in_cache(3,5,1,0,0)
            sage: Lp.series(3,prec=5)
            0
        """
        self.__series[(n, prec, D, eta)] = f

    def _quotient_of_periods_to_twist(self, D):
        r"""
        For a fundamental discriminant `D` of a quadratic number field this
        computes the constant `\eta` such that
        `\sqrt{\vert D\vert }\cdot\Omega_{E_D}^{+} =\eta\cdot \Omega_E^{sign(D)}`.
        As in [MTT1986]_ page 40. This is either 1 or 2 unless the condition
        on the twist is not satisfied, e.g. if we are 'twisting back' to a
        semi-stable curve.

        .. NOTE::

            No check on precision is made, so this may fail for huge `D`.

        EXAMPLES::

            sage: E = EllipticCurve('37b1')
            sage: lp = E.padic_lseries(3)
            sage: lp._quotient_of_periods_to_twist(-20)
            1
            sage: lp._quotient_of_periods_to_twist(-4)
            1
            sage: lp._quotient_of_periods_to_twist(-3)
            1
            sage: lp._quotient_of_periods_to_twist(-8)
            2
            sage: lp._quotient_of_periods_to_twist(8)
            2
            sage: lp._quotient_of_periods_to_twist(5)
            1
            sage: lp._quotient_of_periods_to_twist(12)
            1

            sage: E = EllipticCurve('11a1')
            sage: Et = E.quadratic_twist(-3)
            sage: lpt = Et.padic_lseries(5)
            sage: lpt._quotient_of_periods_to_twist(-3)
            3
        """
        # This function does not depend on p and could be moved out of
        # this file but it is needed only here

        # Note that the number of real components does not change by twisting.
        if D == 1:
            return 1
        Et = self._E.quadratic_twist(D)
        if D > 1:
            qt = Et.period_lattice().basis()[0] / self._E.period_lattice().basis()[0]
            qt *= qt.parent()(D).sqrt()
        else:
            qt = Et.period_lattice().basis()[1].imag() / self._E.period_lattice().basis()[0]
            if Et.real_components() == 1:
                qt *= 2
            qt *= qt.parent()(-D).sqrt()
        verbose('the real approximation is %s' % qt)
        # we know from MTT that the result has a denominator 1
        return QQ((8 * qt).round()) / 8


class pAdicLseriesOrdinary(pAdicLseries):
    def series(self, n=2, quadratic_twist=+1, prec=5, eta=0):
        r"""
        Return the `n`-th approximation to the `p`-adic L-series, in the
        component corresponding to the `\eta`-th power of the Teichmueller
        character, as a power series in `T` (corresponding to `\gamma-1` with
        `\gamma=1+p` as a generator of `1+p\ZZ_p`). Each coefficient is a
        `p`-adic number whose precision is provably correct.

        Here the normalization of the `p`-adic L-series is chosen
        such that `L_p(E,1) = (1-1/\alpha)^2 L(E,1)/\Omega_E`
        where `\alpha` is the unit root of the characteristic
        polynomial of Frobenius on `T_pE` and `\Omega_E` is the
        Néron period of `E`.

        INPUT:

        -  ``n`` - (default: 2) a positive integer
        -  ``quadratic_twist`` - (default: +1) a fundamental discriminant of a
           quadratic field, coprime to the conductor of the curve
        -  ``prec`` - (default: 5) maximal number of terms of the series to
           compute; to compute as many as possible just give a very large
           number for ``prec``; the result will still be correct.
        - ``eta`` (default: 0) an integer (specifying the power of the
          Teichmueller character on the group of roots of unity in
          `\ZZ_p^\times`)

        :meth:`power_series` is identical to ``series``.

        EXAMPLES:

        We compute some `p`-adic L-functions associated to the elliptic
        curve 11a::

            sage: E = EllipticCurve('11a')
            sage: p = 3
            sage: E.is_ordinary(p)
            True
            sage: L = E.padic_lseries(p)
            sage: L.series(3)
            2 + 3 + 3^2 + 2*3^3 + O(3^5) + (1 + 3 + O(3^2))*T + (1 + 2*3 + O(3^2))*T^2 + O(3)*T^3 + O(3)*T^4 + O(T^5)

        Another example at a prime of bad reduction, where the
        `p`-adic L-function has an extra 0 (compared to the non
        `p`-adic L-function)::

            sage: E = EllipticCurve('11a')
            sage: p = 11
            sage: E.is_ordinary(p)
            True
            sage: L = E.padic_lseries(p)
            sage: L.series(2)
            O(11^4) + (10 + O(11))*T + (6 + O(11))*T^2 + (2 + O(11))*T^3 + (5 + O(11))*T^4 + O(T^5)

        We compute a `p`-adic L-function that vanishes to order 2::

            sage: E = EllipticCurve('389a')
            sage: p = 3
            sage: E.is_ordinary(p)
            True
            sage: L = E.padic_lseries(p)
            sage: L.series(1)
            O(T^1)
            sage: L.series(2)
            O(3^4) + O(3)*T + (2 + O(3))*T^2 + O(T^3)
            sage: L.series(3)
            O(3^5) + O(3^2)*T + (2 + 2*3 + O(3^2))*T^2 + (2 + O(3))*T^3 + (1 + O(3))*T^4 + O(T^5)

        Checks if the precision can be changed (:trac:`5846`)::

            sage: L.series(3,prec=4)
            O(3^5) + O(3^2)*T + (2 + 2*3 + O(3^2))*T^2 + (2 + O(3))*T^3 + O(T^4)
            sage: L.series(3,prec=6)
            O(3^5) + O(3^2)*T + (2 + 2*3 + O(3^2))*T^2 + (2 + O(3))*T^3 + (1 + O(3))*T^4 + (1 + O(3))*T^5 + O(T^6)

        Rather than computing the `p`-adic L-function for the curve '15523a1', one can
        compute it as a quadratic_twist::

            sage: E = EllipticCurve('43a1')
            sage: lp = E.padic_lseries(3)
            sage: lp.series(2,quadratic_twist=-19)
            2 + 2*3 + 2*3^2 + O(3^4) + (1 + O(3))*T + (1 + O(3))*T^2 + O(T^3)
            sage: E.quadratic_twist(-19).label()    # optional -- database_cremona_ellcurve
            '15523a1'

        This proves that the rank of '15523a1' is zero, even if ``mwrank`` cannot determine this.

        We calculate the `L`-series in the nontrivial Teichmueller components::

            sage: L = EllipticCurve('110a1').padic_lseries(5, implementation="sage")
            sage: for j in [0..3]: print(L.series(4, eta=j))
            O(5^6) + (2 + 2*5 + 2*5^2 + O(5^3))*T + (5 + 5^2 + O(5^3))*T^2 + (4 + 4*5 + 2*5^2 + O(5^3))*T^3 + (1 + 5 + 3*5^2 + O(5^3))*T^4 + O(T^5)
            4 + 3*5 + 2*5^2 + 3*5^3 + 5^4 + O(5^6) + (1 + 3*5 + 4*5^2 + O(5^3))*T + (3 + 4*5 + 3*5^2 + O(5^3))*T^2 + (3 + 3*5^2 + O(5^3))*T^3 + (1 + 2*5 + 2*5^2 + O(5^3))*T^4 + O(T^5)
            2 + O(5^6) + (1 + 5 + O(5^3))*T + (2 + 4*5 + 3*5^2 + O(5^3))*T^2 + (4 + 5 + 2*5^2 + O(5^3))*T^3 + (4 + O(5^3))*T^4 + O(T^5)
            3 + 5 + 2*5^2 + 5^3 + 3*5^4 + 4*5^5 + O(5^6) + (1 + 2*5 + 4*5^2 + O(5^3))*T + (1 + 4*5 + O(5^3))*T^2 + (3 + 2*5 + 2*5^2 + O(5^3))*T^3 + (5 + 5^2 + O(5^3))*T^4 + O(T^5)

        It should now also work with `p=2` (:trac:`20798`)::

            sage: E = EllipticCurve("53a1")
            sage: lp = E.padic_lseries(2)
            sage: lp.series(7)
            O(2^8) + (1 + 2^2 + 2^3 + O(2^5))*T + (1 + 2^3 + O(2^4))*T^2 + (2^2 + 2^3 + O(2^4))*T^3 + (2 + 2^2 + O(2^3))*T^4 + O(T^5)

            sage: E = EllipticCurve("109a1")
            sage: lp = E.padic_lseries(2)
            sage: lp.series(6)
            2^2 + 2^6 + O(2^7) + (2 + O(2^4))*T + O(2^3)*T^2 + (2^2 + O(2^3))*T^3 + (2 + O(2^2))*T^4 + O(T^5)

        Check that twists by odd Teichmuller characters are ok (:trac:`32258`)::

            sage: E = EllipticCurve("443c1")
            sage: lp = E.padic_lseries(17, implementation="num")
            sage: l8 = lp.series(2,eta=8,prec=3)
            sage: l8.list()[0] - 1/lp.alpha()
            O(17^4)
            sage: lp = E.padic_lseries(2, implementation="num")
            sage: l1 = lp.series(8,eta=1,prec=3)
            sage: l1.list()[0] - 4/lp.alpha()^2
            O(2^9)


        """
        n = ZZ(n)
        if n < 1:
            raise ValueError("n (=%s) must be a positive integer" % n)
        if self._p == 2 and n == 1:
            raise ValueError("n (=%s) must be a at least 2 if p is 2" % n)
        if prec < 1:
            raise ValueError("Insufficient precision (%s)" % prec)

        # check if the conditions on quadratic_twist are satisfied
        eta = ZZ(eta) % (self._p- 1) if self._p != 2 else ZZ(eta) % 2
        D = ZZ(quadratic_twist)
        if D != 1:
            if eta != 0:
                raise NotImplementedError("quadratic twists only implemented for the 0th Teichmueller component")
            if D % 4 == 0:
                d = D//4
                if not d.is_squarefree() or d % 4 == 1:
                    raise ValueError("quadratic_twist (=%s) must be a fundamental discriminant of a quadratic field"%D)
            else:
                if not D.is_squarefree() or D % 4 != 1:
                    raise ValueError("quadratic_twist (=%s) must be a fundamental discriminant of a quadratic field"%D)
            if gcd(D,self._p) != 1:
                raise ValueError("quadratic twist (=%s) must be coprime to p (=%s) "%(D,self._p))
            if gcd(D, self._E.conductor()) != 1:
                for ell in prime_divisors(D):
                    if valuation(self._E.conductor(), ell) > valuation(D, ell):
                        raise ValueError("cannot twist a curve of conductor (=%s) by the quadratic twist (=%s)."%(self._E.conductor(),D))
        p = self._p
        si = 1-2*(eta % 2)

        #verbose("computing L-series for p=%s, n=%s, and prec=%s"%(p,n,prec))

        if prec == 1:
            if eta == 0:
                # trac 15737: if we only ask for the leading term we don't
                # need to do any sum as L_p(E,0) = (1-1/alpha)^2 * m(0) (good case)
                # set prec arbitrary to 20.
                K = Qp(p, 20, print_mode='series')
                R = PowerSeriesRing(K,'T',1)
                L = self.modular_symbol(0, sign=+1, quadratic_twist=D)
                chip = kronecker_symbol(D,p)
                if self._E.conductor() % p == 0:
                    L *= 1 - chip/self.alpha()
                else:
                    L *= (1-chip/self.alpha())**2
                L /= self._quotient_of_periods_to_twist(D)*self._E.real_components()
                L = R(L, 1)
                return L
            else:
                # here we need some sums anyway
                bounds = self._prec_bounds(n,prec,sign=si)
                padic_prec = 20
        else:
            bounds = self._prec_bounds(n,prec,sign=si)
            padic_prec = max(bounds[1:]) + 5

        verbose("using p-adic precision of %s"%padic_prec)

        if p == 2:
            res_series_prec = min(p**(n-2), prec)
        else:
            res_series_prec = min(p**(n-1), prec)
        verbose("using series precision of %s"%res_series_prec)

        ans = self._get_series_from_cache(n, res_series_prec,D,eta)
        if ans is not None:
            verbose("found series in cache")
            return ans

        K = QQ
        R = PowerSeriesRing(K,'T',res_series_prec)
        T = R(R.gen(),res_series_prec )
        L = R(0)
        one_plus_T_factor = R(1)
        gamma_power = K(1)
        teich = self.teichmuller(padic_prec)
        if p == 2:
            teich = [0, 1, -1]
            gamma = K(5)
            p_power = 2**(n-2)
            a_range = 3
        else:
            teich = self.teichmuller(padic_prec)
            gamma = K(1+ p)
            p_power = p**(n-1)
            a_range = p

        verbose("Now iterating over %s summands"%((p-1)*p_power))
        verbose_level = get_verbose()
        count_verb = 0
        for j in range(p_power):
            s = K(0)
            if verbose_level >= 2 and j/p_power*100 > count_verb + 3:
                verbose("%.2f percent done"%(float(j)/p_power*100))
                count_verb += 3
            for a in range(1,a_range):
                b = teich[a] * gamma_power
                s += teich[a]**eta * self.measure(b, n, padic_prec, quadratic_twist=D, sign=si).lift()
            L += s * one_plus_T_factor
            one_plus_T_factor *= 1+T
            gamma_power *= gamma

        verbose("the series before adjusting the precision is %s"%L)
        # Now create series but with each coefficient truncated
        # so it is proven correct:
        K = Qp(p, padic_prec, print_mode='series')
        R = PowerSeriesRing(K,'T',res_series_prec)
        L = R(L,res_series_prec)
        aj = L.list()
        if aj:
            aj = [aj[0].add_bigoh(padic_prec-2)] + \
                 [aj[j].add_bigoh(bounds[j]) for j in range(1,len(aj))]
        L = R(aj,res_series_prec )

        L /= self._quotient_of_periods_to_twist(D)
        if si == +1:
            L /= self._E.real_components()

        self._set_series_in_cache(n, res_series_prec, D, eta, L)

        return L

    power_series = series

    def is_ordinary(self):
        r"""
        Return ``True`` if the elliptic curve that this L-function is attached
        to is ordinary.

        EXAMPLES::

            sage: L = EllipticCurve('11a').padic_lseries(5)
            sage: L.is_ordinary()
            True
        """
        return True

    def is_supersingular(self):
        r"""
        Return ``True`` if the elliptic curve that this L function is attached
        to is supersingular.

        EXAMPLES::

            sage: L = EllipticCurve('11a').padic_lseries(5)
            sage: L.is_supersingular()
            False
        """
        return False

    @cached_method
    def _c_bound(self, sign=+1):
        r"""
        A helper function not designed for direct use.

        It returns an upper bound to the maximal `p`-adic valuation
        of the possible denominators  of the modular symbols appearing
        in the sum for the `p`-adic `L`-function with the given ``sign``.

        If the implementation of modular symbols used is 'sage', this is
        simply the maximum over all modular symbols. For others,
        we rely on the fact that the `p`-adic `L`-function is a sum of
        unitary modular symbols. These cusps are defined over `\QQ` and
        we know only need to find a torsion points on the `X_0`-optimal
        curve and compare the periods.

        EXAMPLES::

            sage: E = EllipticCurve('11a1')
            sage: Lp = E.padic_lseries(5)
            sage: Lp._c_bound()
            1
            sage: EllipticCurve('11a2').padic_lseries(5)._c_bound()
            0
            sage: EllipticCurve('11a3').padic_lseries(5)._c_bound()
            2
            sage: EllipticCurve('11a3').padic_lseries(5, implementation="sage")._c_bound()
            2
            sage: EllipticCurve('50b1').padic_lseries(3)._c_bound()
            0
            sage: EllipticCurve('50b1').padic_lseries(3, implementation="sage")._c_bound()
            1
            sage: l = EllipticCurve("11a1").padic_lseries(5)
            sage: ls = l.series(1,eta=1);
            sage: l._c_bound(sign=-1)
            0
        """
        E = self._E
        p = self._p
        N = self._E.conductor()
        if E.galois_representation().is_irreducible(p):
            return 0

        if self._implementation=="sage":
            m = E.modular_symbol_space(sign=1)
            b = m.boundary_map().codomain()
            C = b._known_cusps()  # all known, since computed the boundary map
            if sign == +1:
                return max([valuation(self.modular_symbol(a).denominator(), p) for a in C])
            else:
                try:
                    m = self._negative_modular_symbol
                except (KeyError, AttributeError):
                    if not hasattr(self, '_modular_symbol_negative'):
                        self._add_negative_space()
                        m = self._negative_modular_symbol
                return max([valuation(m(a).denominator(), p) for a in C])

        # else the same reasoning as in _set_denom in numerical
        # modular symbol. We rely on the fact that p is semistable
        from sage.databases.cremona import CremonaDatabase
        isog =  E.isogeny_class()
        t = 0
        if N <= CremonaDatabase().largest_conductor():
            E0 = E.optimal_curve()
        else:
            # we can't know which is the X_0-optimal curve
            # so we take one of the worst cases
            # if p=2 this may not be unique so we are cautious.
            ff = lambda C: C.period_lattice().complex_area()
            E0 = min(isog.curves, key=ff)
            if p == 2:
                t = 1
        # all modular symbols evaluated in a p-adic L-series
        # have denominator a power of p. Hence they come from
        # unitary cusps if p is semistable. Unitary cusps
        # are defined over Q, so they map to rational
        # torsion points on the X_0-optimal curve.
        if sign == 1:
            t += E.torsion_order().valuation(p)
        else:
            # no torsion point other than 2-torsion
            # can be non-real in the lattice
            if p == 2:
                t += 1
        if p == 2 and E0.real_components() == 1:
                t += 1 # slanted lattice

        # this was the bound for E0 now compare periods
        # to get the bound for E
        L0 = E0.period_lattice().basis()
        L = E.period_lattice().basis()
        if sign == 1:
            om = L[0]
            om0 = L0[0]
        else:
            om = L[1].imag()
            if E.real_components() == 1:
                om *= 2
            om0 = L[1].imag()
            if E0.real_components() == 1:
                om0 *= 2
        m = max(isog.matrix().list())
        q = (om/om0 *m).round()/m
        t += valuation(q,p)
        return max(t,0)

    def _prec_bounds(self, n, prec, sign=+1):
        r"""
        A helper function not designed for direct use.

        It returns the `p`-adic precisions of the approximation
        to the `p`-adic L-function.

        EXAMPLES::

            sage: E = EllipticCurve('11a1')
            sage: Lp = E.padic_lseries(5)
            sage: Lp._prec_bounds(3,10)
            [+Infinity, 1, 1, 1, 1, 0, 0, 0, 0, 0]
            sage: Lp._prec_bounds(3,12)
            [+Infinity, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0]
            sage: Lp._prec_bounds(4,5)
            [+Infinity, 2, 2, 2, 2]
            sage: Lp._prec_bounds(15,10)
            [+Infinity, 13, 13, 13, 13, 12, 12, 12, 12, 12]

            sage: Lp = E.padic_lseries(3)
            sage: Lp._prec_bounds(15,10)
            [+Infinity, 14, 14, 13, 13, 13, 13, 13, 13, 12]
        """
        if self._p == 2:
            e = self._e_bounds(n - 2, prec)
        else:
            e = self._e_bounds(n - 1, prec)
        c = self._c_bound()
        return [e[j] - c for j in range(len(e))]


class pAdicLseriesSupersingular(pAdicLseries):
    def series(self, n=3, quadratic_twist=+1, prec=5, eta=0):
        r"""
        Return the `n`-th approximation to the `p`-adic L-series as a
        power series in `T` (corresponding to `\gamma-1` with
        `\gamma=1+p` as a generator of `1+p\ZZ_p`).  Each
        coefficient is an element of a quadratic extension of the `p`-adic
        number whose precision is provably correct.

        Here the normalization of the `p`-adic L-series is chosen
        such that `L_p(E,1) = (1-1/\alpha)^2 L(E,1)/\Omega_E`
        where `\alpha` is a root of the characteristic
        polynomial of Frobenius on `T_pE` and `\Omega_E` is the
        Néron period of `E`.

        INPUT:

        -  ``n`` - (default: 2) a positive integer
        -  ``quadratic_twist`` - (default: +1) a fundamental discriminant of a
           quadratic field, coprime to the conductor of the curve
        -  ``prec`` - (default: 5) maximal number of terms of the series to
           compute; to compute as many as possible just give a very large
           number for ``prec``; the result will still be correct.
        - ``eta`` (default: 0) an integer (specifying the power of the
          Teichmueller character on the group of roots of unity in
          `\ZZ_p^\times`)

        OUTPUT:

        a power series with coefficients in a quadratic ramified extension of
        the `p`-adic numbers generated by a root `alpha` of the characteristic
        polynomial of Frobenius on `T_pE`.

        ALIAS: power_series is identical to series.

        EXAMPLES:

        A supersingular example, where we must compute to higher precision to see anything::

            sage: e = EllipticCurve('37a')
            sage: L = e.padic_lseries(3); L
            3-adic L-series of Elliptic Curve defined by y^2 + y = x^3 - x over Rational Field
            sage: L.series(2)
            O(T^3)
            sage: L.series(4)         # takes a long time (several seconds)
            O(alpha) + (alpha^-2 + O(alpha^0))*T + (alpha^-2 + O(alpha^0))*T^2 + O(T^5)
            sage: L.alpha(2).parent()
            3-adic Eisenstein Extension Field in alpha defined by x^2 + 3*x + 3

        An example where we only compute the leading term (:trac:`15737`)::

            sage: E = EllipticCurve("17a1")
            sage: L = E.padic_lseries(3)
            sage: L.series(4,prec=1)
            alpha^-2 + alpha^-1 + 2 + 2*alpha + ... + O(alpha^38) + O(T)

        It works also for `p=2`::

            sage: E = EllipticCurve("11a1")
            sage: lp = E.padic_lseries(2)
            sage: lp.series(10)
            O(alpha^-3) + (alpha^-4 + O(alpha^-3))*T + (alpha^-4 + O(alpha^-3))*T^2 + (alpha^-5 + alpha^-4 + O(alpha^-3))*T^3 + (alpha^-4 + O(alpha^-3))*T^4 + O(T^5)
        """
        n = ZZ(n)
        if n < 1:
            raise ValueError("n (=%s) must be a positive integer" % n)
        if self._p == 2 and n == 1:
            raise ValueError("n (=%s) must be at least 2 when p=2" % n)
        if prec < 1:
            raise ValueError("Insufficient precision (%s)" % prec)

        # check if the conditions on quadratic_twist are satisfied
        D = ZZ(quadratic_twist)
        if D != 1:
            if eta != 0:
                raise NotImplementedError("quadratic twists only implemented for the 0th Teichmueller component")
            if D % 4 == 0:
                d = D//4
                if not d.is_squarefree() or d % 4 == 1:
                    raise ValueError("quadratic_twist (=%s) must be a fundamental discriminant of a quadratic field"%D)
            else:
                if not D.is_squarefree() or D % 4 != 1:
                    raise ValueError("quadratic_twist (=%s) must be a fundamental discriminant of a quadratic field" % D)
            if gcd(D, self._E.conductor()) != 1:
                for ell in prime_divisors(D):
                    if valuation(self._E.conductor(), ell) > valuation(D, ell):
                        raise ValueError("cannot twist a curve of conductor (=%s) by the quadratic twist (=%s)." % (self._E.conductor(), D))

        p = self._p
        eta = ZZ(eta) % (p - 1) if p != 2 else ZZ(eta) % 2

        if prec == 1:
            if eta == 0:
                # trac 15737: if we only ask for the leading term we don't
                # need to do any sum as L_p(E,0) = (1-1/alpha)^2 * m(0) (good case)
                # set prec arbitrary to 20.
                alpha = self.alpha(prec=20)
                K = alpha.parent()
                R = PowerSeriesRing(K,'T',1)
                L = self.modular_symbol(0, sign=+1, quadratic_twist=D)
                L *= (1-1/self.alpha())**2
                L /= self._quotient_of_periods_to_twist(D)*self._E.real_components()
                L = R(L, 1)
                return L
            else:
                # here we need some sums anyway
                bounds = self._prec_bounds(n,prec)
                alphaadic_prec = 20
        else:
            prec = min(p**(n-1), prec)
            bounds = self._prec_bounds(n,prec)
            alphaadic_prec = max(bounds[1:]) + 5

        padic_prec = alphaadic_prec//2+1
        verbose("using alpha-adic precision of %s" % padic_prec)
        ans = self._get_series_from_cache(n, prec, quadratic_twist,eta)
        if ans is not None:
            verbose("found series in cache")
            return ans

        alpha = self.alpha(prec=padic_prec)
        K = alpha.parent()
        R = PowerSeriesRing(K,'T',prec)
        T = R(R.gen(), prec)
        L = R(0)
        one_plus_T_factor = R(1)
        gamma_power = 1
        teich = self.teichmuller(padic_prec)
        if p == 2:
            teich = [0, 1,-1]
            gamma = 5
            p_power = 2**(n-2)
            a_range = 3
        else:
            teich = self.teichmuller(padic_prec)
            gamma = 1+ p
            p_power = p**(n-1)
            a_range = p
        si = 1-2*(eta % 2)

        verbose("Now iterating over %s summands"%((p-1)*p_power))
        verbose_level = get_verbose()
        count_verb = 0
        for j in range(p_power):
            s = K(0)
            if verbose_level >= 2 and j/p_power*100 > count_verb + 3:
                verbose("%.2f percent done"%(float(j)/p_power*100))
                count_verb += 3
            for a in range(1,a_range):
                b = teich[a] * gamma_power
                s += teich[a]**eta * self.measure(b, n, padic_prec, quadratic_twist=D, sign=si)
            L += s * one_plus_T_factor
            one_plus_T_factor *= 1+T
            gamma_power *= gamma

        # Now create series but with each coefficient truncated
        # so it is proven correct:
        # the coefficients are now treated as alpha-adic numbers (trac 20254)
        L = R(L,prec)
        aj = L.list()
        if aj:
            bj = [aj[0].add_bigoh(2*(padic_prec-2))]
            j = 1
            while j < len(aj):
                bj.append( aj[j].add_bigoh(bounds[j]) )
                j += 1
            L = R(bj, prec)
        L /= self._quotient_of_periods_to_twist(D)
        if si == +1:
            L /= self._E.real_components()
        self._set_series_in_cache(n, prec, quadratic_twist, eta, L)
        return L

    power_series = series

    def is_ordinary(self):
        r"""
        Return ``True`` if the elliptic curve that this L-function is attached
        to is ordinary.

        EXAMPLES::

            sage: L = EllipticCurve('11a').padic_lseries(19)
            sage: L.is_ordinary()
            False
        """
        return False

    def is_supersingular(self):
        r"""
        Return ``True`` if the elliptic curve that this L function is attached
        to is supersingular.

        EXAMPLES::

            sage: L = EllipticCurve('11a').padic_lseries(19)
            sage: L.is_supersingular()
            True
        """
        return True

    def _prec_bounds(self, n, prec):
        r"""
        A helper function not designed for direct use.

        It returns the `\alpha`-adic precisions of the approximation
        to the `p`-adic L-function.

        EXAMPLES::

            sage: E = EllipticCurve('11a1')
            sage: Lp = E.padic_lseries(19)
            sage: Lp._prec_bounds(3,5)
            [+Infinity, -1, -1, -1, -1]
            sage: Lp._prec_bounds(2,5)
            [+Infinity, -2, -2, -2, -2]
            sage: Lp._prec_bounds(10,5)
            [+Infinity, 6, 6, 6, 6]
        """
        if self._p == 2:
            e = self._e_bounds(n - 2, prec)
        else:
            e = self._e_bounds(n - 1, prec)
        c0 = ZZ(n + 2)
        return [infinity] + [2 * e[j] - c0 for j in range(1, len(e))]

    def _poly(self, a):
        """
        Given an element a in Qp[alpha] this returns the list
        containing the two coordinates in Qp.

        EXAMPLES::

            sage: E = EllipticCurve("14a1")
            sage: lp = E.padic_lseries(5)
            sage: K = lp.alpha().parent()
            sage: a = K(5)
            sage: a
            4*alpha^2 + alpha^4 + O(alpha^42)
            sage: lp._poly(a)
            [5 + O(5^21), O(5^21)]
        """
        # this should be implemented in elements of Eisenstein rings at some point trac 20248

        if a.is_zero():
            return [0,0]
        v, k = a._ntl_rep_abs()
        K = a.base_ring()
        pi = K.uniformiser()
        v0 =  K(v[0]._sage_()) * pi**k
        v1 =  K(v[1]._sage_()) * pi**k
        alpha = a.parent().gen()
        assert v0 + v1*alpha == a
        return [ v0, v1 ]

    def Dp_valued_series(self, n=3, quadratic_twist=+1, prec=5):
        r"""
        Return a vector of two components which are p-adic power series.

        The answer v is such that

            `(1-\varphi)^{-2}\cdot L_p(E,T) =` ``v[1]`` `\cdot \omega +` ``v[2]`` `\cdot \varphi(\omega)`

        as an element of the Dieudonné module `D_p(E) = H^1_{dR}(E/\QQ_p)` where
        `\omega` is the invariant differential and `\varphi` is the Frobenius on `D_p(E)`.

        According to the `p`-adic Birch and Swinnerton-Dyer
        conjecture [BP1993]_ this function has a zero of order
        rank of `E(\QQ)` and it's leading term is contains the order of
        the Tate-Shafarevich group, the Tamagawa numbers, the order of the
        torsion subgroup and the `D_p`-valued `p`-adic regulator.

        INPUT:

        -  ``n`` -- (default: 3) a positive integer
        -  ``prec`` -- (default: 5) a positive integer

        EXAMPLES::

            sage: E = EllipticCurve('14a')
            sage: L = E.padic_lseries(5)
            sage: L.Dp_valued_series(4)  # long time (9s on sage.math, 2011)
            (1 + 4*5 + O(5^2) + (4 + O(5))*T + (1 + O(5))*T^2 + (4 + O(5))*T^3 + (2 + O(5))*T^4 + O(T^5), 5^2 + O(5^3) + O(5^2)*T + (4*5 + O(5^2))*T^2 + (2*5 + O(5^2))*T^3 + (2 + 2*5 + O(5^2))*T^4 + O(T^5))
        """
        E = self._E
        p = self._p
        lps = self.series(n, quadratic_twist=quadratic_twist, prec=prec)

        # now split up the series in two lps = G + H * alpha
        R = lps.base_ring().base_ring() # Qp
        QpT , T = PowerSeriesRing(R, 'T', prec).objgen()
        Gli = []
        Hli = []
        for n in range(lps.prec()):
            v = self._poly(lps[n])
            Gli.append(v[0])
            Hli.append(v[1])
        G = QpT(Gli, prec)
        H = QpT(Hli, prec)

        # now compute phi
        phi = matrix.matrix([[0,-1/p],[1,E.ap(p)/p]])
        lpv = vector([G  + (E.ap(p))*H  , - R(p) * H ])  # this is L_p
        eps = (1-phi)**(-2)
        resu = lpv*eps.transpose()
        return resu

    def frobenius(self, prec=20, algorithm="mw"):
        r"""
        Return a geometric Frobenius `\varphi` on the Dieudonné module `D_p(E)`
        with respect to the basis `\omega`, the invariant differential, and `\eta=x\omega`.

        It satisfies  `\varphi^2 - a_p/p\, \varphi + 1/p = 0`.

        INPUT:

        - ``prec`` - (default: 20) a positive integer

        - ``algorithm`` - either 'mw' (default) for Monsky-Washnitzer
          or 'approx' for the algorithm described by Bernardi and Perrin-Riou
          (much slower and not fully tested)

        EXAMPLES::

            sage: E = EllipticCurve('14a')
            sage: L = E.padic_lseries(5)
            sage: phi = L.frobenius(5)
            sage: phi
            [                  2 + 5^2 + 5^4 + O(5^5)    3*5^-1 + 3 + 5 + 4*5^2 + 5^3 + O(5^4)]
            [      3 + 3*5^2 + 4*5^3 + 3*5^4 + O(5^5) 3 + 4*5 + 3*5^2 + 4*5^3 + 3*5^4 + O(5^5)]
            sage: -phi^2
            [5^-1 + O(5^4)        O(5^4)]
            [       O(5^5) 5^-1 + O(5^4)]
        """
        E = self._E
        p = self._p
        if algorithm != "mw" and algorithm !="approx":
            raise ValueError("Unknown algorithm %s."%algorithm)
        if algorithm == "approx":
            return self.__phi_bpr(prec=prec)
        if p < 4 and algorithm == "mw":
            print("Warning: If this fails try again using algorithm=\"approx\"")
        Ew = E.integral_short_weierstrass_model()
        adjusted_prec = sage.schemes.hyperelliptic_curves.monsky_washnitzer.adjusted_prec(p, prec)
        modprecring = Integers(p**adjusted_prec)
        output_ring = Qp(p, prec)
        R, x = PolynomialRing(modprecring, 'x').objgen()
        Q = x**3 + modprecring(Ew.a4()) * x + modprecring(Ew.a6())
        trace = Ew.ap(p)
        fr = sage.schemes.hyperelliptic_curves.monsky_washnitzer.matrix_of_frobenius(Q, p, adjusted_prec, trace)
        fr = matrix.matrix(output_ring,2,2,fr)

        # return a vector for PARI's ellchangecurve to pass from e1 to e2
        def isom(e1,e2):
            if not e1.is_isomorphic(e2):
                raise ValueError("Curves must be isomorphic.")
            usq = (e1.discriminant()/e2.discriminant()).nth_root(6)
            u = usq.sqrt()
            s = (u   *  e2.a1() - e1.a1() )/ZZ(2)
            r = (usq *  e2.a2() - e1.a2() + s**2 + e1.a1()*s)/ZZ(3)
            t = (u**3 * e2.a3() - e1.a3() - e1.a1()*r)/ZZ(2)
            return [u,r,s,t]

        v = isom(E,Ew)
        u = v[0]
        r = v[1]

        # change basis
        A = matrix.matrix([[u,-r/u],[0,1/u]])
        frn = A * fr * A**(-1)
        return 1/p*frn



    def __phi_bpr(self, prec=0):
        r"""
        This returns a geometric Frobenius `\varphi` on the Dieudonné module `D_p(E)`
        with respect to the basis `\omega`, the invariant differential, and `\eta=x\omega`.
        It satisfies  `\varphi^2 - a_p/p\, \varphi + 1/p = 0`.

        The algorithm used here is described in bernardi-perrin-riou on page 232.

        .. WARNING::

            This function has not been sufficiently tested. It is very slow.

        EXAMPLES::

            sage: E = EllipticCurve('11a1')
            sage: lp = E.padic_lseries(19)
            sage: lp.frobenius(prec=1,algorithm="approx")   #indirect doctest
            [          O(19^0) 4*19^-1 + O(19^0)]
            [       14 + O(19)           O(19^0)]

            sage: E = EllipticCurve('17a1')
            sage: lp = E.padic_lseries(3)
            sage: lp.frobenius(prec=3,algorithm="approx")
            [             O(3) 2*3^-1 + 2 + O(3)]
            [       1 + O(3^2)              O(3)]
            sage: lp.frobenius(prec=5,algorithm="approx")
            [             3 + O(3^2) 2*3^-1 + 2 + 3 + O(3^2)]
            [     1 + 2*3^2 + O(3^3)            2*3 + O(3^2)]


        """
        E = self._E
        p = self._p
        if prec > 10:
            print("Warning: Very large value for the precision.")
        if prec == 0:
            prec = floor((log(10000)/log(p)))
            verbose("prec set to %s"%prec)
        eh = E.formal()
        om = eh.differential(prec = p**prec+3)
        verbose("differential computed")
        xt = eh.x(prec=p**prec + 3)
        et = xt*om
        # c_(p^k) = cs[k] d...
        cs = [om[p**k-1] for k in range(prec + 1)]
        ds = [et[p**k-1] for k in range(prec + 1)]
        delta = 0
        dpr = 0
        gamma = 0
        dga = 0
        for k in range(1,prec+1):
            # this is the equation eq[0]*x+eq[1]*y+eq[2] == 0
            # such that delta_ = delta + d^dpr*x ...
            eq = [(p**dpr*cs[k]) % p**k,(-p**dga*ds[k]) % p**k , (delta*cs[k]-gamma*ds[k]-cs[k-1]) % p**k ]
            verbose("valuations : %s" % ([x.valuation(p) for x in eq]))
            v = min([x.valuation(p) for x in eq])
            if v == infinity:
                verbose("no new information at step k=%s"%k)
            else:
                eq = [ZZ(x/p**v) for x in eq]
                verbose("renormalised eq mod p^%s is now %s"%(k-v,eq))
                if eq[0].valuation(p) == 0:
                    l = min(eq[1].valuation(p),k-v)
                    if l == 0:
                        verbose("not uniquely determined at step k=%s"%k)
                    else:
                        ainv = eq[0].inverse_mod(p**l)
                        delta = delta - eq[2]*ainv*p**dpr
                        dpr = dpr + l
                        delta = delta % p**dpr
                        verbose("delta_prec increased to %s\n delta is now %s"%(dpr,delta))
                elif eq[1].valuation(p) == 0:
                    l = min(eq[0].valuation(p),k-v)
                    ainv = eq[1].inverse_mod(p**l)
                    gamma = gamma - eq[2]*ainv*p**dga
                    dga = dga + l
                    gamma = gamma % p**dga
                    verbose("gamma_prec increased to %s\n gamma is now %s"%(dga,gamma))
                else:
                    raise RuntimeError("Bug: no delta or gamma can exist")

        # end of approximation of delta and gamma
        R = Qp(p,max(dpr,dga)+1)
        delta = R(delta,absprec=dpr)
        gamma = R(gamma,absprec=dga)
        verbose("result delta = %s\n      gamma = %s\n check : %s"%(delta,gamma, [Qp(p,k)(delta * cs[k] - gamma * ds[k] - cs[k-1]) for k in range(1,prec+1)] ))
        a = delta
        c = -gamma
        d = E.ap(p) - a
        b = (-1/p+a*d)/c
        phi = matrix.matrix([[a,b],[c,d]])
        return phi

    def bernardi_sigma_function(self, prec=20):
        r"""
        Return the  `p`-adic sigma function of Bernardi in terms of `z = log(t)`.

        This is the same as ``padic_sigma`` with ``E2 = 0``.

        EXAMPLES::

            sage: E = EllipticCurve('14a')
            sage: L = E.padic_lseries(5)
            sage: L.bernardi_sigma_function(prec=5) # Todo: some sort of consistency check!?
            z + 1/24*z^3 + 29/384*z^5 - 8399/322560*z^7 - 291743/92897280*z^9 + O(z^10)
        """
        E = self._E

        Eh = E.formal()
        lo = Eh.log(prec + 5)
        F = lo.reverse()

        S = LaurentSeriesRing(QQ,'z')
        z = S.gen()
        F = F(z)
        xofF = Eh.x(prec + 2)(F)
        #r =  ( E.a1()**2 + 4*E.a2() ) / ZZ(12)
        g = (1/z**2 - xofF ).power_series()
        h = g.integral().integral()
        sigma_of_z = z.power_series() * h.exp()

        return sigma_of_z

    def Dp_valued_height(self,prec=20):
        r"""
        Return the canonical `p`-adic height with values in the Dieudonné module `D_p(E)`.

        It is defined to be

            `h_{\eta} \cdot \omega - h_{\omega} \cdot \eta`

        where `h_{\eta}` is made out of the sigma function of Bernardi and
        `h_{\omega}` is `log_E^2`.

        The answer ``v`` is given as ``v[1]*omega + v[2]*eta``.
        The coordinates of ``v`` are dependent of the
        Weierstrass equation.

        EXAMPLES::

            sage: E = EllipticCurve('53a')
            sage: L = E.padic_lseries(5)
            sage: h = L.Dp_valued_height(7)
            sage: h(E.gens()[0])
            (3*5 + 5^2 + 2*5^3 + 3*5^4 + 4*5^5 + 5^6 + 5^7 + O(5^8), 5^2 + 4*5^4 + 2*5^7 + 3*5^8 + O(5^9))
        """
        E = self._E
        p = self._p
        Ehat = E.formal()
        elog = Ehat.log(prec + Integer(3))

        # we will have to do it properly with David Harvey's _multiply_point()
        n = arith.LCM(E.tamagawa_numbers())
        n = arith.LCM(n, E.Np(p)) # allowed here because E has good reduction at p

        def height(P,check=True):
            if P.is_finite_order():
                return Qp(p,prec)(0)
            if check:
                assert P.curve() == E, 'the point P must lie on the curve from which the height function was created'

            Q = n * P
            tt = - Q[0]/Q[1]
            R = Qp(p,prec+5)
            tt = R(tt)
            zz = elog(tt)

            homega = -zz**2/n**2

            eQ = denominator(Q[1])/denominator(Q[0])
            si = self.bernardi_sigma_function(prec=prec+4)
            heta =  2 * log(si(zz)/eQ) / n**2

            R = Qp(p,prec)

            return vector([-R(heta),R(homega)])

        return height

    def Dp_valued_regulator(self, prec=20, v1=0, v2=0):
        r"""
        Return the canonical `p`-adic regulator with values in the Dieudonné module `D_p(E)`
        as defined by Perrin-Riou using the `p`-adic height with values in `D_p(E)`.

        The result is written in the basis `\omega`, `\varphi(\omega)`, and hence the
        coordinates of the result are independent of the chosen Weierstrass equation.

        .. NOTE::

            The definition here is corrected with respect to
            Perrin-Riou's article [PR2003]_. See [SW2013]_.

        EXAMPLES::

            sage: E = EllipticCurve('43a')
            sage: L = E.padic_lseries(7)
            sage: L.Dp_valued_regulator(7)
            (5*7 + 6*7^2 + 4*7^3 + 4*7^4 + 7^5 + 4*7^7 + O(7^8), 4*7^2 + 2*7^3 + 3*7^4 + 7^5 + 6*7^6 + 4*7^7 + O(7^8))
        """
        p = self._p
        E = self._E

        h =  self.Dp_valued_height(prec=prec)

        # this is the height_{v} (P) for a v in D_p
        def hv(vec,P):
            hP = h(P)
            return - vec[0]*hP[1] +vec[1]*hP[0]

        #    def hvpairing(vec,P,Q):
        #        return (hv(vec,    P+Q) - hv(vec,P)-hv(vec,Q))/2
        K = Qp(p, prec)

        if v1 == 0 and v2 == 0:
            v1 = vector([K(0), K(1)])  # that is eta
            v2 = vector([K(-1), K(1)])  # and this is eta-omega.
        #                      the rest should not depend on this choice
        #                      as long as it is outside Q_p * omega

        rk = E.rank()
        if rk == 0:
            return vector([K(1), K(0)])

        basis = E.gens()

        def regv(vec):
            M = matrix.matrix(K, rk, rk, 0)
            point_height = [hv(vec, P) for P in basis]
            for i in range(rk):
                for j in range(i+1, rk):
                    M[i, j] = M[j, i] = (hv(vec,basis[i] + basis[j])- point_height[i] - point_height[j] )/2
            for i in range(rk):
                M[i, i] = point_height[i]

            return M.determinant()

        def Dp_pairing(vec1,vec2):
            return (vec1[0]*vec2[1]-vec1[1]*vec2[0])

        omega_vec = vector([K(1),K(0)])

        # note the correction here with respect to Perrin-Riou's definition.
        # only this way the result will be independent of the choice of v1 and v2.
        reg1 = regv(v1) / Dp_pairing(omega_vec, v1)**(rk - 1)

        reg2 = regv(v2) / Dp_pairing(omega_vec, v2)**(rk - 1)

        # the regulator in the basis omega,eta
        reg_oe = (reg1 * v2 - reg2 * v1 ) / Dp_pairing(v2, v1)

        if p < 5:
            phi = self.frobenius(min(6, prec), algorithm="approx")
        else:
            phi = self.frobenius(prec + 2, algorithm="mw")

        c = phi[1, 0]  # this is the 'period' [omega,phi(omega)]
        a = phi[0, 0]

        return vector([reg_oe[0] - a/c*reg_oe[1],reg_oe[1]/c])
