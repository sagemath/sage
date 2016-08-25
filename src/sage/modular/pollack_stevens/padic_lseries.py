# -*- coding: utf-8 -*-
r"""
`p`-adic `L`-series attached to overconvergent eigensymbols
""" ## mm TODO
#*****************************************************************************
#       Copyright (C) 2012 Robert Pollack <rpollack@math.bu.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import print_function
from __future__ import absolute_import
from sage.rings.padics.all import pAdicField
from sage.rings.all import ZZ, QQ
from sage.rings.power_series_ring import PowerSeriesRing
from sage.rings.big_oh import O
from sage.arith.all import binomial, gcd, kronecker
from sage.rings.padics.precision_error import PrecisionError

from sage.structure.sage_object import SageObject
from .sigma0 import Sigma0
from .fund_domain import M2Z


class pAdicLseries(SageObject):
    r"""
    The `p`-adic `L`-series associated to an overconvergent eigensymbol.

    INPUT:

    - ``symb`` -- an overconvergent eigensymbol
    - ``gamma`` -- topological generator of `1 + pZ_p` (default: `1+p` or 5 if `p=2`)
    - ``quadratic_twist`` -- conductor of quadratic twist `\chi` (default: 1)
    - ``precision`` -- if None (default) is specified, the correct precision bound is
      computed and the answer is returned modulo that accuracy

    EXAMPLES::

        sage: E = EllipticCurve('37a')
        sage: p = 5
        sage: prec = 4
        sage: L = E.padic_lseries(p, implementation="pollackstevens", precision=prec) # long time
        sage: L[1]                # long time
        1 + 4*5 + 2*5^2 + O(5^3)
        sage: L.series(prec,3)    # long time
        O(5^4) + (1 + 4*5 + 2*5^2 + O(5^3))*T + (3 + O(5^2))*T^2 + O(T^3)

    ::

        sage: from sage.modular.pollack_stevens.padic_lseries import pAdicLseries
        sage: E = EllipticCurve('20a')
        sage: phi = E.pollack_stevens_modular_symbol()
        sage: Phi = phi.p_stabilize_and_lift(3, 4) # long time
        sage: L = pAdicLseries(Phi)                # long time
        sage: L.series(prec, 4)                    # long time
        2*3 + O(3^4) + (3 + O(3^2))*T + (2 + O(3))*T^2 + O(3^0)*T^3 + O(T^4)

    An example of a `p`-adic `L`-series associated to a modular
    abelian surface. This is not tested as it takes too long.::

        sage: from sage.modular.pollack_stevens.space import ps_modsym_from_simple_modsym_space
        sage: from sage.modular.pollack_stevens.padic_lseries import pAdicLseries
        sage: A = ModularSymbols(103,2,1).cuspidal_submodule().new_subspace().decomposition()[0]
        sage: p = 19
        sage: prec = 4
        sage: phi = ps_modsym_from_simple_modsym_space(A)
        sage: ap = phi.Tq_eigenvalue(p,prec)
        sage: c1,c2 = phi.completions(p,prec)
        sage: phi1,psi1 = c1
        sage: phi2,psi2 = c2
        sage: phi1p = phi1.p_stabilize_and_lift(p,ap = psi1(ap), M = prec) # not tested - too long
        sage: L1 = pAdicLseries(phi1p)                                     # not tested - too long
        sage: phi2p = phi2.p_stabilize_and_lift(p,ap = psi2(ap), M = prec) # not tested - too long
        sage: L2  = pAdicLseries(phi2p)                                    # not tested - too long
        sage: L1[1]*L2[1]                                                  # not tested - too long
        13 + 9*19 + 18*19^2 + O(19^3)
    """

    def __init__(self, symb, gamma=None, quadratic_twist=1, precision=None):
        r"""
        Initialize the class

        EXAMPLE::

            sage: from sage.modular.pollack_stevens.padic_lseries import pAdicLseries
            sage: E = EllipticCurve('11a3')
            sage: phi = E.pollack_stevens_modular_symbol()
            sage: p = 11
            sage: prec = 3
            sage: Phi = phi.lift(p, prec,eigensymbol=True) # long time
            sage: L = pAdicLseries(Phi)                    # long time
            sage: L.series(3, prec=3)                      # long time
            O(11^3) + (2 + 5*11 + O(11^2))*T + (10 + O(11))*T^2 + O(T^3)

            sage: TestSuite(L).run()                       # long time
        """
        self._coefficients = {}

        if symb.parent().prime() is None:
            raise ValueError("Not a p-adic overconvergent modular symbol.")

        self._symb = symb

        if gamma is None:
            p = self._symb.parent().prime()
            if p == 2:
                gamma = 1 + 4
            else:
                gamma = 1 + self._symb.parent().prime()

        self._gamma = gamma
        self._quadratic_twist = quadratic_twist
        self._precision = precision
        self._cinf = ZZ(1) # is set when called for an elliptic curve

    def __getitem__(self, n):
        r"""
        Return the `n`-th coefficient of the `p`-adic `L`-series

        EXAMPLES::

            sage: E = EllipticCurve('14a5')
            sage: L = E.padic_lseries(7,implementation="pollackstevens",precision=5) # long time
            sage: L[0]                                   # long time
            O(7^5)
            sage: L[1]                                   # long time
            5 + 5*7 + 2*7^2 + 2*7^3 + O(7^4)
         """

        if n in self._coefficients:
            return self._coefficients[n]
        else:
            p = self.prime()
            symb = self.symbol()
            # ap = symb.Tq_eigenvalue(p)
            gamma = self._gamma
            precision = self._precision

            S = QQ[['z']]
            z = S.gen()
            M = symb.precision_relative()
            K = pAdicField(p, M)
            dn = 0
            if n == 0:
                precision = M
                lb = [1] + [0 for a in range(M - 1)]
            else:
                lb = log_gamma_binomial(p, gamma, z, n, 2 * M)
                if precision is None:
                    precision = min([j + lb[j].valuation(p)
                                     for j in range(M, len(lb))])
                lb = [lb[a] for a in range(M)]

            for j in range(len(lb)):
                cjn = lb[j]
                temp = sum((ZZ(K.teichmuller(a)) ** (-j))
                           * self._basic_integral(a, j) for a in range(1, p))
                dn = dn + cjn * temp
            self._coefficients[n] = dn.add_bigoh(precision)
            self._coefficients[n] /= self._cinf
            return self._coefficients[n]

    def __cmp__(self, other):
        r"""
        Compare ``self`` and ``other``.

        EXAMPLE::

            sage: E = EllipticCurve('11a')
            sage: L = E.padic_lseries(11,implementation="pollackstevens",precision=6) # long time
            sage: L == loads(dumps(L)) # indirect doctest long time
            True
        """
        return (cmp(type(self), type(other))
                or cmp(self._symb, other._symb)
                or cmp(self._quadratic_twist, other._quadratic_twist)
                or cmp(self._gamma, other._gamma)
                or cmp(self._precision, other._precision))

    def symbol(self):
        r"""
        Return the overconvergent modular symbol

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.padic_lseries import pAdicLseries
            sage: E = EllipticCurve('21a4')
            sage: phi = E.pollack_stevens_modular_symbol()
            sage: Phi = phi.p_stabilize_and_lift(2,5)   # long time
            sage: L = pAdicLseries(Phi)                 # long time
            sage: L.symbol()                              # long time
            Modular symbol of level 42 with values in Space of 2-adic distributions with k=0 action and precision cap 15
            sage: L.symbol() is Phi                       # long time
            True
        """
        return self._symb

    def prime(self):
        r"""
        Return the prime `p` as in `p`-adic `L`-series.

        EXAMPLES::

            sage: E = EllipticCurve('19a')
            sage: L = E.padic_lseries(19, implementation="pollackstevens",precision=6) # long time
            sage: L.prime()                   # long time
            19
        """
        return self._symb.parent().prime()

    def quadratic_twist(self):
        r"""
        Return the discriminant of the quadratic twist

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.padic_lseries import pAdicLseries
            sage: E = EllipticCurve('37a')
            sage: phi = E.pollack_stevens_modular_symbol()
            sage: Phi = phi.lift(37,4)
            sage: L = pAdicLseries(Phi, quadratic_twist=-3)
            sage: L.quadratic_twist()
            -3
        """
        return self._quadratic_twist

    def _repr_(self):
        r"""
        Return the string representation

        EXAMPLES::

            sage: E = EllipticCurve('14a2')
            sage: L = E.padic_lseries(3, implementation="pollackstevens", precision=4)  # long time
            sage: L._repr_()                           # long time
            '3-adic L-series of Modular symbol of level 42 with values in Space of 3-adic distributions with k=0 action and precision cap 8'
        """
        return "%s-adic L-series of %s" % (self.prime(), self.symbol())

    def series(self, n, prec=5):
        r"""
        Return the `n`-th approximation to the `p`-adic `L`-series
        associated to self, as a power series in `T` (corresponding to
        `\gamma-1` with `\gamma` the chosen generator of `1+p\ZZ_p`).

        INPUT:

        - ``n`` -- ## mm TODO

        - ``prec`` -- (default 5) is the precision of the power series

        EXAMPLES::

            sage: E = EllipticCurve('14a2')
            sage: p = 3
            sage: prec = 6
            sage: L = E.padic_lseries(p,implementation="pollackstevens",precision=prec) # long time
            sage: L.series(prec, 4)          # long time
            2*3 + 3^4 + 3^5 + O(3^6) + (2*3 + 3^2 + O(3^4))*T + (2*3 + O(3^2))*T^2 + (3 + O(3^2))*T^3 + O(T^4)

            sage: E = EllipticCurve("15a3")
            sage: L = E.padic_lseries(5,implementation="pollackstevens",precision=15)  # long time
            sage: L.series(10, 3)            # long time
            O(5^15) + (2 + 4*5^2 + 3*5^3 + 5^5 + 2*5^6 + 3*5^7 + 3*5^8 + 2*5^9 + 2*5^10 + 3*5^11 + 5^12 + O(5^13))*T + (4*5 + 4*5^3 + 3*5^4 + 4*5^5 + 3*5^6 + 2*5^7 + 5^8 + 4*5^9 + 3*5^10 + O(5^11))*T^2 + O(T^3)

            sage: E = EllipticCurve("79a1")
            sage: L = E.padic_lseries(2,implementation="pollackstevens",precision=10) # not tested
            sage: L.series(10, 4)  # not tested
            O(2^9) + (2^3 + O(2^4))*T + O(2^0)*T^2 + (O(2^-3))*T^3 + O(T^4)
        """
        p = self.prime()
        M = self.symbol().precision_relative()
        K = pAdicField(p, M)
        R = PowerSeriesRing(K, names='T')
        T = R.gens()[0]
        R.set_default_prec(n)
        return (sum(self[i] * T ** i for i in range(prec))).add_bigoh(prec)

    def interpolation_factor(self, ap, chip=1, psi=None):
        r"""
        Return the interpolation factor associated to self.
        This is the `p`-adic multiplier that which appears in
        the interpolation formula of the `p`-adic `L`-function.

        INPUT:

        - ``ap`` -- ## mm TODO

        - ``chip`` --

        - ``psi`` --

        OUTPUT: a `p`-adic number

        EXAMPLES::

            sage: E = EllipticCurve('19a2')
            sage: L = E.padic_lseries(3,implementation="pollackstevens",precision=6)  # long time
            sage: ap = E.ap(3)               # long time
            sage: L.interpolation_factor(ap) # long time
            3^2 + 3^3 + 2*3^5 + 2*3^6 + O(3^7)

        Comparing against a different implementation::

            sage: L = E.padic_lseries(3)
            sage: (1-1/L.alpha(prec=4))^2
            3^2 + 3^3 + O(3^5)
        """
        M = self.symbol().precision_relative()
        p = self.prime()
        if p == 2:
            R = pAdicField(2, M + 1)
        else:
            R = pAdicField(p, M)
        if psi is not None:
            ap = psi(ap)
        ap = ap * chip
        sdisc = R(ap ** 2 - 4 * p).sqrt()
        v0 = (R(ap) + sdisc) / 2
        v1 = (R(ap) - sdisc) / 2
        if v0.valuation() > 0:
            v0, v1 = v1, v0
        alpha = v0
        return (1 - 1 / alpha) ** 2

    def _basic_integral(self, a, j):
        r"""
        Return `\int_{a+pZ_p} (z-{a})^j d\Phi(0-infty)`
        -- see formula [Pollack-Stevens, sec 9.2]

        INPUT:

        - ``a`` -- integer in range(p)
        - ``j`` -- integer in range(self.symbol().precision_relative())

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.padic_lseries import pAdicLseries
            sage: E = EllipticCurve('11a3')
            sage: L = E.padic_lseries(5, implementation="pollackstevens", precision=4) #long time
            sage: L._basic_integral(1,2) # long time
            2*5^2 + 5^3 + O(5^4)
        """
        symb = self.symbol()
        M = symb.precision_relative()
        if j > M:
            raise PrecisionError("Too many moments requested")
        p = self.prime()
        ap = symb.Tq_eigenvalue(p)
        D = self._quadratic_twist
        ap = ap * kronecker(D, p)
        K = pAdicField(p, M)
        symb_twisted = symb.evaluate_twisted(a, D)
        return ( sum(binomial(j, r)
                 * ((a - ZZ(K.teichmuller(a))) ** (j - r))
                 * (p ** r)
                 * symb_twisted.moment(r) for r in range(j + 1))
                 / ap  )


def log_gamma_binomial(p, gamma, z, n, M):
    r"""
    Return the list of coefficients in the power series
    expansion (up to precision `M`) of `{\log_p(z)/\log_p(\gamma) \choose n}`

    INPUT:

    - ``p`` --  prime
    - ``gamma`` -- topological generator, e.g. `1+p`
    - ``z`` -- variable
    - ``n`` -- nonnegative integer
    - ``M`` -- precision

    OUTPUT:

    The list of coefficients in the power series expansion of
    `{\log_p(z)/\log_p(\gamma) \choose n}`

    EXAMPLES::

        sage: R.<z> = QQ['z']
        sage: from sage.modular.pollack_stevens.padic_lseries import log_gamma_binomial
        sage: log_gamma_binomial(5,1+5,z,2,4)
        [0, -3/205, 651/84050, -223/42025]
        sage: log_gamma_binomial(5,1+5,z,3,4)
        [0, 2/205, -223/42025, 95228/25845375]
    """
    L = sum([ZZ(-1) ** j / j * z ** j for j in range(1, M)])  # log_p(1+z)
    loggam = L / (L(gamma - 1))   # log_{gamma}(1+z)= log_p(1+z)/log_p(gamma)
    return z.parent()(binomial(loggam, n)).truncate(M).list()
