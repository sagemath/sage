# sage_setup: distribution = sagemath-schemes
# sage.doctest: needs sage.rings.padics
#
# All these methods are imported in EllipticCurve_rational_field,
# so there is no reason to add this module to the documentation.
r"""
Miscellaneous `p`-adic methods
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

import math

from sage.arith.functions import lcm as LCM
from sage.arith.misc import valuation
from sage.matrix.constructor import Matrix as matrix
from sage.misc.cachefunc import cached_method
from sage.misc.lazy_import import lazy_import
from sage.misc.misc import newton_method_sizes
from sage.rings.big_oh import O
from sage.rings.finite_rings.finite_field_constructor import FiniteField as GF
from sage.rings.finite_rings.integer_mod_ring import IntegerModRing as Integers
from sage.rings.integer import Integer
from sage.rings.laurent_series_ring import LaurentSeriesRing
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.power_series_ring import PowerSeriesRing
from sage.rings.rational_field import RationalField
import sage.schemes.hyperelliptic_curves.hypellfrob
import sage.schemes.hyperelliptic_curves.monsky_washnitzer

lazy_import('sage.rings.padics.factory', ['Qp', 'Zp'])


from . import padic_lseries as plseries

sqrt = math.sqrt


def __check_padic_hypotheses(self, p):
    r"""
    Helper function that determines if `p`
    is an odd prime of good ordinary reduction.

    EXAMPLES::

        sage: E = EllipticCurve('11a1')
        sage: from sage.schemes.elliptic_curves.padics import __check_padic_hypotheses
        sage: __check_padic_hypotheses(E,5)
        5
        sage: __check_padic_hypotheses(E,29)
        Traceback (most recent call last):
        ...
        ArithmeticError: p must be a good ordinary prime
    """
    p = Integer(p)
    if not p.is_prime():
        raise ValueError("p = (%s) must be prime" % p)
    if p == 2:
        raise ValueError("p must be odd")
    if self.conductor() % p == 0 or self.ap(p) % p == 0:
        raise ArithmeticError("p must be a good ordinary prime")
    return p


def _normalize_padic_lseries(self, p, normalize, implementation, precision):
    r"""
    Normalize parameters for :meth:`padic_lseries`.

    TESTS::

        sage: from sage.schemes.elliptic_curves.padics import _normalize_padic_lseries
        sage: u = _normalize_padic_lseries(None, 5, None, 'sage', 10)
        sage: v = _normalize_padic_lseries(None, 5, "L_ratio", 'sage', 10)
        sage: u == v
        True
    """
    if implementation == 'eclib':
        if normalize is None:
            normalize = "L_ratio"
    elif implementation == 'sage':
        if normalize is None:
            normalize = "L_ratio"
    elif implementation == 'pollackstevens':
        if precision is None:
            raise ValueError("Must specify precision when using 'pollackstevens'")
        if normalize is not None:
            raise ValueError("The 'normalize' parameter is not used for Pollack-Stevens' overconvergent modular symbols")
    elif implementation == "num":
        if normalize is not None:
            raise ValueError("The 'normalize' parameter is not used for numerical modular symbols")
    else:
        raise ValueError("Implementation should be one of  'sage', 'eclib', 'num' or 'pollackstevens'")
    return (p, normalize, implementation, precision)


@cached_method(key=_normalize_padic_lseries)
def padic_lseries(self, p, normalize=None, implementation='eclib',
                  precision=None):
    r"""
    Return the `p`-adic `L`-series of ``self`` at
    `p`, which is an object whose approx method computes
    approximation to the true `p`-adic `L`-series to
    any desired precision.

    INPUT:

    - ``p`` -- prime

    - ``normalize`` -- 'L_ratio' (default), 'period' or 'none';
      this is describes the way the modular symbols
      are normalized. See modular_symbol for
      more details.

    - ``implementation`` -- 'eclib' (default), 'sage', 'num' or 'pollackstevens';
      Whether to use John Cremona's eclib, the Sage implementation,
      numerical modular symbols
      or Pollack-Stevens' implementation of overconvergent
      modular symbols.

    EXAMPLES::

        sage: E = EllipticCurve('37a')
        sage: L = E.padic_lseries(5); L
        5-adic L-series of Elliptic Curve defined by y^2 + y = x^3 - x over Rational Field
        sage: type(L)
        <class 'sage.schemes.elliptic_curves.padic_lseries.pAdicLseriesOrdinary'>

    We compute the `3`-adic `L`-series of two curves of
    rank `0` and in each case verify the interpolation property
    for their leading coefficient (i.e., value at 0)::

        sage: e = EllipticCurve('11a')
        sage: ms = e.modular_symbol()
        sage: [ms(1/11), ms(1/3), ms(0), ms(oo)]
        [0, -3/10, 1/5, 0]
        sage: ms(0)
        1/5
        sage: L = e.padic_lseries(3)
        sage: P = L.series(5)
        sage: P(0)
        2 + 3 + 3^2 + 2*3^3 + 2*3^5 + 3^6 + O(3^7)
        sage: alpha = L.alpha(9); alpha
        2 + 3^2 + 2*3^3 + 2*3^4 + 2*3^6 + 3^8 + O(3^9)
        sage: R.<x> = QQ[]
        sage: f = x^2 - e.ap(3)*x + 3
        sage: f(alpha)
        O(3^9)
        sage: r = e.lseries().L_ratio(); r
        1/5
        sage: (1 - alpha^(-1))^2 * r
        2 + 3 + 3^2 + 2*3^3 + 2*3^5 + 3^6 + 3^7 + O(3^9)
        sage: P(0)
        2 + 3 + 3^2 + 2*3^3 + 2*3^5 + 3^6 + O(3^7)

    Next consider the curve 37b::

        sage: e = EllipticCurve('37b')
        sage: L = e.padic_lseries(3)
        sage: P = L.series(5)
        sage: alpha = L.alpha(9); alpha
        1 + 2*3 + 3^2 + 2*3^5 + 2*3^7 + 3^8 + O(3^9)
        sage: r = e.lseries().L_ratio(); r
        1/3
        sage: (1 - alpha^(-1))^2 * r
        3 + 3^2 + 2*3^4 + 2*3^5 + 2*3^6 + 3^7 + O(3^9)
        sage: P(0)
        3 + 3^2 + 2*3^4 + 2*3^5 + O(3^6)

    We can use Sage modular symbols instead to compute the `L`-series::

        sage: e = EllipticCurve('11a')
        sage: L = e.padic_lseries(3, implementation = 'sage')
        sage: L.series(5,prec=10)
        2 + 3 + 3^2 + 2*3^3 + 2*3^5 + 3^6 + O(3^7) + (1 + 3 + 2*3^2 + 3^3 + O(3^4))*T + (1 + 2*3 + O(3^4))*T^2 + (3 + 2*3^2 + O(3^3))*T^3 + (2*3 + 3^2 + O(3^3))*T^4 + (2 + 2*3 + 2*3^2 + O(3^3))*T^5 + (1 + 3^2 + O(3^3))*T^6 + (2 + 3^2 + O(3^3))*T^7 + (2 + 2*3 + 2*3^2 + O(3^3))*T^8 + (2 + O(3^2))*T^9 + O(T^10)

    Also the numerical modular symbols can be used.
    This may allow for much larger conductor in some instances::

        sage: E = EllipticCurve([101,103])
        sage: L = E.padic_lseries(5, implementation='num')
        sage: L.series(2)
        O(5^4) + (3 + O(5))*T + (1 + O(5))*T^2 + (3 + O(5))*T^3 + O(5)*T^4 + O(T^5)

    Finally, we can use the overconvergent method of Pollack-Stevens.::

        sage: e = EllipticCurve('11a')
        sage: L = e.padic_lseries(3, implementation = 'pollackstevens', precision = 6)
        sage: L.series(5)
        2 + 3 + 3^2 + 2*3^3 + 2*3^5 + O(3^6) + (1 + 3 + 2*3^2 + 3^3 + O(3^4))*T + (1 + 2*3 + O(3^2))*T^2 + (3 + O(3^2))*T^3 + O(3^0)*T^4 + O(T^5)
        sage: L[3]
        3 + O(3^2)

    Another example with a semistable prime.::

        sage: E = EllipticCurve("11a1")
        sage: L = E.padic_lseries(11, implementation = 'pollackstevens', precision=3)
        sage: L[1]
        10 + 3*11 + O(11^2)
        sage: L[3]
        O(11^0)
    """
    p, normalize, implementation, precision = self._normalize_padic_lseries(p,
        normalize, implementation, precision)

    if implementation in ['sage', 'eclib', 'num']:
        if self.ap(p) % p != 0:
            Lp = plseries.pAdicLseriesOrdinary(self, p,
                                  normalize=normalize, implementation=implementation)
        else:
            Lp = plseries.pAdicLseriesSupersingular(self, p,
                                  normalize=normalize, implementation=implementation)
    else:
        phi = self.pollack_stevens_modular_symbol(sign=0)
        if phi.parent().level() % p == 0:
            Phi = phi.lift(p, precision, eigensymbol=True)
        else:
            Phi = phi.p_stabilize_and_lift(p, precision, eigensymbol=True)
        Lp = Phi.padic_lseries()  #mm TODO should this pass precision on too ?
        Lp._cinf = self.real_components()
    return Lp


def padic_regulator(self, p, prec=20, height=None, check_hypotheses=True):
    r"""
    Compute the cyclotomic `p`-adic regulator of this curve.
    The model of the curve needs to be integral and minimal at `p`.
    Moreover the reduction at `p` should not be additive.

    INPUT:

    - ``p`` -- prime >= 5

    - ``prec`` -- answer will be returned modulo `p^{\mathrm{prec}}`

    - ``height`` -- precomputed height function; if not supplied, this function
      will call ``padic_height`` to compute it

    - ``check_hypotheses`` -- boolean; whether to check
      that this is a curve for which the `p`-adic height makes sense

    OUTPUT: the `p`-adic cyclotomic regulator of this curve, to the
    requested precision

    If the rank is 0, we output 1.

    AUTHORS:

    - Liang Xiao: original implementation at the 2006 MSRI
      graduate workshop on modular forms

    - David Harvey (2006-09-13): cleaned up and integrated into Sage,
      removed some redundant height computations

    - Chris Wuthrich (2007-05-22): added multiplicative and
      supersingular cases

    - David Harvey (2007-09-20): fixed some precision loss that was
      occurring

    EXAMPLES::

        sage: E = EllipticCurve("37a")
        sage: E.padic_regulator(5, 10)
        5 + 5^2 + 5^3 + 3*5^6 + 4*5^7 + 5^9 + O(5^10)

    An anomalous case::

        sage: E.padic_regulator(53, 10)
        26*53^-1 + 30 + 20*53 + 47*53^2 + 10*53^3 + 32*53^4 + 9*53^5 + 22*53^6 + 35*53^7 + 30*53^8 + O(53^9)

    An anomalous case where the precision drops some::

        sage: E = EllipticCurve("5077a")
        sage: E.padic_regulator(5, 10)
        5 + 5^2 + 4*5^3 + 2*5^4 + 2*5^5 + 2*5^6 + 4*5^7 + 2*5^8 + 5^9 + O(5^10)

    Check that answers agree over a range of precisions::

        sage: max_prec = 30    # make sure we get past p^2    # long time
        sage: full = E.padic_regulator(5, max_prec)           # long time
        sage: for prec in range(1, max_prec):                 # long time
        ....:     assert E.padic_regulator(5, prec) == full

    A case where the generator belongs to the formal group already
    (:issue:`3632`)::

        sage: E = EllipticCurve([37,0])
        sage: E.padic_regulator(5,10)
        2*5^2 + 2*5^3 + 5^4 + 5^5 + 4*5^6 + 3*5^8 + 4*5^9 + O(5^10)

    The result is not dependent on the model for the curve::

        sage: E = EllipticCurve([0,0,0,0,2^12*17])
        sage: Em = E.minimal_model()
        sage: E.padic_regulator(7) == Em.padic_regulator(7)
        True

    Allow a python int as input::

        sage: E = EllipticCurve('37a')
        sage: E.padic_regulator(int(5))
        5 + 5^2 + 5^3 + 3*5^6 + 4*5^7 + 5^9 + 5^10 + 3*5^11 + 3*5^12 + 5^13 + 4*5^14 + 5^15 + 2*5^16 + 5^17 + 2*5^18 + 4*5^19 + O(5^20)
    """
    p = Integer(p)  # this is assumed in code below
    if check_hypotheses:
        if not p.is_prime():
            raise ValueError("p = (%s) must be prime" % p)
        if p == 2:
            raise ValueError("p must be odd")   # todo
        if self.conductor() % (p**2) == 0:
            raise ArithmeticError("p must be a semi-stable prime")

    if self.conductor() % p == 0:
        Eq = self.tate_curve(p)
        reg = Eq.padic_regulator(prec=prec)
        return reg
    elif self.ap(p) % p == 0:
        lp = self.padic_lseries(p)
        reg = lp.Dp_valued_regulator(prec=prec)
        return reg
    else:
        if self.rank() == 0:
            return Qp(p,prec)(1)
        if height is None:
            height = self.padic_height(p, prec, check_hypotheses=False)
        d = self.padic_height_pairing_matrix(p=p, prec=prec, height=height, check_hypotheses=False)
        return d.determinant()


def padic_height_pairing_matrix(self, p, prec=20, height=None, check_hypotheses=True):
    r"""
    Compute the cyclotomic `p`-adic height pairing matrix of
    this curve with respect to the basis ``self.gens()`` for the
    Mordell-Weil group for a given odd prime `p` of good ordinary
    reduction.
    The model needs to be integral and minimal at `p`.

    INPUT:

    - ``p`` -- prime >= 5

    - ``prec`` -- answer will be returned modulo `p^{\mathrm{prec}}`

    - ``height`` -- precomputed height function; if not supplied, this function
      will call ``padic_height`` to compute it

    - ``check_hypotheses`` -- boolean; whether to check that this is a curve
      for which the `p`-adic height makes sense

    OUTPUT: the `p`-adic cyclotomic height pairing matrix of this curve
    to the given precision

    AUTHORS:

    - David Harvey, Liang Xiao, Robert Bradshaw, Jennifer
      Balakrishnan: original implementation at the 2006 MSRI graduate
      workshop on modular forms

    - David Harvey (2006-09-13): cleaned up and integrated into Sage,
      removed some redundant height computations

    EXAMPLES::

        sage: E = EllipticCurve("37a")
        sage: E.padic_height_pairing_matrix(5, 10)
        [5 + 5^2 + 5^3 + 3*5^6 + 4*5^7 + 5^9 + O(5^10)]

    A rank two example::

        sage: e =EllipticCurve('389a')
        sage: e._set_gens([e(-1, 1), e(1,0)])  # avoid platform dependent gens
        sage: e.padic_height_pairing_matrix(5,10)
        [                      3*5 + 2*5^2 + 5^4 + 5^5 + 5^7 + 4*5^9 + O(5^10) 5 + 4*5^2 + 5^3 + 2*5^4 + 3*5^5 + 4*5^6 + 5^7 + 5^8 + 2*5^9 + O(5^10)]
        [5 + 4*5^2 + 5^3 + 2*5^4 + 3*5^5 + 4*5^6 + 5^7 + 5^8 + 2*5^9 + O(5^10)                         4*5 + 2*5^4 + 3*5^6 + 4*5^7 + 4*5^8 + O(5^10)]

    An anomalous rank 3 example::

        sage: e = EllipticCurve("5077a")
        sage: e._set_gens([e(-1,3), e(2,0), e(4,6)])
        sage: e.padic_height_pairing_matrix(5,4)
        [4 + 3*5 + 4*5^2 + 4*5^3 + O(5^4)       4 + 4*5^2 + 2*5^3 + O(5^4)       3*5 + 4*5^2 + 5^3 + O(5^4)]
        [      4 + 4*5^2 + 2*5^3 + O(5^4)   3 + 4*5 + 3*5^2 + 5^3 + O(5^4)                 2 + 4*5 + O(5^4)]
        [      3*5 + 4*5^2 + 5^3 + O(5^4)                 2 + 4*5 + O(5^4)     1 + 3*5 + 5^2 + 5^3 + O(5^4)]
    """
    if check_hypotheses:
        p = __check_padic_hypotheses(self, p)

    K = Qp(p, prec=prec)

    rank = self.rank()
    M = matrix(K, rank, rank, 0)
    if rank == 0:
        return M

    basis = self.gens()

    if height is None:
        height = self.padic_height(p, prec, check_hypotheses=False)

    # Use <P, Q> =1/2*( h(P + Q) - h(P) - h(Q) )

    for i in range(rank):
        M[i,i] = height(basis[i])
    for i in range(rank):
        for j in range(i+1, rank):
            M[i, j] = ( height(basis[i] + basis[j]) - M[i,i] - M[j,j] ) / 2
            M[j, i] = M[i, j]

    return M


def _multiply_point(E, R, P, m):
    r"""
    Compute coordinates of a multiple of `P` with entries in a ring.

    INPUT:

    - ``E`` -- elliptic curve over `\QQ` with integer
      coefficients

    - ``P`` -- a rational point on `P` that reduces to a
      non-singular point at all primes

    - ``R`` -- a ring in which 2 is invertible (typically
      `\ZZ/L\ZZ` for some positive odd integer `L`)

    - ``m`` -- integer (default: 1)

    OUTPUT:

    A triple `(a', b', d')` such that if the point
    `mP` has coordinates `(a/d^2, b/d^3)`, then we have
    `a' \equiv a`, `b' \equiv \pm b`,
    `d' \equiv \pm d` all in `R` (i.e. modulo
    `L`).

    Note the ambiguity of signs for `b'` and `d'`.
    There's not much one can do about this, but at least one can say
    that the sign for `b'` will match the sign for
    `d'`.

    ALGORITHM: Proposition 9 in [Har2009]_.

    Complexity is soft-`O(\log L \log m + \log^2 m)`.

    AUTHORS:

    - David Harvey (2008-01): replaced _DivPolyContext with
      _multiply_point

    EXAMPLES:

    37a has trivial Tamagawa numbers so all points have nonsingular
    reduction at all primes::

        sage: E = EllipticCurve("37a")
        sage: P = E([0, -1]); P
        (0 : -1 : 1)
        sage: 19*P
        (-59997896/67387681 : 88075171080/553185473329 : 1)
        sage: R = Integers(625)
        sage: from sage.schemes.elliptic_curves.padics import _multiply_point
        sage: _multiply_point(E, R, P, 19)
        (229, 170, 541)
        sage: -59997896 % 625
        229
        sage: -88075171080 % 625         # note sign is flipped
        170
        sage: -67387681.sqrt() % 625     # sign is flipped here too
        541

    Trivial cases (:issue:`3632`)::

        sage: _multiply_point(E, R, P, 1)
        (0, 624, 1)
        sage: _multiply_point(E, R, 19*P, 1)
        (229, 455, 84)
        sage: (170 + 455) % 625       # note the sign did not flip here
        0
        sage: (541 + 84) % 625
        0

    Test over a range of `n` for a single curve with fairly
    random coefficients::

        sage: R = Integers(625)
        sage: E = EllipticCurve([4, -11, 17, -8, -10])
        sage: P = E.gens()[0] * LCM(E.tamagawa_numbers())
        sage: from sage.schemes.elliptic_curves.padics import _multiply_point
        sage: Q = P
        sage: for n in range(1, 25):
        ....:      naive = R(Q[0].numerator()), \
        ....:              R(Q[1].numerator()), \
        ....:              R(Q[0].denominator().sqrt())
        ....:      triple = _multiply_point(E, R, P, n)
        ....:      assert (triple[0] == naive[0]) and ( \
        ....:        (triple[1] == naive[1] and triple[2] == naive[2]) or \
        ....:        (triple[1] == -naive[1] and triple[2] == -naive[2])), \
        ....:           "_multiply_point() gave an incorrect answer"
        ....:      Q = Q + P
    """
    assert m >= 1

    alpha = R(P[0].numerator())
    beta = R(P[1].numerator())
    d = R(P[0].denominator().sqrt())
    if m == 1:
        return alpha, beta, d

    a1 = R(E.a1()) * d
    a3 = R(E.a3()) * d**3

    b2 = R(E.b2()) * d**2
    b4 = R(E.b4()) * d**4
    b6 = R(E.b6()) * d**6
    b8 = R(E.b8()) * d**8

    B4 = 6*alpha**2 + b2*alpha + b4
    B6 = 4*alpha**3 + b2*alpha**2 + 2*b4*alpha + b6
    B6_sqr = B6*B6
    B8 = 3*alpha**4 + b2*alpha**3 + 3*b4*alpha**2 + 3*b6*alpha + b8

    T = 2*beta + a1*alpha + a3

    # make a list of disjoint intervals [a[i], b[i]) such that we need to
    # compute g(k) for all a[i] <= k <= b[i] for each i
    intervals = []
    interval = (m - 2, m + 3)
    while interval[0] < interval[1]:
        intervals.append(interval)
        interval = max((interval[0] - 3) >> 1, 0), \
                   min((interval[1] + 5) >> 1, interval[0])

    # now walk through list and compute g(k)
    g = {0 : R(0), 1 : R(1), 2 : R(-1), 3 : B8, 4 : B6**2 - B4*B8}
    for i in reversed(intervals):
        k = i[0]
        while k < i[1]:
            if k > 4:
                j = k >> 1
                if k & 1:
                    t1 = g[j]
                    t2 = g[j+1]
                    prod1 = g[j+2] * t1*t1*t1
                    prod2 = g[j-1] * t2*t2*t2
                    g[k] = prod1 - B6_sqr * prod2 if j & 1 else B6_sqr * prod1 - prod2
                else:
                    t1 = g[j-1]
                    t2 = g[j+1]
                    g[k] = g[j] * (g[j-2] * t2*t2 - g[j+2] * t1*t1)
            k = k + 1

    if m & 1:
        psi_m = g[m]
        psi_m_m1 = g[m-1] * T
        psi_m_p1 = g[m+1] * T
    else:
        psi_m = g[m] * T
        psi_m_m1 = g[m-1]
        psi_m_p1 = g[m+1]

    theta = alpha * psi_m * psi_m - psi_m_m1 * psi_m_p1
    t1 = g[m-2] * g[m+1] * g[m+1] - g[m+2] * g[m-1] * g[m-1]
    if m & 1:
        t1 = t1 * T
    omega = (t1 + (a1 * theta + a3 * psi_m * psi_m) * psi_m) / -2

    return theta, omega, psi_m * d


def _multiple_to_make_good_reduction(E):
    r"""
    Return the integer `n_2` such that for all points `P` in `E(\QQ)`
    `n_2*P` has good reduction at all primes.
    If the model is globally minimal the lcm of the
    Tamagawa numbers will do, otherwise we have to take into
    account the change of the model.

    INPUT:

    - ``E`` -- an elliptic curve over `\QQ`

    OUTPUT: a positive integer ``n2``

    EXAMPLE (:issue:`34790`)::

        sage: from sage.schemes.elliptic_curves.padics import _multiple_to_make_good_reduction
        sage: E = EllipticCurve([-1728,-100656])
        sage: _multiple_to_make_good_reduction(E)
        30

    The number ``n_2`` is not always optimal but it is in this example.
    The first multiple of the generator `P` with good reduction in this
    non-minimal model is `30 P`.

    TESTS::

        sage: from sage.schemes.elliptic_curves.padics import _multiple_to_make_good_reduction
        sage: E = EllipticCurve([1/2,1/3])
        sage: _multiple_to_make_good_reduction(E)
        Traceback (most recent call last):
        ...
        NotImplementedError: This only implemented for integral models. Please change the model first.
        sage: E = EllipticCurve([0,3])
        sage: _multiple_to_make_good_reduction(E)
        1
        sage: E = EllipticCurve([0,5^7]) # min eq is additive
        sage: _multiple_to_make_good_reduction(E)
        5
        sage: E = EllipticCurve([7,0,0,0,7^7]) # min eq is split mult
        sage: _multiple_to_make_good_reduction(E)
        6
        sage: E = EllipticCurve([0,-3^2,0,0,3^7]) # min eq is non-split mult
        sage: _multiple_to_make_good_reduction(E)
        4
    """
    if not E.is_integral():
        st = ("This only implemented for integral models. "
              "Please change the model first.")
        raise NotImplementedError(st)
    if E.is_minimal():
        n2 = LCM(E.tamagawa_numbers())
    else:
        # generalising to number fields one can get the u from local_data
        Emin = E.global_minimal_model()
        iota = E.isomorphism_to(Emin)
        u = Integer(iota.u)
        ps = u.prime_divisors()
        li = []
        for p in ps:
            np = u.valuation(p)
            if Emin.discriminant() % p != 0:
                li.append(Emin.Np(p) * p**(np-1))
            elif Emin.has_additive_reduction(p):
                li.append(E.tamagawa_number(p) * p**np)
            elif E.has_split_multiplicative_reduction(p):
                li.append(E.tamagawa_number(p) * (p-1) * p**(np-1))
            else: # non split
                li.append(E.tamagawa_number(p) * (p+1) * p**(np-1))
        otherbad = Integer(Emin.discriminant()).prime_divisors()
        otherbad = [p for p in otherbad if u % p != 0 ]
        li += [E.tamagawa_number(p) for p in otherbad]
        n2 = LCM(li)
    return n2


def padic_height(self, p, prec=20, sigma=None, check_hypotheses=True):
    r"""
    Compute the cyclotomic `p`-adic height.

    The equation of the curve must be integral and minimal at `p`.

    INPUT:

    - ``p`` -- prime >= 5 for which the curve has semi-stable reduction

    - ``prec`` -- integer >= 1 (default: 20); desired precision of result

    - ``sigma`` -- precomputed value of sigma; if not supplied, this function
      will call ``padic_sigma`` to compute it

    - ``check_hypotheses`` -- boolean; whether to check that this is a curve
      for which the `p`-adic height makes sense

    OUTPUT:

    A function that accepts two parameters:

    - a `\QQ`-rational point on the curve whose height should be computed

    - optional boolean flag 'check': if ``False``, it skips some input
      checking, and returns the `p`-adic height of that point to the
      desired precision.

    - The normalization (sign and a factor 1/2 with respect to some other
      normalizations that appear in the literature) is chosen in such a way
      as to make the `p`-adic Birch Swinnerton-Dyer conjecture hold as stated
      in [MTT1986]_.

    AUTHORS:

    - Jennifer Balakrishnan: original code developed at the 2006 MSRI
      graduate workshop on modular forms

    - David Harvey (2006-09-13): integrated into Sage, optimised to
      speed up repeated evaluations of the returned height function,
      addressed some thorny precision questions

    - David Harvey (2006-09-30): rewrote to use division polynomials
      for computing denominator of `nP`.

    - David Harvey (2007-02): cleaned up according to algorithms in
      "Efficient Computation of p-adic Heights"

    - Chris Wuthrich (2007-05): added supersingular and multiplicative heights

    EXAMPLES::

        sage: E = EllipticCurve("37a")
        sage: P = E.gens()[0]
        sage: h = E.padic_height(5, 10)
        sage: h(P)
        5 + 5^2 + 5^3 + 3*5^6 + 4*5^7 + 5^9 + O(5^10)

    An anomalous case::

        sage: h = E.padic_height(53, 10)
        sage: h(P)
        26*53^-1 + 30 + 20*53 + 47*53^2 + 10*53^3 + 32*53^4 + 9*53^5 + 22*53^6 + 35*53^7 + 30*53^8 + 17*53^9 + O(53^10)

    Boundary case::

        sage: E.padic_height(5, 3)(P)
        5 + 5^2 + O(5^3)

    A case that works the division polynomial code a little harder::

        sage: E.padic_height(5, 10)(5*P)
        5^3 + 5^4 + 5^5 + 3*5^8 + 4*5^9 + O(5^10)

    Check that answers agree over a range of precisions::

        sage: max_prec = 30    # make sure we get past p^2    # long time
        sage: full = E.padic_height(5, max_prec)(P)           # long time
        sage: for prec in range(1, max_prec):                 # long time
        ....:     assert E.padic_height(5, prec)(P) == full

    A supersingular prime for a curve::

        sage: E = EllipticCurve('37a')
        sage: E.is_supersingular(3)
        True
        sage: h = E.padic_height(3, 5)
        sage: h(E.gens()[0])
        (3 + 3^3 + O(3^6), 2*3^2 + 3^3 + 3^4 + 3^5 + 2*3^6 + O(3^7))
        sage: E.padic_regulator(5)
        5 + 5^2 + 5^3 + 3*5^6 + 4*5^7 + 5^9 + 5^10 + 3*5^11 + 3*5^12 + 5^13 + 4*5^14 + 5^15 + 2*5^16 + 5^17 + 2*5^18 + 4*5^19 + O(5^20)
        sage: E.padic_regulator(3, 5)
        (3 + 2*3^2 + 3^3 + O(3^4), 3^2 + 2*3^3 + 3^4 + O(3^5))

    A torsion point in both the good and supersingular cases::

        sage: E = EllipticCurve('11a')
        sage: P = E.torsion_subgroup().gen(0).element(); P
        (5 : 5 : 1)
        sage: h = E.padic_height(19, 5)
        sage: h(P)
        0
        sage: h = E.padic_height(5, 5)
        sage: h(P)
        0

    The result is not dependent on the model for the curve::

        sage: E = EllipticCurve([0,0,0,0,2^12*17])
        sage: Em = E.minimal_model()
        sage: P = E.gens()[0]
        sage: Pm = Em.gens()[0]
        sage: h = E.padic_height(7)
        sage: hm = Em.padic_height(7)
        sage: h(P) == hm(Pm)
        True

    TESTS:

    Check that issue :issue:`20798` is solved::

        sage: E = EllipticCurve("91b")
        sage: h = E.padic_height(7,10)
        sage: P = E.gen(0)
        sage: h(P)
        2*7 + 7^2 + 5*7^3 + 6*7^4 + 2*7^5 + 3*7^6 + 7^7 + 4*7^9 + 5*7^10 + O(7^11)
        sage: h(P+P)
        7 + 5*7^2 + 6*7^3 + 5*7^4 + 4*7^5 + 6*7^6 + 5*7^7 + 2*7^9 + 7^10 + O(7^11)
    """
    if check_hypotheses:
        if not p.is_prime():
            raise ValueError("p = (%s) must be prime" % p)
        if p == 2:
            raise ValueError("p must be odd")   # todo
        if self.conductor() % (p**2) == 0:
            raise ArithmeticError("p must be a semi-stable prime")

    prec = int(prec)
    if prec < 1:
        raise ValueError("prec (=%s) must be at least 1" % prec)

    if self.conductor() % p == 0:
        Eq = self.tate_curve(p)
        return Eq.padic_height(prec=prec)
    elif self.ap(p) % p == 0:
        lp = self.padic_lseries(p)
        return lp.Dp_valued_height(prec=prec)

    # else good ordinary case

    # For notation and definitions, see [Har2009]_.

    n1 = self.change_ring(GF(p)).cardinality()
    n2 = _multiple_to_make_good_reduction(self)
    n = LCM(n1, n2)
    m = int(n / n2)

    adjusted_prec = prec + 2 * valuation(n, p)   # this is M'
    R = Integers(p ** adjusted_prec)

    if sigma is None:
        sigma = self.padic_sigma(p, adjusted_prec, check_hypotheses=False)

    # K is the field for the final result
    K = Qp(p, prec=adjusted_prec-1)
    E = self

    def height(P, check=True):
        if P.is_finite_order():
            return K(0)

        if check:
            assert P.curve() == E, "the point P must lie on the curve " \
                   "from which the height function was created"

        Q = n2 * P
        alpha, beta, d = _multiply_point(E, R, Q, m)

        assert beta.lift() % p != 0, "beta should be a unit!"
        assert d.lift() % p == 0, "d should not be a unit!"

        t = -d * alpha / beta

        total = R(1)
        t_power = t
        for k in range(2, adjusted_prec + 1):
            total = total + t_power * sigma[k].lift()
            t_power = t_power * t
        total = (-alpha / beta) * total

        L = Qp(p, prec=adjusted_prec)
        total = L(total.lift(), adjusted_prec)   # yuck... get rid of this lift!

        # changed sign to make it correct for p-adic bsd
        answer = -total.log() * 2 / n**2

        if check:
            assert answer.precision_absolute() >= prec, "we should have got an " \
                   "answer with precision at least prec, but we didn't."
        return K(answer)

    # (man... I love python's local function definitions...)
    return height


def padic_height_via_multiply(self, p, prec=20, E2=None, check_hypotheses=True):
    r"""
    Compute the cyclotomic `p`-adic height.

    The equation of the curve must be minimal at `p`.

    INPUT:

    - ``p`` -- prime >= 5 for which the curve has good ordinary reduction

    - ``prec`` -- integer >= 2 (default: 20); desired precision of result

    - ``E2`` -- precomputed value of E2. If not supplied,
      this function will call padic_E2 to compute it. The value supplied
      must be correct mod `p^{prec-2}` (or slightly higher in the
      anomalous case; see the code for details).

    - ``check_hypotheses`` -- boolean; whether to check
      that this is a curve for which the `p`-adic height makes sense

    OUTPUT:

    A function that accepts two parameters:

    - a `\QQ`-rational point on the curve whose height should be computed

    - optional boolean flag 'check': if ``False``, it skips some input
      checking, and returns the `p`-adic height of that point to the
      desired precision.

    - The normalization (sign and a factor 1/2 with respect to some other
      normalizations that appear in the literature) is chosen in such a way
      as to make the `p`-adic Birch Swinnerton-Dyer conjecture hold as stated
      in [MTT1986]_.

    AUTHORS:

    - David Harvey (2008-01): based on the padic_height() function,
      using the algorithm of [Har2009]_.

    EXAMPLES::

        sage: E = EllipticCurve("37a")
        sage: P = E.gens()[0]
        sage: h = E.padic_height_via_multiply(5, 10)
        sage: h(P)
        5 + 5^2 + 5^3 + 3*5^6 + 4*5^7 + 5^9 + O(5^10)

    An anomalous case::

        sage: h = E.padic_height_via_multiply(53, 10)
        sage: h(P)
        26*53^-1 + 30 + 20*53 + 47*53^2 + 10*53^3 + 32*53^4 + 9*53^5 + 22*53^6 + 35*53^7 + 30*53^8 + 17*53^9 + O(53^10)

    Supply the value of E2 manually::

        sage: E2 = E.padic_E2(5, 8)
        sage: E2
        2 + 4*5 + 2*5^3 + 5^4 + 3*5^5 + 2*5^6 + O(5^8)
        sage: h = E.padic_height_via_multiply(5, 10, E2=E2)
        sage: h(P)
        5 + 5^2 + 5^3 + 3*5^6 + 4*5^7 + 5^9 + O(5^10)

    Boundary case::

        sage: E.padic_height_via_multiply(5, 3)(P)
        5 + 5^2 + O(5^3)

    Check that answers agree over a range of precisions::

        sage: max_prec = 30    # make sure we get past p^2    # long time
        sage: full = E.padic_height(5, max_prec)(P)           # long time
        sage: for prec in range(2, max_prec):                 # long time
        ....:     assert E.padic_height_via_multiply(5, prec)(P) == full
    """
    if check_hypotheses:
        if not p.is_prime():
            raise ValueError("p = (%s) must be prime" % p)
        if p == 2:
            raise ValueError("p must be odd")   # todo
        if self.conductor() % p == 0:
            raise ArithmeticError("must have good reduction at p")
        if self.ap(p) % p == 0:
            raise ArithmeticError("must be ordinary at p")

    prec = int(prec)
    if prec < 1:
        raise ValueError("prec (=%s) must be at least 1" % prec)

    # For notation and definitions, [Har2009]_

    n1 = self.change_ring(GF(p)).cardinality()
    n2 = _multiple_to_make_good_reduction(self)
    n = LCM(n1, n2)
    m = int(n / n2)

    lamb = int(math.floor(math.sqrt(prec)))

    adjusted_prec = prec + 2 * valuation(n, p)   # this is M'
    R = Integers(p ** (adjusted_prec + 2*lamb))

    sigma = self.padic_sigma_truncated(p, N=adjusted_prec, E2=E2, lamb=lamb)

    # K is the field for the final result
    K = Qp(p, prec=adjusted_prec-1)
    E = self

    def height(P, check=True):
        if P.is_finite_order():
            return K(0)

        if check:
            assert P.curve() == E, "the point P must lie on the curve " \
                   "from which the height function was created"

        Q = n2 * P
        alpha, beta, d = _multiply_point(E, R, Q, m * p**lamb)

        assert beta.lift() % p != 0, "beta should be a unit!"
        assert d.lift() % p == 0, "d should not be a unit!"

        t = -d * alpha / beta

        total = R(1)
        t_power = t
        for k in range(2, sigma.prec()):
            total = total + t_power * sigma[k].lift()
            t_power = t_power * t
        total = (-alpha / beta) * total

        L = Qp(p, prec=adjusted_prec + 2*lamb)
        total = L(total.lift(), adjusted_prec + 2*lamb)

        # changed sign to make it correct for p-adic bsd
        answer = -total.log() * 2 / (n * p**lamb)**2

        if check:
            assert answer.precision_absolute() >= prec, "we should have got an " \
                   "answer with precision at least prec, but we didn't."
        return K(answer)

    # (man... I love python's local function definitions...)
    return height


def padic_sigma(self, p, N=20, E2=None, check=False, check_hypotheses=True):
    r"""
    Compute the `p`-adic sigma function with respect to the standard
    invariant differential `dx/(2y + a_1 x + a_3)`, as
    defined by Mazur and Tate in [MT1991]_, as a power series in the usual
    uniformiser `t` at the origin.

    The equation of the curve must be minimal at `p`.

    INPUT:

    - ``p`` -- prime >= 5 for which the curve has good ordinary reduction

    - ``N`` -- integer >= 1 (default: 20); precision of result,
      see OUTPUT section for description

    - ``E2`` -- precomputed value of E2. If not supplied,
      this function will call padic_E2 to compute it. The value supplied
      must be correct mod `p^{N-2}`.

    - ``check`` -- boolean; whether to perform a
      consistency check (i.e. verify that the computed sigma satisfies
      the defining

    - ``differential equation`` -- note that this does NOT guarantee
      correctness of all the returned digits, but it comes pretty close

    - ``check_hypotheses`` -- boolean; whether to check that this is a curve
      for which the `p`-adic sigma function makes sense

    OUTPUT: a power series `t + \cdots` with coefficients in `\ZZ_p`

    The output series will be truncated at `O(t^{N+1})`, and
    the coefficient of `t^n` for `n \geq 1` will be
    correct to precision `O(p^{N-n+1})`.

    In practice this means the following. If `t_0 = p^k u`,
    where `u` is a `p`-adic unit with at least
    `N` digits of precision, and `k \geq 1`, then the
    returned series may be used to compute `\sigma(t_0)`
    correctly modulo `p^{N+k}` (i.e. with `N` correct
    `p`-adic digits).

    ALGORITHM: Described in "Efficient Computation of p-adic Heights"
    (David Harvey) [Har2009]_ which is basically an optimised version of the
    algorithm from "p-adic Heights and Log Convergence" (Mazur, Stein,
    Tate) [MST2006]_.

    Running time is soft-`O(N^2 \log p)`, plus whatever time is
    necessary to compute `E_2`.

    AUTHORS:

    - David Harvey (2006-09-12)

    - David Harvey (2007-02): rewrote

    EXAMPLES::

        sage: EllipticCurve([-1, 1/4]).padic_sigma(5, 10)
        O(5^11) + (1 + O(5^10))*t + O(5^9)*t^2 + (3 + 2*5^2 + 3*5^3 + 3*5^6 + 4*5^7 + O(5^8))*t^3 + O(5^7)*t^4 + (2 + 4*5^2 + 4*5^3 + 5^4 + 5^5 + O(5^6))*t^5 + O(5^5)*t^6 + (2 + 2*5 + 5^2 + 4*5^3 + O(5^4))*t^7 + O(5^3)*t^8 + (1 + 2*5 + O(5^2))*t^9 + O(5)*t^10 + O(t^11)

    Run it with a consistency check::

        sage: EllipticCurve("37a").padic_sigma(5, 10, check=True)
        O(5^11) + (1 + O(5^10))*t + O(5^9)*t^2 + (3 + 2*5^2 + 3*5^3 + 3*5^6 + 4*5^7 + O(5^8))*t^3 + (3 + 2*5 + 2*5^2 + 2*5^3 + 2*5^4 + 2*5^5 + 2*5^6 + O(5^7))*t^4 + (2 + 4*5^2 + 4*5^3 + 5^4 + 5^5 + O(5^6))*t^5 + (2 + 3*5 + 5^4 + O(5^5))*t^6 + (4 + 3*5 + 2*5^2 + O(5^4))*t^7 + (2 + 3*5 + 2*5^2 + O(5^3))*t^8 + (4*5 + O(5^2))*t^9 + (1 + O(5))*t^10 + O(t^11)

    Boundary cases::

        sage: EllipticCurve([1, 1, 1, 1, 1]).padic_sigma(5, 1)
         (1 + O(5))*t + O(t^2)
        sage: EllipticCurve([1, 1, 1, 1, 1]).padic_sigma(5, 2)
         (1 + O(5^2))*t + (3 + O(5))*t^2 + O(t^3)

    Supply your very own value of E2::

        sage: X = EllipticCurve("37a")
        sage: my_E2 = X.padic_E2(5, 8)
        sage: my_E2 = my_E2 + 5**5    # oops!!!
        sage: X.padic_sigma(5, 10, E2=my_E2)
        O(5^11) + (1 + O(5^10))*t + O(5^9)*t^2 + (3 + 2*5^2 + 3*5^3 + 4*5^5 + 2*5^6 + 3*5^7 + O(5^8))*t^3 + (3 + 2*5 + 2*5^2 + 2*5^3 + 2*5^4 + 2*5^5 + 2*5^6 + O(5^7))*t^4 + (2 + 4*5^2 + 4*5^3 + 5^4 + 3*5^5 + O(5^6))*t^5 + (2 + 3*5 + 5^4 + O(5^5))*t^6 + (4 + 3*5 + 2*5^2 + O(5^4))*t^7 + (2 + 3*5 + 2*5^2 + O(5^3))*t^8 + (4*5 + O(5^2))*t^9 + (1 + O(5))*t^10 + O(t^11)

    Check that sigma is "weight 1".

    ::

        sage: f = EllipticCurve([-1, 3]).padic_sigma(5, 10)
        sage: g = EllipticCurve([-1*(2**4), 3*(2**6)]).padic_sigma(5, 10)
        sage: t = f.parent().gen()
        sage: f(2*t)/2
        (1 + O(5^10))*t + (4 + 3*5 + 3*5^2 + 3*5^3 + 4*5^4 + 4*5^5 + 3*5^6 + 5^7 + O(5^8))*t^3 + (3 + 3*5^2 + 5^4 + 2*5^5 + O(5^6))*t^5 + (4 + 5 + 3*5^3 + O(5^4))*t^7 + (4 + 2*5 + O(5^2))*t^9 + O(5)*t^10 + O(t^11)
        sage: g
        O(5^11) + (1 + O(5^10))*t + O(5^9)*t^2 + (4 + 3*5 + 3*5^2 + 3*5^3 + 4*5^4 + 4*5^5 + 3*5^6 + 5^7 + O(5^8))*t^3 + O(5^7)*t^4 + (3 + 3*5^2 + 5^4 + 2*5^5 + O(5^6))*t^5 + O(5^5)*t^6 + (4 + 5 + 3*5^3 + O(5^4))*t^7 + O(5^3)*t^8 + (4 + 2*5 + O(5^2))*t^9 + O(5)*t^10 + O(t^11)
        sage: f(2*t)/2 -g
        O(t^11)

    Test that it returns consistent results over a range of precision::

        sage: # long time
        sage: max_N = 30   # get up to at least p^2
        sage: E = EllipticCurve([1, 1, 1, 1, 1])
        sage: p = 5
        sage: E2 = E.padic_E2(5, max_N)
        sage: max_sigma = E.padic_sigma(p, max_N, E2=E2)
        sage: for N in range(3, max_N):
        ....:    sigma = E.padic_sigma(p, N, E2=E2)
        ....:    assert sigma == max_sigma
    """
    if check_hypotheses:
        p = __check_padic_hypotheses(self, p)

    # todo: implement the p == 3 case
    # NOTE: If we ever implement p == 3, it's necessary to check over
    # the precision loss estimates (below) very carefully; I think it
    # may become necessary to compute E2 to an even higher precision.
    if p < 5:
        raise NotImplementedError("p (=%s) must be at least 5" % p)

    N = int(N)

    # a few special cases for small N
    if N < 1:
        raise ValueError("N (=%s) must be at least 1" % N)

    if N == 1:
        # return simply t + O(t^2)
        K = Qp(p, 2)
        return PowerSeriesRing(K, "t")([K(0), K(1, 1)], prec=2)

    if N == 2:
        # return t + a_1/2 t^2 + O(t^3)
        K = Qp(p, 3)
        return PowerSeriesRing(K, "t")([K(0), K(1, 2),
                                        K(self.a1()/2, 1)], prec=3)

    if self.discriminant().valuation(p) != 0:
        raise NotImplementedError("equation of curve must be minimal at p")

    if E2 is None:
        E2 = self.padic_E2(p, N-2, check_hypotheses=False)
    elif E2.precision_absolute() < N-2:
        raise ValueError("supplied E2 has insufficient precision")

    QQt = LaurentSeriesRing(RationalField(), "x")

    R = Integers(p**(N-2))
    X = self.change_ring(R)
    c = (X.a1()**2 + 4*X.a2() - R(E2)) / 12

    f = X.formal_group().differential(N+2)   # f = 1 + ... + O(t^{N+2})
    x = X.formal_group().x(N)                # x = t^{-2} + ... + O(t^N)

    Rt = x.parent()

    A = (x + c) * f
    # do integral over QQ, to avoid divisions by p
    A = Rt(QQt(A).integral())
    A = (-X.a1()/2 - A) * f

    # Convert to a power series and remove the -1/x term.
    # Also we artificially bump up the accuracy from N-2 to N-1 digits;
    # the constant term needs to be known to N-1 digits, so we compute
    # it directly
    assert A.valuation() == -1 and A[-1] == 1
    A = A - A.parent().gen() ** (-1)
    A = A.power_series().list()
    R = Integers(p**(N-1))
    A = [R(u) for u in A]
    A[0] = self.change_ring(R).a1()/2     # fix constant term
    A = PowerSeriesRing(R, "x")(A, len(A))

    theta = _brent(A, p, N)
    sigma = theta * theta.parent().gen()

    # Convert the answer to power series over p-adics; drop the precision
    # of the t^k coefficient to p^(N-k+1).
    # [Note: there are actually more digits available, but it's a bit
    # tricky to figure out exactly how many, and we only need p^(N-k+1)
    # for p-adic height purposes anyway]
    K = Qp(p, N + 1)

    sigma = sigma.padded_list(N+1)

    sigma[0] = K(0, N + 1)
    sigma[1] = K(1, N)
    for n in range(2, N+1):
        sigma[n] = K(sigma[n].lift(), N - n + 1)

    S = PowerSeriesRing(K, "t", N+1)
    sigma = S(sigma, N+1)

    # if requested, check that sigma satisfies the appropriate
    # differential equation
    if check:
        R = Integers(p**N)
        X = self.change_ring(R)
        x = X.formal_group().x(N+5)       # few extra terms for safety
        f = X.formal_group().differential(N+5)
        c = (X.a1()**2 + 4*X.a2() - R(E2)) / 12

        # convert sigma to be over Z/p^N
        s = f.parent()(sigma)
        sinv = s**(-1)
        finv = f**(-1)

        # apply differential equation
        temp = (s.derivative() * sinv * finv).derivative() * finv + c + x

        # coefficient of t^k in the result should be zero mod p^(N-k-2)
        for k in range(N-2):
            assert temp[k].lift().valuation(p) >= N - k - 2, \
                        "sigma correctness check failed!"

    return sigma


def padic_sigma_truncated(self, p, N=20, lamb=0, E2=None, check_hypotheses=True):
    r"""
    Compute the `p`-adic sigma function with respect to the standard
    invariant differential `dx/(2y + a_1 x + a_3)`, as
    defined by Mazur and Tate in [MT1991]_, as a power series in the usual
    uniformiser `t` at the origin.

    The equation of the curve must be minimal at `p`.

    This function differs from :func:`padic_sigma` in the precision profile
    of the returned power series; see OUTPUT below.

    INPUT:

    - ``p`` -- prime >= 5 for which the curve has good ordinary reduction

    - ``N`` -- integer >= 2 (default: 20); precision of result,
      see OUTPUT section for description

    - ``lamb`` -- integer >= 0; see OUTPUT section for description

    - ``E2`` -- precomputed value of E2. If not supplied,
      this function will call padic_E2 to compute it. The value supplied
      must be correct mod `p^{N-2}`.

    - ``check_hypotheses`` -- boolean; whether to check that this is a curve
      for which the `p`-adic sigma function makes sense

    OUTPUT: a power series `t + \cdots` with coefficients in `\ZZ_p`

    The coefficient of `t^j` for `j \geq 1` will be
    correct to precision `O(p^{N - 2 + (3 - j)(lamb + 1)})`.

    ALGORITHM: Described in "Efficient Computation of p-adic Heights"
    [Har2009]_, which is basically an
    optimised version of the algorithm from "p-adic Heights and Log
    Convergence" (Mazur, Stein, Tate) [MST2006]_.

    Running time is soft-`O(N^2 \lambda^{-1} \log p)`, plus
    whatever time is necessary to compute `E_2`.

    AUTHORS:

    - David Harvey (2008-01): wrote based on previous
      :func:`padic_sigma` function

    EXAMPLES::

        sage: E = EllipticCurve([-1, 1/4])
        sage: E.padic_sigma_truncated(5, 10)
        O(5^11) + (1 + O(5^10))*t + O(5^9)*t^2 + (3 + 2*5^2 + 3*5^3 + 3*5^6 + 4*5^7 + O(5^8))*t^3 + O(5^7)*t^4 + (2 + 4*5^2 + 4*5^3 + 5^4 + 5^5 + O(5^6))*t^5 + O(5^5)*t^6 + (2 + 2*5 + 5^2 + 4*5^3 + O(5^4))*t^7 + O(5^3)*t^8 + (1 + 2*5 + O(5^2))*t^9 + O(5)*t^10 + O(t^11)

    Note the precision of the `t^3` coefficient depends only on
    `N`, not on lamb::

        sage: E.padic_sigma_truncated(5, 10, lamb=2)
        O(5^17) + (1 + O(5^14))*t + O(5^11)*t^2 + (3 + 2*5^2 + 3*5^3 + 3*5^6 + 4*5^7 + O(5^8))*t^3 + O(5^5)*t^4 + (2 + O(5^2))*t^5 + O(t^6)

    Compare against plain padic_sigma() function over a dense range of
    N and lamb

    ::

        sage: E = EllipticCurve([1, 2, 3, 4, 7])                            # long time
        sage: E2 = E.padic_E2(5, 50)                                        # long time
        sage: for N in range(2, 10):                                        # long time
        ....:    for lamb in range(10):
        ....:       correct = E.padic_sigma(5, N + 3*lamb, E2=E2)
        ....:       compare = E.padic_sigma_truncated(5, N=N, lamb=lamb, E2=E2)
        ....:       assert compare == correct
    """
    if check_hypotheses:
        p = __check_padic_hypotheses(self, p)

    # todo: implement the p == 3 case
    # NOTE: If we ever implement p == 3, it's necessary to check over
    # the precision loss estimates (below) very carefully; I think it
    # may become necessary to compute E2 to an even higher precision.
    if p < 5:
        raise NotImplementedError("p (=%s) must be at least 5" % p)

    N = int(N)
    lamb = int(lamb)

    if lamb < 0:
        raise ValueError("lamb (=%s) must be at least 0" % lamb)

    # a few special cases for small N
    if N <= 1:
        raise ValueError("N (=%s) must be at least 2" % N)

    if N == 2:
        # return t + a_1/2 t^2 + O(t^3)
        K = Qp(p, 3*(lamb+1))
        return PowerSeriesRing(K, "t")([K(0), K(1, 2*(lamb+1)),
                                        K(self.a1()/2, lamb+1)], prec=3)

    if self.discriminant().valuation(p) != 0:
        raise NotImplementedError("equation of curve must be minimal at p")

    if E2 is None:
        E2 = self.padic_E2(p, N-2, check_hypotheses=False)
    elif E2.precision_absolute() < N-2:
        raise ValueError("supplied E2 has insufficient precision")

    # The main part of the algorithm is exactly the same as
    # for padic_sigma(), but we truncate all the series earlier.
    # Want the answer O(t^(trunc+1)) instead of O(t^(N+1)) like in padic_sigma().
    trunc = (Integer(N-2) / (lamb + 1)).ceil() + 2

    QQt = LaurentSeriesRing(RationalField(), "x")

    R = Integers(p**(N-2))
    X = self.change_ring(R)
    c = (X.a1()**2 + 4*X.a2() - R(E2)) / 12

    f = X.formal_group().differential(trunc+2)   # f = 1 + ... + O(t^{trunc+2})
    x = X.formal_group().x(trunc)                # x = t^{-2} + ... + O(t^trunc)

    Rt = x.parent()

    A = (x + c) * f
    # do integral over QQ, to avoid divisions by p
    A = Rt(QQt(A).integral())
    A = (-X.a1()/2 - A) * f

    # Convert to a power series and remove the -1/x term.
    # Also we artificially bump up the accuracy from N-2 to N-1+lamb digits;
    # the constant term needs to be known to N-1+lamb digits, so we compute
    # it directly
    assert A.valuation() == -1 and A[-1] == 1
    A = A - A.parent().gen() ** (-1)
    A = A.power_series().list()
    R = Integers(p**(N-1+lamb))
    A = [R(u) for u in A]
    A[0] = self.change_ring(R).a1()/2     # fix constant term
    A = PowerSeriesRing(R, "x")(A, len(A))

    theta = _brent(A, p, trunc)
    sigma = theta * theta.parent().gen()

    # Convert the answer to power series over p-adics; drop the precision
    # of the t^j coefficient to p^{N - 2 + (3 - j)(lamb + 1)}).
    K = Qp(p, N - 2 + 3*(lamb+1))

    sigma = sigma.padded_list(trunc+1)

    sigma[0] = K(0, N - 2 + 3*(lamb+1))
    sigma[1] = K(1, N - 2 + 2*(lamb+1))
    for j in range(2, trunc+1):
        sigma[j] = K(sigma[j].lift(), N - 2 + (3 - j)*(lamb+1))

    S = PowerSeriesRing(K, "t", trunc + 1)
    sigma = S(sigma, trunc+1)

    return sigma


def padic_E2(self, p, prec=20, check=False, check_hypotheses=True, algorithm='auto'):
    r"""
    Return the value of the `p`-adic modular form `E2`
    for `(E, \omega)` where `\omega` is the usual
    invariant differential `dx/(2y + a_1 x + a_3)`.

    INPUT:

    - ``p`` -- prime (= 5) for which `E` is good and ordinary

    - ``prec`` -- (relative) `p`-adic precision (= 1) for result

    - ``check`` -- boolean; whether to perform a consistency check. This will
      slow down the computation by a constant factor 2. (The consistency check
      is to compute the whole matrix of frobenius on Monsky-Washnitzer
      cohomology, and verify that its trace is correct to the specified
      precision. Otherwise, the trace is used to compute one column from the
      other one (possibly after a change of basis).)

    - ``check_hypotheses`` -- boolean; whether to check that this is a curve
      for which the `p`-adic sigma function makes sense

    - ``algorithm`` -- one of ``'standard'``, ``'sqrtp'``, or
      ``'auto'``. This selects which version of Kedlaya's algorithm is used.
      The ``'standard'`` one is the one described in Kedlaya's paper. The
      ``'sqrtp'`` one has better performance for large `p`, but only
      works when `p > 6N` (`N=` ``prec``). The ``'auto'`` option
      selects ``'sqrtp'`` whenever possible.

      Note that if the ``'sqrtp'`` algorithm is used, a consistency check
      will automatically be applied, regardless of the setting of the
      ``check`` flag.

    OUTPUT: `p`-adic number to precision ``prec``

    .. NOTE::

        If the discriminant of the curve has nonzero valuation at p,
        then the result will not be returned mod `p^\text{prec}`,
        but it still *will* have ``prec`` *digits* of precision.

    .. TODO::

        Once we have a better implementation of the "standard"
        algorithm, the algorithm selection strategy for "auto" needs to be
        revisited.

    AUTHORS:

    - David Harvey (2006-09-01): partly based on code written
      by Robert Bradshaw at the MSRI 2006 modular forms workshop

    ACKNOWLEDGMENT: - discussion with Eyal Goren that led to the trace
    trick.

    EXAMPLES: Here is the example discussed in the paper "Computation
    of p-adic Heights and Log Convergence" (Mazur, Stein, Tate) [MST2006]_::

        sage: EllipticCurve([-1, 1/4]).padic_E2(5)
        2 + 4*5 + 2*5^3 + 5^4 + 3*5^5 + 2*5^6 + 5^8 + 3*5^9 + 4*5^10 + 2*5^11 + 2*5^12 + 2*5^14 + 3*5^15 + 3*5^16 + 3*5^17 + 4*5^18 + 2*5^19 + O(5^20)

    Let's try to higher precision (this is the same answer the MAGMA
    implementation gives)::

        sage: EllipticCurve([-1, 1/4]).padic_E2(5, 100)
        2 + 4*5 + 2*5^3 + 5^4 + 3*5^5 + 2*5^6 + 5^8 + 3*5^9 + 4*5^10 + 2*5^11 + 2*5^12 + 2*5^14 + 3*5^15 + 3*5^16 + 3*5^17 + 4*5^18 + 2*5^19 + 4*5^20 + 5^21 + 4*5^22 + 2*5^23 + 3*5^24 + 3*5^26 + 2*5^27 + 3*5^28 + 2*5^30 + 5^31 + 4*5^33 + 3*5^34 + 4*5^35 + 5^36 + 4*5^37 + 4*5^38 + 3*5^39 + 4*5^41 + 2*5^42 + 3*5^43 + 2*5^44 + 2*5^48 + 3*5^49 + 4*5^50 + 2*5^51 + 5^52 + 4*5^53 + 4*5^54 + 3*5^55 + 2*5^56 + 3*5^57 + 4*5^58 + 4*5^59 + 5^60 + 3*5^61 + 5^62 + 4*5^63 + 5^65 + 3*5^66 + 2*5^67 + 5^69 + 2*5^70 + 3*5^71 + 3*5^72 + 5^74 + 5^75 + 5^76 + 3*5^77 + 4*5^78 + 4*5^79 + 2*5^80 + 3*5^81 + 5^82 + 5^83 + 4*5^84 + 3*5^85 + 2*5^86 + 3*5^87 + 5^88 + 2*5^89 + 4*5^90 + 4*5^92 + 3*5^93 + 4*5^94 + 3*5^95 + 2*5^96 + 4*5^97 + 4*5^98 + 2*5^99 + O(5^100)

    Check it works at low precision too::

        sage: EllipticCurve([-1, 1/4]).padic_E2(5, 1)
        2 + O(5)
        sage: EllipticCurve([-1, 1/4]).padic_E2(5, 2)
        2 + 4*5 + O(5^2)
        sage: EllipticCurve([-1, 1/4]).padic_E2(5, 3)
        2 + 4*5 + O(5^3)

    TODO: With the old(-er), i.e., = sage-2.4 `p`-adics we got
    `5 + O(5^2)` as output, i.e., relative precision 1, but
    with the newer `p`-adics we get relative precision 0 and absolute
    precision 1.

    ::

        sage: EllipticCurve([1, 1, 1, 1, 1]).padic_E2(5, 1)
        O(5)

    Check it works for different models of the same curve (37a), even
    when the discriminant changes by a power of p (note that E2 depends
    on the differential too, which is why it gets scaled in some of the
    examples below)::

        sage: X1 = EllipticCurve([-1, 1/4])
        sage: X1.j_invariant(), X1.discriminant()
         (110592/37, 37)
        sage: X1.padic_E2(5, 10)
         2 + 4*5 + 2*5^3 + 5^4 + 3*5^5 + 2*5^6 + 5^8 + 3*5^9 + O(5^10)

    ::

        sage: X2 = EllipticCurve([0, 0, 1, -1, 0])
        sage: X2.j_invariant(), X2.discriminant()
         (110592/37, 37)
        sage: X2.padic_E2(5, 10)
         2 + 4*5 + 2*5^3 + 5^4 + 3*5^5 + 2*5^6 + 5^8 + 3*5^9 + O(5^10)

    ::

        sage: X3 = EllipticCurve([-1*(2**4), 1/4*(2**6)])
        sage: X3.j_invariant(), X3.discriminant() / 2**12
         (110592/37, 37)
        sage: 2**(-2) * X3.padic_E2(5, 10)
         2 + 4*5 + 2*5^3 + 5^4 + 3*5^5 + 2*5^6 + 5^8 + 3*5^9 + O(5^10)

    ::

        sage: X4 = EllipticCurve([-1*(7**4), 1/4*(7**6)])
        sage: X4.j_invariant(), X4.discriminant() / 7**12
         (110592/37, 37)
        sage: 7**(-2) * X4.padic_E2(5, 10)
         2 + 4*5 + 2*5^3 + 5^4 + 3*5^5 + 2*5^6 + 5^8 + 3*5^9 + O(5^10)

    ::

        sage: X5 = EllipticCurve([-1*(5**4), 1/4*(5**6)])
        sage: X5.j_invariant(), X5.discriminant() / 5**12
         (110592/37, 37)
        sage: 5**(-2) * X5.padic_E2(5, 10)
         2 + 4*5 + 2*5^3 + 5^4 + 3*5^5 + 2*5^6 + 5^8 + 3*5^9 + O(5^10)

    ::

        sage: X6 = EllipticCurve([-1/(5**4), 1/4/(5**6)])
        sage: X6.j_invariant(), X6.discriminant() * 5**12
         (110592/37, 37)
        sage: 5**2 * X6.padic_E2(5, 10)
         2 + 4*5 + 2*5^3 + 5^4 + 3*5^5 + 2*5^6 + 5^8 + 3*5^9 + O(5^10)

    Test check=True vs check=False::

        sage: EllipticCurve([-1, 1/4]).padic_E2(5, 1, check=False)
        2 + O(5)
        sage: EllipticCurve([-1, 1/4]).padic_E2(5, 1, check=True)
        2 + O(5)
        sage: EllipticCurve([-1, 1/4]).padic_E2(5, 30, check=False)
        2 + 4*5 + 2*5^3 + 5^4 + 3*5^5 + 2*5^6 + 5^8 + 3*5^9 + 4*5^10 + 2*5^11 + 2*5^12 + 2*5^14 + 3*5^15 + 3*5^16 + 3*5^17 + 4*5^18 + 2*5^19 + 4*5^20 + 5^21 + 4*5^22 + 2*5^23 + 3*5^24 + 3*5^26 + 2*5^27 + 3*5^28 + O(5^30)
        sage: EllipticCurve([-1, 1/4]).padic_E2(5, 30, check=True)
        2 + 4*5 + 2*5^3 + 5^4 + 3*5^5 + 2*5^6 + 5^8 + 3*5^9 + 4*5^10 + 2*5^11 + 2*5^12 + 2*5^14 + 3*5^15 + 3*5^16 + 3*5^17 + 4*5^18 + 2*5^19 + 4*5^20 + 5^21 + 4*5^22 + 2*5^23 + 3*5^24 + 3*5^26 + 2*5^27 + 3*5^28 + O(5^30)

    Here's one using the `p^{1/2}` algorithm::

        sage: EllipticCurve([-1, 1/4]).padic_E2(3001, 3, algorithm='sqrtp')
        1907 + 2819*3001 + 1124*3001^2 + O(3001^3)
    """
    if self.conductor() % p == 0:
        if not self.conductor() % (p**2) == 0:
            eq = self.tate_curve(p)
            return eq.E2(prec=prec)

    X = self.minimal_model().short_weierstrass_model()
    frob_p = X.matrix_of_frobenius(p, prec, check, check_hypotheses, algorithm).change_ring(Integers(p**prec))

    frob_p_n = frob_p**prec

    # todo: think about the sign of this. Is it correct?
    output_ring = Qp(p, prec)

    E2_of_X = output_ring( (-12 * frob_p_n[0,1] / frob_p_n[1,1]).lift() ) \
              + O(p**prec)

    # Take into account the coordinate change.
    fudge_factor = (X.discriminant() / self.discriminant()).nth_root(6)
    # todo: here I should be able to write:
    #  return E2_of_X / fudge_factor
    # However, there is a bug in Sage (#51 on trac) which makes this
    # crash sometimes when prec == 1. For example,
    #    EllipticCurve([1, 1, 1, 1, 1]).padic_E2(5, 1)
    # makes it crash. I haven't figured out exactly what the bug
    # is yet, but for now I use the following workaround:
    fudge_factor_inverse = Qp(p, prec=(E2_of_X.precision_absolute() + 1))(1 / fudge_factor)
    return output_ring(E2_of_X * fudge_factor_inverse)


def matrix_of_frobenius(self, p, prec=20, check=False, check_hypotheses=True, algorithm='auto'):
    r"""
    Return the matrix of Frobenius on the Monsky Washnitzer cohomology of
    the short Weierstrass model of the minimal model of the elliptic curve.

    INPUT:

    - ``p`` -- prime (>= 3) for which `E` is good and ordinary

    - ``prec`` -- (relative) `p`-adic precision for result (default: 20)

    - ``check`` -- boolean (default: ``False``); whether to perform a
      consistency check. This will slow down the computation by a
      constant factor 2. (The consistency check is to verify
      that its trace is correct to the specified precision. Otherwise,
      the trace is used to compute one column from the other one
      (possibly after a change of basis).)

    - ``check_hypotheses`` -- boolean; whether to check that this is a curve
      for which the `p`-adic sigma function makes sense

    - ``algorithm`` -- one of ``'standard'``, ``'sqrtp'``, or
      ``'auto'``. This selects which version of Kedlaya's algorithm is used.
      The ``'standard'`` one is the one described in Kedlaya's paper. The
      ``'sqrtp'`` one has better performance for large `p`, but only
      works when `p > 6N` (`N=` prec). The ``'auto'`` option
      selects ``'sqrtp'`` whenever possible.

      Note that if the ``'sqrtp'`` algorithm is used, a consistency check
      will automatically be applied, regardless of the setting of the
      ``check`` flag.

    OUTPUT: a matrix of `p`-adic number to precision ``prec``

    See also the documentation of padic_E2.

    EXAMPLES::

        sage: E = EllipticCurve('37a1')
        sage: E.matrix_of_frobenius(7)
        [             2*7 + 4*7^2 + 5*7^4 + 6*7^5 + 6*7^6 + 7^8 + 4*7^9 + 3*7^10 + 2*7^11 + 5*7^12 + 4*7^14 + 7^16 + 2*7^17 + 3*7^18 + 4*7^19 + 3*7^20 + O(7^21)                                   2 + 3*7 + 6*7^2 + 7^3 + 3*7^4 + 5*7^5 + 3*7^7 + 7^8 + 3*7^9 + 6*7^13 + 7^14 + 7^16 + 5*7^17 + 4*7^18 + 7^19 + O(7^20)]
        [    2*7 + 3*7^2 + 7^3 + 3*7^4 + 6*7^5 + 2*7^6 + 3*7^7 + 5*7^8 + 3*7^9 + 2*7^11 + 6*7^12 + 5*7^13 + 4*7^16 + 4*7^17 + 6*7^18 + 6*7^19 + 4*7^20 + O(7^21) 6 + 4*7 + 2*7^2 + 6*7^3 + 7^4 + 6*7^7 + 5*7^8 + 2*7^9 + 3*7^10 + 4*7^11 + 7^12 + 6*7^13 + 2*7^14 + 6*7^15 + 5*7^16 + 4*7^17 + 3*7^18 + 2*7^19 + O(7^20)]
        sage: M = E.matrix_of_frobenius(11,prec=3); M
        [   9*11 + 9*11^3 + O(11^4)          10 + 11 + O(11^3)]
        [     2*11 + 11^2 + O(11^4) 6 + 11 + 10*11^2 + O(11^3)]
        sage: M.det()
        11 + O(11^4)
        sage: M.trace()
        6 + 10*11 + 10*11^2 + O(11^3)
        sage: E.ap(11)
        -5
        sage: E = EllipticCurve('83a1')
        sage: E.matrix_of_frobenius(3,6)
        [                      2*3 + 3^5 + O(3^6)             2*3 + 2*3^2 + 2*3^3 + O(3^6)]
        [              2*3 + 3^2 + 2*3^5 + O(3^6) 2 + 2*3^2 + 2*3^3 + 2*3^4 + 3^5 + O(3^6)]
    """
    p = Integer(p)
    prec = int(prec)

    if p < 3:
        raise NotImplementedError("p (=%s) must be at least 3" % p)
    if prec < 1:
        raise ValueError("prec (=%s) must be at least 1" % prec)

    if check_hypotheses:
        p = __check_padic_hypotheses(self, p)

    if algorithm == "auto":
        algorithm = "standard" if p < 6*prec else "sqrtp"
    elif algorithm == "sqrtp" and p < 6*prec:
        raise ValueError("sqrtp algorithm is only available when p > 6*prec")

    if algorithm not in ["standard", "sqrtp"]:
        raise ValueError("unknown algorithm '%s'" % algorithm)

    # for p = 3, we create the corresponding hyperelliptic curve
    # and call matrix of frobenius on it
    if p == 3:
        from sage.schemes.hyperelliptic_curves.constructor import HyperellipticCurve
        f,g = self.hyperelliptic_polynomials()
        return HyperellipticCurve(f + (g/2)**2).matrix_of_frobenius(p,prec)

    # To run matrix_of_frobenius(), we need to have the equation in the
    # form y^2 = x^3 + ax + b, whose discriminant is invertible mod p.
    # When we change coordinates like this, we might scale the invariant
    # differential, so we need to account for this. We do this by
    # comparing discriminants: if the discriminants differ by u^12,
    # then the differentials differ by u. There's a sign ambiguity here,
    # but it doesn't matter because E2 changes by u^2 :-)

    # todo: In fact, there should be available a function that returns
    # exactly *which* coordinate change was used. If I had that I could
    # avoid the discriminant circus at the end.

    # todo: The following strategy won't work at all for p = 2, 3.

    # TODO change the basis back to the original equation.
    X = self.minimal_model().short_weierstrass_model()

    assert X.discriminant().valuation(p) == 0, "Something's gone wrong. " \
           "The discriminant of the Weierstrass model should be a unit " \
           " at p."

    if algorithm == "standard":
        # Need to increase precision a little to compensate for precision
        # losses during the computation. (See monsky_washnitzer.py
        # for more details.)
        adjusted_prec = sage.schemes.hyperelliptic_curves.monsky_washnitzer.adjusted_prec(p, prec)

        if check:
            trace = None
        else:
            trace = self.ap(p)

        base_ring = Integers(p**adjusted_prec)

        R, x = PolynomialRing(base_ring, 'x').objgen()
        Q = x**3 + base_ring(X.a4()) * x + base_ring(X.a6())
        frob_p = sage.schemes.hyperelliptic_curves.monsky_washnitzer.matrix_of_frobenius(
                         Q, p, adjusted_prec, trace)

    else:   # algorithm == "sqrtp"
        p_to_prec = p**prec
        R = PolynomialRing(Integers(), "x")
        Q = R([X.a6() % p_to_prec, X.a4() % p_to_prec, 0, 1])
        frob_p = sage.schemes.hyperelliptic_curves.hypellfrob.hypellfrob(p, prec, Q)

        # let's force a trace-check since this algorithm is fairly new
        # and we don't want to get caught with our pants down...
        trace = self.ap(p)
        check = True

    if check:
        trace_of_frobenius = frob_p.trace().lift() % p**prec
        correct_trace = self.ap(p) % p**prec
        assert trace_of_frobenius == correct_trace, \
                "Consistency check failed! (correct = %s, actual = %s)" % \
                (correct_trace, trace_of_frobenius)

    return frob_p.change_ring(Zp(p, prec))


def _brent(F, p, N):
    r"""
    This is an internal function; it is used by padic_sigma().

    `F` is a assumed to be a power series over `R = \ZZ/p^{N-1}\ZZ`.

    It solves the differential equation `G'(t)/G(t) = F(t)`
    using Brent's algorithm, with initial condition `G(0) = 1`.
    It is assumed that the solution `G` has `p`-integral coefficients.

    More precisely, suppose that `f(t)` is a power series with
    genuine `p`-adic coefficients, and suppose that
    `g(t)` is an exact solution to `g'(t)/g(t) = f(t)`.
    Let `I` be the ideal
    `(p^N, p^{N-1} t, \ldots, p t^{N-1}, t^N)`. The input
    `F(t)` should be a finite-precision approximation to
    `f(t)`, in the sense that `\int (F - f) dt` should
    lie in `I`. Then the function returns a series
    `G(t)` such that `(G - g)(t)` lies in `I`.

    Complexity should be about `O(N^2 \log^2 N \log p)`, plus
    some log-log factors.

    For more information, and a proof of the precision guarantees, see
    Lemma 4 in [Har2009]_.

    AUTHORS:

    - David Harvey (2007-02)

    EXAMPLES: Carefully test the precision guarantees::

        sage: brent = sage.schemes.elliptic_curves.padics._brent
        sage: for p in [2, 3, 5]:
        ....:   for N in [2, 3, 10, 50]:
        ....:     R = Integers(p**(N-1))
        ....:     Rx.<x> = PowerSeriesRing(R, "x")
        ....:     for _ in range(5):
        ....:       g = [R.random_element() for i in range(N)]
        ....:       g[0] = R(1)
        ....:       g = Rx(g, len(g))
        ....:       f = g.derivative() / g
        ....:       # perturb f by something whose integral is in I
        ....:       err = [R.random_element() * p**(N-i) for i in range(N+1)]
        ....:       err = Rx(err, len(err))
        ....:       err = err.derivative()
        ....:       F = f + err
        ....:       # run brent() and compare output modulo I
        ....:       G = brent(F, p, N)
        ....:       assert G.prec() >= N, "not enough output terms"
        ....:       err = (G - g).list()
        ....:       for i in range(len(err)):
        ....:         assert err[i].lift().valuation(p) >= (N - i), \
        ....:                "incorrect precision output"
    """
    Rx = F.parent()           # Rx = power series ring over Z/p^{N-1} Z
    Qx = PowerSeriesRing(RationalField(), "x")

    # initial approximation:
    G = Rx.one()

    # loop over an appropriate increasing sequence of lengths s
    for s in newton_method_sizes(N):
        # zero-extend to s terms
        # todo: there has to be a better way in Sage to do this...
        G = Rx(G.list(), s)

        # extend current approximation to be correct to s terms
        H = G.derivative() / G - F
        # Do the integral of H over QQ[x] to avoid division by p problems
        H = Rx(Qx(H).integral())
        G = G * (1 - H)

    return G
