r"""
General extensions of p-adic rings and fields; the base ring may also be an
extension.

These are implemented as proxy parents, backed by an absolute extension.

EXAMPLES:

A trivial extension::

    sage: L.<a> = Qp(2).extension(x)
    sage: L
    2-adic Trivial Extension Field in a defined by x
    sage: a == 0
    True

A trivial extension of a trivial extension::

    sage: R.<b> = L[]
    sage: M.<b> = L.extension(b - a)
    sage: M
    2-adic Trivial Extension Field in b defined by b
    sage: b == a
    True

An unramified extension::

    sage: L.<a> = Qp(2).extension(x^2 + 2*x + 4)
    sage: L
    2-adic Unramified Extension Field in a defined by x^2 + 2*x + 4
    sage: a^2 + 2*a + 4 == 0
    True
    sage: L.absolute_f()
    2

An unramified extension given by a non-monic defining polynomial (currently, not supported, see #33362)::

    sage: L.<a> = Qp(2).extension(4*x^2 + 2*x + 1)
    Traceback (most recent call last):
    ...
    ValueError: G must be integral
    sage: a^2 + 2*a + 4 == 0  # not tested
    True
    sage: L.absolute_f()
    2

An unramified extension given by a non-integral defining polynomial (currently, not supported, see #33362)::

    sage: L.<a> = Qp(2).extension(x^2 + x/4 + 1/16)
    Traceback (most recent call last):
    ...
    ValueError: G must be integral
    sage: a^2 + 2*a + 4 == 0  # not tested
    True
    sage: L.absolute_f()
    2

A trivial extension of an unramified extension::

    sage: L.<a> = Qp(2).extension(x^2 + 2*x + 4)
    sage: R.<b> = L[]
    sage: M.<b> = L.extension(b - 2)
    sage: M
    2-adic Trivial Extension Field in b defined by b - 2 over its base field
    sage: M.relative_f()
    1
    sage: M.absolute_f()
    2

An unramified extension of a trivial extension::

    sage: L.<a> = Qp(2).extension(x + 1)
    sage: R.<b> = L[]
    sage: M.<b> = L.extension(b^2 - b - a)
    sage: M
    2-adic Unramified Extension Field in b defined by b^2 - b + 1
    sage: M.absolute_f()
    2

An unramified extension of an unramified extension::

    sage: L.<a> = Qp(2).extension(x^2 + x + 1)
    sage: R.<b> = L[]
    sage: M.<b> = L.extension(b^2 + b + a)
    sage: M
    2-adic Unramified Extension Field in b defined by b^2 + b + a over its base field
    sage: M.relative_f()
    2
    sage: M.absolute_f()
    4

::

    sage: L.<a> = Qp(2).extension(x^2 + 2*x + 4)
    sage: R.<b> = L[]
    sage: M.<b> = L.extension(b^2 + a*b + 4)
    sage: M
    2-adic Unramified Extension Field in b defined by b^2 + a*b + 4 over its base field
    sage: M.relative_f()
    2
    sage: M.absolute_f()
    4

A totally ramified extension not given by an Eisenstein polynomial::

    sage: L.<a> = Qp(2).extension(x^2 + 8)
    sage: L.absolute_e()
    2

A trivial extension of a totally ramified extension::

    sage: R.<b> = L[]
    sage: M.<b> = L.extension(b - a)

A totally ramified extension of a trivial extension::

    sage: L.<a> = Qp(2).extension(x - 2)
    sage: R.<b> = L[]
    sage: M.<b> = L.extension(b^2 - a)
    sage: M.absolute_e(), M.absolute_f()
    (2, 1)

A totally ramified extension of an unramified extension::

    sage: L.<a> = Qp(2).extension(x^2 + 2*x + 4)
    sage: R.<b> = L[]
    sage: M.<b> = L.extension(b^2 - 8)  # long time, 1s in early 2022
    sage: M.absolute_e(), M.absolute_f()  # long time
    (2, 2)

An unramified extension of a totally ramified extension::

    sage: L.<a> = Qp(2).extension(x^2 + 8)
    sage: R.<b> = L[]
    sage: M.<b> = L.extension(b^2 + a*b + a^2)
    sage: M.absolute_e(), M.absolute_f()
    (2, 2)

A totally ramified extension of a totally ramified extension::

    sage: L.<a> = Qp(2).extension(x^2 + 8)
    sage: R.<b> = L[]
    sage: M.<b> = L.extension(b^2 + 2*a)
    sage: M.absolute_e(), M.absolute_f()
    (4, 1)

A mixed case::

    sage: L.<a> = Qp(2).extension(x^4 + 8*x^2 + 64)
    sage: L.absolute_e(), L.absolute_f()
    (2, 2)

A trivial extension of a mixed extension::

    sage: L.<a> = Qp(2).extension(x^4 + 8*x^2 + 64)
    sage: R.<b> = L[]
    sage: M.<b> = L.extension(b - a)
    sage: M.absolute_e(), M.absolute_f()
    (2, 2)

An unramified extension of a mixed extension::

    sage: L.<a> = Qp(2).extension(x^4 + 8*x^2 + 64)
    sage: R.<b> = L[]
    sage: M.<b> = L.extension(b^2 - a^2/8*b - a^2/8)  # long time, 5s in early 2022
    sage: M.absolute_e(), M.absolute_f()  # long time
    (2, 4)

An Eisenstein extension of a mixed extension::

    sage: L.<a> = Qp(2).extension(x^4 + 8*x^2 + 64)
    sage: R.<b> = L[]
    sage: M.<b> = L.extension(b^2 - L.uniformizer())
    sage: M.absolute_e(), M.absolute_f()
    (4, 2)

A mixed extension of a trivial extension::

    sage: L.<a> = Qp(2).extension(x - 2)
    sage: R.<b> = L[]
    sage: M.<b> = L.extension(x^4 + 8*x^2 + 64)
    sage: M.absolute_e(), M.absolute_f()
    (2, 2)

A mixed extension of an unramified extension::

    sage: L.<a> = Qp(2).extension(x^2 + x + 1)
    sage: R.<b> = L[]
    sage: M.<b> = L.extension(b^4 + 2*b^2 + 4*a)
    sage: M.absolute_e(), M.absolute_f()
    (2, 4)

::

    sage: L.<a> = Qp(2).extension(x^2 + 2*x + 4)
    sage: R.<b> = L[]
    sage: M.<b> = L.extension(b^4 - a*b^2 - 2*a)
    sage: M.absolute_e(), M.absolute_f()
    (2, 4)

A mixed extension of an Eisenstein extension::

    sage: L.<a> = Qp(2).extension(x^3 - 2)
    sage: R.<b> = L[]
    sage: M.<b> = L.extension(x^4 + 8*x^2 + 64)
    sage: M.absolute_e(), M.absolute_f()
    (6, 2)

A mixed extension of a mixed extension::

    TODO

A tower of mixed extensions::

    TODO

"""
# ****************************************************************************
#       Copyright (C)      2019 David Roe <roed.math@gmail.com>
#                     2019-2022 Julian Rüth <julian.rueth@fsfe.org>
#                     2019-2022 Xavier Caruso <xavier.caruso@normalesup.org>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

from .padic_general_extension_element import pAdicGeneralRingExtensionElement, pAdicGeneralFieldExtensionElement
from .padic_extension_generic import pAdicExtensionGeneric
from sage.rings.ring_extension import RingExtensionWithGen, RingExtension_generic
from sage.rings.ring_extension_conversion import backend_parent
from sage.rings.padics.pow_computer import PowComputer_class
from sage.rings.morphism import RingMap
from sage.categories.homset import Hom

from sage.rings.padics.precision_error import PrecisionError


# Helper functions
##################

def prem(P, Q):
    r"""
    Return the pseudo-remainder of `P` by `Q`.

    EXAMPLES::

        sage: from sage.rings.padics.padic_general_extension import prem
        sage: S.<x> = ZZ[]
        sage: P = 2*x^4 + 3*x^2 - x + 7
        sage: Q = 5*x^2 - 2
        sage: prem(P, Q)
        -125*x + 1065
    """
    lc = Q.leading_coefficient()
    p = P.degree()
    q = Q.degree()
    for i in range(p, q-1, -1):
        P = lc*P - (Q*P[i] << (i-q))
        P = P[:i]
    return P


def subresultant(P, Q, s=0):
    r"""
    Return the `s`-th subresultant of `P` and `Q` up to a sign.

    INPUT:

    - `P` and `Q` -- two univariate `p`-adic polynomials

    - `s` -- an integer (default: `0`)

    ALGORITHM:

    We use the subresultant PRS algorithm (see for instance [Bro1978]_).

    Besides, we lift all intermediate results to maximal precision in
    order to avoid numerical instability.
    This can lead to mathematically wrong results but the inaccuracy
    remains under control (see [Car2017]_) and we actually do not care here
    about precision because the result will be (truncated and) validated
    afterwards.

    EXAMPLES::

        sage: from sage.rings.padics.padic_general_extension import subresultant
        sage: S.<x> = Qp(5, print_mode='digits')[]
        sage: P = -x^5 - 2*x^4 - 12*x^2 + 1/5*x - 17
        sage: Q = -5*x^5 + 2*x^4 + x - 1
        sage: subresultant(P, Q)  # here s=0, so it's the resultant
        ...0012244212301202.1231

    We observe on this example that some digits at high precision are wrong::

        sage: P.resultant(Q)
        ...4412244212301202.1231

    ::

        sage: P1 = (x + 1) * P
        sage: Q1 = (x + 1) * Q
        sage: R = subresultant(P1, Q1, s=1)
        sage: R.monic()
        ...00000000000000000001*x + ...12100000000000000001

    """
    K = P.base_ring()
    if P.degree() > Q.degree():
        Gi, Gj = P, Q
    else:
        Gi, Gj = Q, P
    first = True
    while True:
        if first:
            h = K.one()
        else:
            g = Gi.leading_coefficient()
            h = (g**d / h**(d-1)).lift_to_precision()
        if Gj.degree() < s:
            break
        d = Gi.degree() - Gj.degree()
        if first:
            R = prem(Gi, Gj).map_coefficients(lambda x: x.lift_to_precision())
        else:
            c = g * h**d
            R = prem(Gi, Gj).map_coefficients(lambda x: (x / c).lift_to_precision())
        Gi, Gj = Gj, R
        first = False
    if Gi.degree() > s:
        return K.zero()
    else:
        return Gi


def resultant_bivariate(P, Q):
    r"""
    Compute the resultant of `P` and `Q` up to a sign.

    INPUT:

    - ``P`` and ``Q`` -- two bivariate polynomials in
      `K[y][x]` over a two-step extension `K` of `\QQ_p`.

    We moreover assume that `P` and `Q` have particular shapes
    (this is our use case):

    - `P` does not depend on `y`, i.e. it is a univariate
      polynomial in `x`

    - the `x`-degree of `Q` is strictly less than the `x`-degree
      of `P`.

    The resultant is computed with respect to the variable `x`.

    EXAMPLES::

        sage: from sage.rings.padics.padic_general_extension import resultant_bivariate
        sage: K = Qp(3, prec=5, print_mode='digits')
        sage: A.<y> = K[]
        sage: S.<x> = A[]
        sage: P = x^3 + 2*x^2 + 5*x + 3
        sage: Q = (x + 1)*y^2 + (x^2 + 4)*y + (2*x^2 + x)
        sage: resultant_bivariate(P, Q)
        ...00001*y^6 + ...00120*y^5 + ...21111*y^4 + ...21020*y^3 + ...01002*y^2 + ...01011*y + ...200210

    An example over an extension::

        sage: R.<w> = K[]
        sage: L.<w> = K.extension(w^2 + 3*w - 6)
        sage: A.<y> = K[]
        sage: S.<x> = A[]
        sage: P = x^2 + w*x + 1
        sage: Q = (x + 1)*y^2 + x*y + 1
        sage: resultant_bivariate(P, Q)
        ...1120021022*y^4 + ...1120021022*y^3 + ...2112002220*y^2 + ...1120021020*y + ...0000000001
    """
    Ky = P.base_ring()
    K = Ky.base_ring()
    name = Ky.variable_name()  # 'y' in our notations
    pi = K.uniformizer()
    e = K.absolute_e()

    Ku = K.maximal_unramified_subextension()
    S = PolynomialRing(Ku, names=name)

    # Instead of working over K[y], we work in the totally
    # ramified extension L = K[y]/(y^d - pi) for a suitable d
    d = P.degree() * max(c.degree() for c in Q.list()) + 1
    if e == 1:
        modulus = S.gen()**d - pi
        L = Ku.extension(modulus, names=name)
        base_map = None
    else:
        modulus = K.defining_polynomial()(S.gen()**d)
        L = Ku.extension(modulus, names=name)
        base_map = K.hom([L.gen()**d])
    f = Ky.hom([L.gen()], base_map=base_map)
    P = P.map_coefficients(f)
    Q = Q.map_coefficients(f)

    # We compute the resultant
    R = subresultant(P, Q)

    # We convert back the resultant to the correct ring
    # and return it
    if R.degree() > 0:
        return Ky.zero()
    coeffs = R[0].polynomial().list()
    return sum(S(coeffs[i*d:(i+1)*d]) * pi**i for i in range(e))


def newton_lift(P, x, prec=None):
    r"""
    Apply the Newton scheme to lift the approximate root `x` of `P`
    to an actual root.

    A ``PrecisionError`` is raised if the scheme does not seem to
    converge.

    INPUT:

    - ``P`` -- a `p`-adic polynomial

    - ``x`` -- an element in the base ring

    - ``prec`` -- an integer or ``None`` (default: ``None``); if
      given, stop the lifting once ``P(x)`` has attained valuation
      at least ``prec``

    EXAMPLES::

        sage: from sage.rings.padics.padic_general_extension import newton_lift
        sage: K = Qp(5, prec=10, print_mode='digits')
        sage: S.<x> = K[]
        sage: P = (x - 1/2) * (x - 1/3)
        sage: a = newton_lift(P, K(3)); a
        ...2222222223
        sage: a == 1/2
        True

        sage: newton_lift(P, K(3), prec=3)
        ...4013032223

    An example where the Newton scheme does not converge because the two roots
    are too close to each other::

        sage: Q = (x - 1/2) * (x + 9/2)
        sage: a = newton_lift(Q, K(3)); a
        Traceback (most recent call last):
        ...
        PrecisionError

    """
    v = 0
    dP = P.derivative()
    while True:
        x = x.lift_to_precision()
        num = P(x)
        if num == 0:
            return x
        vn = num.valuation()
        if prec is not None and vn >= prec:
            return x
        if vn > v:
            v = vn
        else:
            raise PrecisionError
        denom = dP(x)
        x -= num/denom


def factor_eisenstein(P, wK, e):
    r"""
    Return a factor of ``P``, together with its multiplicity.

    A ``PrecisionError`` is raised if the scheme does not seem to
    converge.

    INPUT:

    - ``P`` -- a polynomial over a `p`-adic field ``K``
      which is a product of conjugated Eisenstein polynomials

    - ``wK`` -- a uniformizer of ``K``

    - ``e`` -- the degree of the unramified extension
      defined by a factor of ``P``

    EXAMPLES::

        sage: from sage.rings.padics.padic_general_extension import factor_eisenstein
        sage: F = Qp(5)
        sage: K.<a> = F.extension(x^2 + 2)
        sage: S.<x> = K[]
        sage: P = x^4 - 25*a^2
        sage: factor_eisenstein(P, 5, 2)
        ((1 + O(5^20))*x^2 + a*5 + O(5^21), 1)

    An example with multiplicity:

        sage: factor_eisenstein(P^2, 5, 2)
        ((1 + O(5^20))*x^2 + a*5 + O(5^21), 2)

    Two examples where the input is corrupted::

        sage: R = x^4 - 5*a^2
        sage: factor_eisenstein(R, 5, 2)
        Traceback (most recent call last):
        ...
        PrecisionError

        sage: R = x^4 - 25*a
        sage: factor_eisenstein(R, 5, 2)
        Traceback (most recent call last):
        ...
        PrecisionError

    """
    K = P.base_ring()
    d = P.degree()

    # Basic check: if the slopes are not correct,
    # we raise an error
    if P.newton_slopes(repetition=False) != [ 1/e ]:
        raise PrecisionError
    n = d // e  # the number of irreducible factors

    # We compute roots in the residue field
    S0 = PolynomialRing(K.residue_field(), names='xe')
    P0 = S0([ (P[i*e] >> (n-i)).residue() for i in range(n+1) ])
    roots = P0.roots()
    nbroots = len(roots)
    if nbroots == 0:
        raise PrecisionError

    # We rule out multiplicities
    multiplicity = roots[0][1]
    if multiplicity * nbroots < n or any(m != multiplicity for _, m in roots):
        raise PrecisionError
    if multiplicity > 1:
        P = P // subresultant(P, P.derivative(), d - e*nbroots)

    # and lift it using a Newton iteration
    F = P.parent().gen()**e - K(roots[0][0]).lift_to_precision() * wK
    v = 0
    while True:
        Q, R = P.quo_rem(F)
        if R == 0:
            break
        vn = min(R[i].valuation() + i/e for i in range(e))
        if vn > v:
            v = vn
        else:
            raise PrecisionError
        _, _, C = F.xgcd(Q)  # can probably be improved
        F += (C*P) % F

    return F, multiplicity


def krasner_reduce(E):
    r"""
    Return an Eisenstein polynomial ``Ered`` which is small
    and defines the same extension as ``E``.

    INPUT:

    - ``E`` -- a monic Eisenstein polynomial over an unramified
      extension

    EXAMPLES::

        sage: from sage.rings.padics.padic_general_extension import krasner_reduce
        sage: R = Zp(3, print_mode='terse')
        sage: S.<x> = R[]
        sage: E = x^5 + 3 + 9 * S.random_element(degree=5)
        sage: krasner_reduce(E)
        (1 + O(3^20))*x^5 + 3 + O(3^21)

    An example conducting to a widly ramified extension::

        sage: E = cyclotomic_polynomial(81)(1+x)
        sage: E
        (1 + O(3^20))*x^54 + ... (2698955337 + O(3^20))*x^27 + ... + 3 + O(3^20)
        sage: krasner_reduce(E)
        (1 + O(3^20))*x^54 + ... (21 + O(3^21))*x^27 + ... + 3 + O(3^21)

    """
    d = E.degree()
    K = E.base_ring()
    p = K.prime()
    vco = [ c.valuation() for c in E.list() ]
    # valuation of the factorials
    vfa = [ 0 ]
    for i in range(1, d+1):
        vfa.append(vfa[i-1] + ZZ(i).valuation(p))
    for i in range(1,d):
        v = min(vfa[j] - vfa[i] - vfa[j-i] + vco[j] + (j-i)/d
                for j in range(i, d+1))
        if i == 1:
            valder = v
            slope = v / (d - 1)
        else:
            s = (valder - v) / (i - 1)
            if s > slope:
                slope = s
    val = valder + slope
    coeffs = [ ]
    from sage.functions.other import floor
    for i in range(d):
        prec = floor(val - i/d) + 1
        coeffs.append(E[i].add_bigoh(prec).lift_to_precision())
    coeffs.append(K.one())
    return E.parent()(coeffs)


# Classes
#########

class pAdicGeneralExtension(pAdicExtensionGeneric):
    r"""
    Shared base class for a general extension of a p-adic ring such as a
    relative extension or an extension not given by an unramified polynomial or
    an Eisenstein polynomial.

    EXAMPLES:

        sage: L.<a> = Qp(2).extension(x)

    """
    def __init__(self, exact_modulus, approx_modulus, prec, print_mode, shift_seed, names, element_class, category=None):
        r"""

        TESTS::

            sage: L.<a> = Qp(2).extension(x)
            sage: from sage.rings.padics.padic_general_extension import pAdicGeneralExtension
            sage: isinstance(L, pAdicGeneralExtension)
            True

        """
        base = approx_modulus.base_ring()
        self._exact_modulus = exact_modulus
        self._shift_seed = shift_seed
        self._implementation = 'proxy'
        self._prec_type = base._prec_type
        self.prime_pow = PowComputer_general(base.prime(), cache_limit=0, prec_cap=prec, ram_prec_cap=prec, in_field=base.is_field(), poly=approx_modulus)
        category = category or base.category()

        pAdicExtensionGeneric.__init__(self, exact_modulus, approx_modulus, prec, print_mode, names, element_class, category=category)

        if prec != self.base_ring().precision_cap():
            raise NotImplementedError("cannot change precision in general extension yet")

        if not self._exact_modulus.is_monic():
            raise NotImplementedError(f"defining modulus must be monic but {exact_modulus} is not")

    def relative_e(self):
        r"""
        Return the ramification degree of this ring over its base ring.

        EXAMPLES:

        A trivial extension::

            sage: L.<a> = Qp(2).extension(x + 2)
            sage: L.relative_e()
            1

        An unramified extension::

            sage: L.<a> = Qp(2).extension(x^2 + 2*x + 4)
            sage: L.relative_e()
            1

        An Eisenstein extension of an Eisenstein extension::

            sage: L.<a> = Zp(2).extension(x^2 - 2)
            sage: R.<b> = L[]
            sage: M.<b> = L.extension(b^3 - a)
            sage: M.relative_e()
            3

        """
        return self._e

    def relative_f(self):
        r"""
        Return the degree of the residual extension over its base ring.

        EXAMPLES:

        An unramified extension::

            sage: K.<a> = Qq(3^5)
            sage: K.relative_f()
            5

        ::

            sage: L.<a> = Qp(2).extension(x^2 + 2*x + 4)
            sage: L.relative_f()
            2

        A totally ramified extension::

            sage: L.<pi> = Qp(3).extension(x^2 - 3)
            sage: L.relative_f()
            1

        A trivial extension::

            sage: L.<a> = Qp(2).extension(x - 2)
            sage: L.relative_f()
            1

        """
        return self._f

    def absolute_ring(self, map=False):
        r"""
        Return an absolute extension of the absolute base isomorphic to this
        field.

        Note that this might not be a simple extension. It might be a p-adic
        base ring for a trivial extension or a two step extension, i.e., a
        totally ramified extension given by an Eisenstein polynomial over an
        unramified extension.

        EXAMPLES:

        A trivial extension::

            sage: L.<a> = Qp(2).extension(x + 2)
            sage: L.absolute_ring()
            2-adic Field with capped relative precision 20

        An unramified extension::

            sage: L.<a> = Qp(2).extension(x^2 + 2*x + 4)
            sage: L.absolute_ring()
            2-adic Unramified Extension Field in a_u defined by x^2 + x + 1

        Optionally, maps from and to the absolute extension are provided::

            sage: M, M_to_L, L_to_M = L.absolute_ring(map=True)
            sage: M_to_L(L_to_M(L.gen())) == L.gen()
            True

        """
        return backend_parent(self, map=map)

    def teichmuller(self, x, prec=None):
        R = self._backend
        x = R(x) if prec is None else R(x, prec)
        return self(R.teichmuller(x))

    def _prec_type(self):
        return self._backend._prec_type()

    def random_element(self, **kwds):
        return self(self._backend.random_element(**kwds))

    def residue_ring(self, n):
        raise NotImplementedError

    @cached_method
    def residue_class_field(self):
        r"""
        Return the residue class field of this ring.

        EXAMPLES:

        A trivial extension::

            sage: L.<a> = Qp(2).extension(x + 2)
            sage: L.residue_class_field()
            Trivial extension of Finite Field of size 2

        A trivial extension of a trivial extension::

            sage: R.<b> = L[]
            sage: M.<b> = L.extension(b - a)
            sage: M.residue_field()
            Trivial extension of Trivial extension of Finite Field of size 2

        An unramified extension::

            sage: L.<a> = Qp(2).extension(x^2 + 2*x + 4)
            sage: L.residue_class_field()
            Finite Field in z2 of size 2^2 over its base

        A trivial extension of an unramified extension::

            sage: R.<b> = L[]
            sage: M.<b> = L.extension(b - 2)
            sage: M.residue_field()
            Trivial extension of Finite Field in z2 of size 2^2 over its base

        An unramified extension of an unramified extension::

            sage: L.<a> = Qp(2).extension(x^2 + 2*x + 4)
            sage: R.<b> = L[]
            sage: M.<b> = L.extension(b^2 + a*b + 4)
            sage: m = M.residue_field()
            sage: m
            Finite Field in z4 of size 2^4 over its base
            sage: m.base_ring() is L.residue_field()
            True
            sage: m.modulus()
            x^2 + x + z2

        """
        return self.base_ring().residue_class_field().extension(self.relative_f(), absolute=False, implementation="GF", backend=self._backend.residue_class_field())

    def uniformizer(self):
        backend, from_backend, _ = backend_parent(self, map=True)
        return from_backend(backend.uniformizer())

    def uniformizer_pow(self, n):
        backend, from_backend, _ = backend_parent(self, map=True)
        return from_backend(backend.uniformizer_pow(n))

    def _uniformizer_print(self):
        return self._backend._uniformizer_print()

    def gen_unram(self):
        backend, from_backend, _ = backend_parent(self, map=True)
        return from_backend(backend.gen_unram())

    def _unram_print(self):
        return self._backend._unram_print()

    def has_pth_root(self):
        return self._backend.has_pth_root()

    def has_root_of_unity(self, n):
        return self._backend.has_root_of_unity(self, n)

    def construction(self, forbid_frac_field=None):
        # Prefer AlgebraicExtensionFunctor for pushout since FractionField
        # functor often does not work because there is no integer_ring.
        if forbid_frac_field is None:
            forbid_frac_field = True

        # TODO: Change prec of AlgebraicExtensionFunctor
        construction = pAdicExtensionGeneric.construction(self, forbid_frac_field=forbid_frac_field)
        return construction

    def inertia_subring(self):
        r"""
        Return the inertia subring of this extension over its base.

        EXAMPLES:

        The inertia subring of an unramified extension is the ring itself::

            sage: L.<a> = Zp(2).extension(x^2 + 2*x + 4)
            sage: L.inertia_subring() is L
            True

        A trivial extension is added to the inertia subring::

            sage: R.<b> = L[]
            sage: M.<b> = L.extension(b - a)
            sage: M.inertia_subring() is M
            True

        A ramified extension of an unramified extension::

            sage: R.<b> = L[]
            sage: M.<b> = L.extension(b^2 - 2)
            sage: M.inertia_subring() is L
            True

        A general extension. The inertia subring is an extension of the ground
        field::

            sage: L.<a> = Zp(2).extension(x^4 + 8*x^2 + 64)
            sage: M = L.inertia_subring(); M
            2-adic Unramified Extension Ring in a_u defined by x^2 + x + 1
            sage: M.base_ring() is Zp(2)
            True

        """
        if self.relative_e() == 1:
            return self
        if self.relative_f() == 1:
            return self.base_ring()

        backend, from_backend, to_backend = backend_parent(self, map=True)

        backend_inertia_subring = self._backend.inertia_subring()
        inertia_subring_generator = from_backend(backend_inertia_subring.gen())
        inertia_subring = self.base_ring().extension(inertia_subring_generator.minpoly(), names=backend_inertia_subring._unram_print())
        return inertia_subring


class pAdicGeneralRingExtension(pAdicGeneralExtension, RingExtension_generic):
    r"""
    A general extension of a p-adic ring such as a relative extension or an
    extension not given by an unramified polynomial or an Eisenstein
    polynomial.

    This class models rings that are not fields.

    .. SEEALSO:: :class:`pAdicGeneralFieldExtension`

    """
    def __init__(self, exact_modulus, approx_modulus, prec, print_mode, shift_seed, names, implementation='FLINT', category=None):
        base = approx_modulus.base_ring()
        category = category or base.category()

        pAdicGeneralExtension.__init__(self, exact_modulus=exact_modulus, approx_modulus=approx_modulus, prec=prec, print_mode=print_mode, shift_seed=shift_seed, names=names, element_class=pAdicGeneralRingExtensionElement, category=category)

        _, _, to_fraction_field_backend = backend_parent(self.fraction_field(), map=True)

        defining_morphism = pAdicGeneralRingMap_FractionField(self.base_ring(), self.fraction_field()._backend.integer_ring(), to_fraction_field_backend)
        self._backend = defining_morphism.codomain()
        self._prec = self.fraction_field()._prec
        self._e = self.fraction_field()._e
        self._f = self.fraction_field()._f

        RingExtension_generic.__init__(self, defining_morphism=defining_morphism, import_methods=False, category=category, check=False)

        pAdicGeneralMap_Backend(self.fraction_field(), self).register_as_conversion()

    def is_field(self, *args, **kwds):
        r"""
        Return whether this ring is a field.

        EXAMPLES::

            sage: L.<a> = Zp(2).extension(x + 3)
            sage: L.is_field()
            False

        """
        return False

    def gen(self, i=0):
        r"""
        Return the generator of this ring's fraction field over its base.

        EXAMPLES::
        
            sage: L.<a> = Zp(2).extension(x^2 + 2*x + 4)
            sage: L.gen()
            (a_u + 1)*2 + (a_u + 1)*2^2 + (a_u + 1)*2^3 + (a_u + 1)*2^4 + (a_u + 1)*2^5 + (a_u + 1)*2^6 + (a_u + 1)*2^7 + (a_u + 1)*2^8 + (a_u + 1)*2^9 + (a_u + 1)*2^10 + (a_u + 1)*2^11 + (a_u + 1)*2^12 + (a_u + 1)*2^13 + (a_u + 1)*2^14 + (a_u + 1)*2^15 + (a_u + 1)*2^16 + (a_u + 1)*2^17 + (a_u + 1)*2^18 + (a_u + 1)*2^19 + (a_u + 1)*2^20 + O(2^21)
            sage: L.gen().parent() is L
            True

        """
        return self(self.fraction_field().gen(i))

    def integer_ring(self):
        # TODO: Is this the default implementation anyway?
        return self

    def exact_ring(self):
        # TODO: The naive approach of the base class is typically not correct
        # here anymore.
        raise NotImplementedError

    @cached_method
    def fraction_field(self):
        # TODO: (see Merge Request 32 - integer_ring)
        return self.construction()[0](self.base_ring().fraction_field())


class pAdicGeneralFieldExtension(pAdicGeneralExtension, RingExtensionWithGen):
    def __init__(self, exact_modulus, approx_modulus, prec, print_mode, shift_seed, names, implementation='FLINT', category=None):
        base = approx_modulus.base_ring()
        category = category or base.category()

        pAdicGeneralExtension.__init__(self, exact_modulus=exact_modulus, approx_modulus=approx_modulus, prec=prec, print_mode=print_mode, shift_seed=shift_seed, names=names, element_class=pAdicGeneralFieldExtensionElement, category=category)

        # TODO: (see Merge Request 32 - implementation)
        defining_morphism, gen = self._create_backend()

        self._backend = gen.parent()

        # Patch prec which was set not knowing the ramification index.
        self._prec = prec * self.relative_e()

        RingExtensionWithGen.__init__(self, defining_morphism=defining_morphism, gen=gen, names=[self.variable_name()], category=category, import_methods=False, check=False)

    degree = RingExtensionWithGen.degree
    # Use the implementation of __reduce__ from the factory and ignore RingExtensionWithGen's override.
    __reduce__ = pAdicGeneralExtension.__reduce__

    def is_field(self, *args, **kwds):
        r"""
        Return whether this ring is a field.

        EXAMPLES::

            sage: L.<a> = Qp(2).extension(x)
            sage: L.is_field()
            True

        """
        return True

    @cached_method
    def integer_ring(self):
        # TODO: (see Merge Request 32 - integer_ring)
        return self.construction()[0](self.base_ring().integer_ring())

    def fraction_field(self):
        # TODO: Is this the default implementation anyway?
        return self

    def _coerce_map_from_(self, R):
        if isinstance(R, pAdicGeneralRingExtension) and R is self.integer_ring():
            return pAdicGeneralMap_Backend(R, self)
        return RingExtensionWithGen._coerce_map_from_(self, R)

    def _create_backend_exact(self):
        r"""
        Return a backend for this extension, i.e., a p-adic ring that is not a
        general extension itself.
        """
        if not self._exact_modulus.is_squarefree():
            # We only check squarefreeness here. Irreducibility is checked
            # automatically, when the extensions of the valuations on base to
            # the ring are constructed. (If there is more than one extension,
            # i.e., the polynomial is not irreducible, exact_valuation() is
            # going to complain.)
            raise ValueError("polynomial must be irreducible but %r is not"%(polynomial,))

        val = self.exact_valuation()
        self._e = val.E()
        self._f = val.F()

        if self._f == 1 and self._e == 1:
            # This is a trivial extension. The best backend is base ring
            # (possibly rewritten as an absolute extension.)
            assert self._exact_modulus.degree() == 1

            (backend, backend_to_base, base_to_backend) = self.base_ring().absolute_ring(map=True)
            defining_morphism = base_to_backend
            gen = defining_morphism(self.base_ring()(-self._exact_modulus[0]))
        else:
            # The underlying Zp or Qp
            backend_base = self.ground_ring_of_tower()

            # The unramified part of this extension.
            if self.absolute_f() == 1:
                backend_unramified = backend_base
            else:
                backend_unramified = self.ground_ring_of_tower().change(q=self.prime()**self.absolute_f(), names=self._printer.unram_name)

            # The totally ramified part of this extension.
            if self.absolute_e() == 1:
                backend = backend_unramified
            else:
                # We construct the charpoly of the uniformizer and factor it
                # over the unramified part. Currently, we do this completely
                # naively in the corresponding number field which is terribly
                # slow.
                charpoly = self.exact_field().absolute_field('x').valuation(self.prime()).uniformizer().charpoly()

                assert charpoly.degree() == self.absolute_e() * self.absolute_f()

                charpoly = charpoly.change_ring(backend_unramified.exact_field())
                # The charpoly is a power p-adically. Typically, it's
                # irreducible in the number field but it might not actually be.
                # The Montes factorization assumes that it's squarefree so we
                # make sure that's the case and focus of any of the equivalent
                # factors.
                charpoly = charpoly.squarefree_decomposition()[0][0]
                # We factor the charpoly over the number field to the required
                # precision (note that we do not divide the precision with the
                # absolute e of this ring since the rescaling of _prec in
                # __init__ has not been performed yet.)
                # We could do much better here since we do not need all factors
                # but just one of them.
                # Also, we do not need an actual fator of charpoly but just an
                # approximation that singles out the correct totally ramified
                # extension, i.e., we could apply some Krasner bound argument.
                charpoly = backend_unramified.exact_field().valuation(self.prime()).montes_factorization(charpoly, assume_squarefree=True, required_precision=self._prec // self.base_ring().absolute_e())

                assert all(f.degree() == self.absolute_e() for f,e in charpoly), f"charpoly of uniformizer did not factor as an approximate {self.absolute_f()}th power: {charpoly}"

                minpoly = charpoly[0][0]

                backend = backend_unramified.extension(minpoly, names='pi')

            def create_defining_morphism(base=self.base_ring()):
                if base is base.ground_ring_of_tower():
                    return base.hom(backend)

                base_map = create_defining_morphism(base.base_ring())

                # TODO: The poly.change_ring() might not have enough precision.
                # TODO: The any_root() might not have enough precision.
                modulus = base.modulus().change_ring(base_map)
                return base.hom([modulus.any_root()], codomain=backend, base_map=base_map)

            defining_morphism = create_defining_morphism()

            # TODO: The poly.change_ring() might not have enough precision.
            # TODO: The any_root() might not have enough precision.
            gen = self.defining_polynomial().change_ring(defining_morphism).any_root()

        if backend is not backend.absolute_ring():
            raise NotImplementedError("relative backends are not supported for general p-adic extensions yet")

        return defining_morphism, gen

    def _create_backend_padic(self, reduce=True, verbose=False):  # eventually, we will remove verbose
        r"""
        Return a backend for this extension, i.e.,
        a `p`-adic ring that is not a general extension itself.

        ALGORITHM:

        Let `K` be the backend of the base and
        let `K^u` be its maximum unramified subextension.

        Let `P` be the defining polynomial of this extension
        and let `a` be the corresponding generator (that is,
        a root of `P`).

        1. We compute a uniformizer `\pi` of this extension,
           together with its relative ramification index `e`
           and its relative residual degree `f`

        2. We create the unramified extension of `K^u` of degree
           `f` and compute the embedding `f^u : K^u -> L^u`

        3. We create the compositum of `K` and `L^u`: it is just
           the extension of `L^u` defined by the same Eisenstein
           polynomial which defined the extension `K/K^u`.

        4. If `e = 1`, we have finished: `K L^u` is our backend!

        5. Otherwise, we compute the characteristic polynomial
           `\chi` of `\pi` over `K`. It is given as a resultant:
           if `\Pi` is the polynomial such that `\pi = \Pi(a)``,
           then `\chi(Y) = \text{Res}_X(P(X), Y - \Pi(X))`

        6. We compute the minimal polynomial `\mu` of `\pi` over
           `K L^u`. This is done by factoring `\chi` over `K L^u`.

        7. We compute the minimal polynomial `E` of `\pi` over
           `L^u`. Again, it is given a resultant.
           Indeed, let `E_K` be  the Eisenstein polynomial defining
           `K/K^u` and let `\pi_K` be the corresponding uniformizer.
           Let also `\Chi(X,Y)` be a bivariate polynomial over `K^u`
           such that `\chi(Y) = \Chi(\pi_K, Y)`.
           Then `E(Y) = \text{Res}_X(E_K(X), \Chi(X,Y))`.

        8. We replace `E` by a polynomial `E_{\text{red}}` defining
           the same extension, but exhibiting smaller coefficients.
           This can be done using Krasner's lemma.

        9. We create the Eisentein extension `L` of `L^u`
           defined by the polynomial `E_{\text{red}}`:
           it is our backend!

        EXAMPLES::

            sage: Ku.<a> = Qq(5^3)
            sage: S.<x> = Ku[]
            sage: K.<w> = Ku.extension(x^4 - 5*a)
            sage: S.<x> = K[]
            sage: P = x^4 + w*x^3 + a*w*x^2 + (2*a^2 + 3)*w^2*x + w^2
            sage: L.<b> = K.extension(P)  # indirect doctest

            sage: L._backend
            5-adic Eisenstein Extension Field in b_p defined by y^8 + 5*b_u^5 + 5*b_u^4 + 15*b_u^3 + 15*b_u + 10 over its base field

            sage: L._backend.base_ring()
            5-adic Unramified Extension Field in b_u defined by x^6 + x^4 + 4*x^3 + x^2 + 2

        """
        from sage.misc.misc import walltime
        t0 = walltime()

        if not self._exact_modulus.is_squarefree():
            # We only check squarefreeness here. Irreducibility is checked
            # automatically, when the extensions of the valuations on base to
            # the ring are constructed. (If there is more than one extension,
            # i.e., the polynomial is not irreducible, exact_valuation() is
            # going to complain.)
            raise ValueError("polynomial must be irreducible but %r is not" % self._given_poly)

        p = self.prime()
        name = self._names[0]
        K, from_K, to_K = backend_parent(self._base, map=True)
        F = K.ground_ring_of_tower()
        Ku = K.maximal_unramified_subextension()
        if isinstance(self._base, RingExtension_generic):
            P = self._given_poly.map_coefficients(to_K)
        else:
            P = self._given_poly.change_ring(K)
        Pex = P.map_coefficients(lambda x: x.lift_to_precision())
        is_base_unramified = (K.absolute_e() == 1)

        # Step 1: uniformizer, e, f
        t = walltime()
        vals = K.valuation().mac_lane_approximants(Pex, assume_squarefree=True)
        if len(vals) > 1:
            raise ValueError("polynomial must be irreducible but %r is not" % self._given_poly)
        val = vals[0]
        pi = val.uniformizer()
        e_rel = self._e = val.E()
        f_rel = self._f = val.F()
        if verbose:
            print("# uniformizer, e and f computed in %.3fs [e = %s, f = %s]" % (walltime(t), e_rel, f_rel))

        # Step 2: Lu and embedding fu : Ku -> Lu
        t = walltime()
        k = Ku.residue_field()
        if f_rel == 1:
            Lu = Ku
            l = k
            fu = None
        else:
            # Lu
            l = k.extension(ZZ(f_rel), absolute=True)
            modulus = l.modulus().change_ring(ZZ).change_ring(F)
            Lu = F.extension(modulus, names = name + '_u', absolute=True)
            # Fu : Ku -> Lu
            if Ku.absolute_f() == 1:
                fu = Lu.coerce_map_from(Ku)
            else:
                U = Ku.modulus().change_ring(Lu)
                U0 = U.change_ring(Lu.residue_field())
                a0 = U0.any_root()
                a = U.hensel_lift(Lu(a0).lift_to_precision())
                fu = Ku.hom([a])
        if verbose:
            print("# embedding Ku -> Lu computed in %.3fs" % walltime(t))

        # Step 3: KLu and embedding g : K -> KLu
        t = walltime()
        if is_base_unramified:
            KLu = Lu
            wK = KLu(p)
            g = fu
        else:
            EK = K.defining_polynomial()
            if fu is not None:
                EK = EK.map_coefficients(fu)
            KLu = Lu.extension(EK, names='wK', absolute=True)
            wK = KLu.uniformizer()
            g = K.hom([wK], base_map=fu)
        if verbose:
            print("# embedding K -> KLu computed in %.3fs" % walltime(t))

        if e_rel == 1:

            # Step 4: unramified extension
            f = g
            if f is None:  # case of trivial extension
                from sage.categories.homset import End
                f = End(K).identity()

        else:

            # Step 5: characteristic polynomial of pi over K
            t = walltime()
            S = PolynomialRing(K, names='y')
            T = PolynomialRing(S, names = P.variable_name())
            charpoly = resultant_bivariate(P.change_ring(S), S.gen() - pi.change_ring(S))
            if verbose:
                print("# characteristic polynomial of pi over K computed in %.3fs" % walltime(t))

            # Step 6: minimal polynomial of pi over KLu
            t = walltime()
            if g is not None:
                charpoly = charpoly.map_coefficients(g)
            minpoly, multiplicity = factor_eisenstein(charpoly, wK, e_rel)
            if verbose:
                print("# minimal polynomial of pi over KLu computed in %.3fs" % walltime(t))

            # Step 7: minimal polynomial of pi over Lu
            t = walltime()
            S = PolynomialRing(Lu, name='y')
            if is_base_unramified:
                # In this case Lu = KLu and there is nothing to do
                E = minpoly
            else:
                Kname = EK.variable_name()
                T = PolynomialRing(S, names=Kname)
                y = S.gen()
                Q = sum(minpoly[i].polynomial(Kname).change_ring(S) * y**i
                        for i in range(e_rel + 1))
                E = resultant_bivariate(EK.change_ring(S), Q).monic()
            if verbose:
                print("# minimal polynomial of pi over Lu computed in %.3fs" % walltime(t))

            if reduce:
                # Step 8: reduce the Eisenstein polynomial
                t = walltime()
                Ered = krasner_reduce(E)
                if verbose:
                    print("# minimal polynomial of pi over Lu reduced in %.3fs" % walltime(t))
            else:
                E = Ered

            # Step 9: L and embedding f : K -> L
            t = walltime()
            L = Lu.extension(Ered, names = name + '_p')
            wL = L.uniformizer()
            if is_base_unramified:
                f = L.coerce_map_from(Lu)
                if fu is not None:
                    f = f * fu
            else:
                # Since we have reduced E, the canonical uniformizer of L
                # (namely L.uniformizer()) is not a root of E.
                # We first compute a root wL of E using a Newton scheme
                # starting from L.uniformizer().
                # NOTE: It seems actually that this costly step is not
                # needed; so we just skip it.
                #if reduce:
                #    wL = newton_lift(E.change_ring(L), wL)
                # We then compute the image wK of K.uniformizer() in L by
                # solving the equation Q(x, wL) = 0 (of which wK is a solution)
                # Actually, because of inaccuracy issues (and also possibly
                # because we do not use the correct wL), this will not be enough
                # to give the correct answer and we shall need to apply a last
                # Newton scheme afterwards. For this reason, it's more clever
                # to lift at limited precision right now, just to ensure that
                # the final Newton scheme will converge.
                QwL = Q.map_coefficients(lambda C: C(wL), new_base_ring=L)  # it's Q(x, wL)
                eK = K.absolute_e()
                valder = min(i + eK * (ZZ(i+1).valuation(p) + EK[i+1].valuation())
                             for i in range(EK.degree()))
                prec = 2*e_rel*valder + 1
                wK = newton_lift(QwL, L.zero(), prec)
                # Final Newton lift
                wK = newton_lift(EK.change_ring(L), wK)
                # We are finally ready to define f
                f = K.hom([wK], base_map=fu, check=False)  # already checked in newton_lift
            if verbose:
                print("# embedding K -> L computed in %.3fs" % walltime(t))

        # Computation of the generator and final validation
        t = walltime()
        A = P.map_coefficients(f)
        if e_rel > 1:
            # In this case, we have the generator is solution of
            # two equations, namely:
            # - P(gen) = 0
            # - Pi(gen) = wL
            # So, in order to lower the degree, we take the gcd
            # of these two polynomials.
            # It turns out that the gcd method is not numerically
            # stable. However, in our setting, we know in advance
            # the degree of the gcd (it's `multiplicity`), so we
            # can instead compute a subresultant.
            B = pi.map_coefficients(f) - wL
            R = subresultant(A, B, multiplicity)
            if R.degree() > multiplicity:
                raise PrecisionError
            try:
                gen = R.any_root()
            except ValueError:
                raise PrecisionError
            # We correct the generator in case of numerical inaccuracy
            gen = newton_lift(A, gen)
        else:
            # Here the previous trick does not apply
            # So we just compute a root of the given polynomial
            try:
                gen = A.any_root()
            except ValueError:
                raise PrecisionError
        if verbose:
            print("# generator computed in %.3fs" % walltime(t))
            print("# total walltime: %.3fs" % walltime(t0))
        self._backend_defining_morphism = f
        return f, gen

    def _create_backend(self):
        r"""
        Return a backend for this extension, i.e.,
        a `p`-adic ring that is not a general extension itself.
        """
        # TODO: (see Merge Request 32 - backend)
        try:
            return self._create_backend_padic()
        except RuntimeError:
            from warnings import warn
            warn("computation failed p-adically because precision was not enough, we're running an exact computation now (might be slow)")
            return self._create_backend_exact()


class PowComputer_general(PowComputer_class):
    pass


class pAdicGeneralMap_Backend(RingMap):
    def _init__(self, domain, codomain):
        # TODO: Use a better category.
        from sage.categories.sets_with_partial_maps import SetsWithPartialMaps

        RingMap.__init__(self, Hom(domain, codomain, SetsWithPartialMaps()))

    def _call_(self, x):
        _, _, to_domain_backend = backend_parent(self.domain(), map=True)
        _, from_codomain_backend, _ = backend_parent(self.codomain(), map=True)
        return from_codomain_backend(to_domain_backend(x))

    def _call_with_args(self, x, args=(), kwds={}):
        return self.codomain()._element_constructor_(self._call_(x), *args, **kwds)


class pAdicGeneralRingMap_FractionField(RingMap):
    def __init__(self, domain, codomain, fraction_field_map):
        self._fraction_field_map = fraction_field_map

        RingMap.__init__(self, Hom(domain, codomain))

    def _call_(self, x):
        x = self.domain().fraction_field()(x)
        y = self._fraction_field_map(x)
        return self.codomain()(y)

    def _call_with_args(self, x, args=(), kwds={}):
        return self.codomain()._element_constructor_(self._call_(x), *args, **kwds)
