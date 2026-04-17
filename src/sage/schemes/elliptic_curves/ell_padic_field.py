# sage.doctest: needs sage.rings.padics
"""
Elliptic curves over `p`-adic fields

AUTHORS:

- Robert Bradshaw
- William Stein
- Miljan Brakovevic
- Ralf Gerkmann
- Kiran Kedlaya
- Jennifer Balakrishnan
- Francis Clarke
"""
# ***************************************************************************
#       Copyright (C) 2007-2010 Robert Bradshaw <robertwb@math.washington.edu>
#                     2007      William Stein   <wstein@gmail.com>
#                     2007      Miljan Brakovevic
#                     2007      Ralf Gerkmann
#                     2007-2008 Kiran Kedlaya <kedlaya@math.mit.edu>
#                     2007-2010 Jennifer Balakrishnan
#                     2012      Francis Clarke
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
# ***************************************************************************

from .ell_field import EllipticCurve_field
from .constructor import EllipticCurve
from . import ell_point
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.misc.lazy_import import lazy_import
from sage.rings.power_series_ring import PowerSeriesRing
from sage.modules.free_module_element import vector
from sage.matrix.constructor import matrix
from sage.modules.free_module import VectorSpace
from sage.rings.infinity import Infinity
from sage.rings.laurent_series_ring import LaurentSeriesRing
from sage.rings.rational_field import QQ
from sage.rings.big_oh import O
from sage.rings.rational_field import RationalField
from sage.rings.integer_ring import ZZ

lazy_import("sage.rings.padics.factory", "Qp", as_="pAdicField")
lazy_import("sage.schemes.hyperelliptic_curves", "monsky_washnitzer")


class EllipticCurve_padic_field(EllipticCurve_field):
    """
    Elliptic curve over a `p`-adic field.

    EXAMPLES::

        sage: Qp = pAdicField(17)
        sage: E = EllipticCurve(Qp,[2,3]); E
        Elliptic Curve defined by y^2  = x^3 + (2+O(17^20))*x + (3+O(17^20))
         over 17-adic Field with capped relative precision 20
        sage: E == loads(dumps(E))
        True
    """

    _point = ell_point.EllipticCurvePoint_field

    _genus = 1

    def frobenius(self, P=None):
        """
        Return the Frobenius as a function on the group of points of
        this elliptic curve.

        EXAMPLES::

            sage: Qp = pAdicField(13)
            sage: E = EllipticCurve(Qp,[1,1])
            sage: type(E.frobenius())
            <... 'function'>
            sage: point = E(0,1)
            sage: E.frobenius(point)
            (0 : 1 + O(13^20) : 1 + O(13^20))

        Check that :issue:`29709` is fixed::

            sage: Qp = pAdicField(13)
            sage: E = EllipticCurve(Qp,[0,0,1,0,1])
            sage: E.frobenius(E(1,1))
            Traceback (most recent call last):
            ...
            NotImplementedError: Curve must be in weierstrass normal form.
            sage: E = EllipticCurve(Qp,[0,1,0,0,1])
            sage: E.frobenius(E(0,1))
            (0 : 1 + O(13^20) : 1 + O(13^20))
        """
        if not hasattr(self, "_frob"):
            a1, a2, a3, a4, a6 = self.a_invariants()
            if a1 or a3:
                raise NotImplementedError("Curve must be in weierstrass normal form.")

            K = self.base_field()
            p = K.prime()
            x = PolynomialRing(K, "x").gen(0)

            f = x**3 + a2 * x**2 + a4 * x + a6
            h = f(x**p) - f**p

            def _frob(P):
                x0 = P[0]
                y0 = P[1]
                uN = (1 + h(x0) / y0 ** (2 * p)).sqrt()
                yres = y0**p * uN
                xres = x0**p
                if (yres - y0).valuation() == 0:
                    yres = -yres
                return self.point([xres, yres, K(1)])

            self._frob = _frob

        if P is None:
            return self._frob
        return self._frob(P)

    # =============================================================================
    # The functions below were prototyped at the 2007 Arizona Winter School by
    # Robert Bradshaw and Ralf Gerkmann, working with Miljan Brakovevic and
    # Kiran Kedlaya
    #
    # All of the below is with respect to the Monsky Washnitzer cohomology.
    #
    # NOTE: these functions were then taken from the generic HyperellipticCurve
    # and specialised to genus one in preparation of the hyperelliptic curve code
    # moving to the weighted projective model.
    # =============================================================================

    def local_coordinates_at_nonweierstrass(self, P, prec=20, name="t"):
        """
        For a non-Weierstrass point `P = (a,b)` on the elliptic
        curve `y^2 = f(x)`, return `(x(t), y(t))` such that `(y(t))^2 = f(x(t))`,
        where `t = x - a` is the local parameter.

        INPUT:

        - ``P = (a, b)`` -- a non-Weierstrass point on ``self``
        - ``prec`` -- desired precision of the local coordinates
        - ``name`` -- gen of the power series ring (default: ``t``)

        OUTPUT:

        `(x(t),y(t))` such that `y(t)^2 = f(x(t))` and `t = x - a`
        is the local parameter at `P`

        EXAMPLES::

            sage: Qp = pAdicField(13)
            sage: E = EllipticCurve(Qp,[5,22])
            sage: P = E([-1,4])
            sage: xt,yt = E.local_coordinates_at_nonweierstrass(P, prec=5)
            sage: xt
            12 + 12*13 + 12*13^2 + 12*13^3 + 12*13^4 + 12*13^5 + 12*13^6 + 12*13^7 + 12*13^8 + 12*13^9 + 12*13^10 + 12*13^11 + 12*13^12 + 12*13^13 + 12*13^14 + 12*13^15 + 12*13^16 + 12*13^17 + 12*13^18 + 12*13^19 + O(13^20) + (1 + O(13^20))*t + O(t^5)
            sage: yt
            4 + O(13^20) + (1 + O(13^20))*t + (6 + 6*13 + 6*13^2 + 6*13^3 + 6*13^4 + 6*13^5 + 6*13^6 + 6*13^7 + 6*13^8 + 6*13^9 + 6*13^10 + 6*13^11 + 6*13^12 + 6*13^13 + 6*13^14 + 6*13^15 + 6*13^16 + 6*13^17 + 6*13^18 + 6*13^19 + O(13^20))*t^2 + (10 + 9*13 + 9*13^2 + 9*13^3 + 9*13^4 + 9*13^5 + 9*13^6 + 9*13^7 + 9*13^8 + 9*13^9 + 9*13^10 + 9*13^11 + 9*13^12 + 9*13^13 + 9*13^14 + 9*13^15 + 9*13^16 + 9*13^17 + 9*13^18 + 9*13^19 + O(13^20))*t^3 + (6 + 4*13 + 9*13^2 + 7*13^3 + 12*13^4 + 10*13^5 + 2*13^6 + 13^7 + 6*13^8 + 4*13^9 + 9*13^10 + 7*13^11 + 12*13^12 + 10*13^13 + 2*13^14 + 13^15 + 6*13^16 + 4*13^17 + 9*13^18 + 7*13^19 + O(13^20))*t^4 + O(t^5)

        AUTHOR:

        - Jennifer Balakrishnan (2007-12)
        """
        if P[2].is_zero():
            raise ValueError(f"P = {P} is the point at infinity. Use local_coordinates_at_infinity instead!")

        d = P[1]
        if d.is_zero():
            raise ValueError(
                f"P = {P} is a Weierstrass point. Use local_coordinates_at_weierstrass instead!"
            )

        pol = self.hyperelliptic_polynomials()[0]
        L = PowerSeriesRing(self.base_ring(), name, default_prec=prec)
        t = L.gen()
        K = PowerSeriesRing(L, "x")
        pol = K(pol)
        b = P[0]
        f = pol(t + b)
        for _ in range(prec.bit_length()):
            d = (d + f / d) / 2
        return t + b + O(t**prec), d + O(t**prec)

    def local_coordinates_at_weierstrass(self, P, prec=20, name="t"):
        """
        For a finite Weierstrass point on the elliptic
        curve `y^2 = f(x)`, return `(x(t), y(t))` such that
        `(y(t))^2 = f(x(t))`, where `t = y` is the local parameter.

        INPUT:

        - ``P`` -- a finite Weierstrass point on ``self``
        - ``prec`` -- desired precision of the local coordinates
        - ``name`` -- gen of the power series ring (default: `t`)

        OUTPUT:

        `(x(t),y(t))` such that `y(t)^2 = f(x(t))` and `t = y`
        is the local parameter at `P`

        EXAMPLES::

            sage: Q5 = pAdicField(5)
            sage: E = EllipticCurve(Q5, [1,0])
            sage: P = E([0,0])
            sage: xt,yt = E.local_coordinates_at_weierstrass(P, prec=5)
            sage: xt
            (1 + O(5^20))*t^2 + O(t^5)
            sage: yt
            (1 + O(5^20))*t + O(t^5)
            sage: O = E([0,1,0])
            sage: xt,yt = E.local_coordinates_at_weierstrass(O, prec=5)
            Traceback (most recent call last):
            ...
            ValueError:  P = (0 : 1 + O(5^20) : 0) is not a finite Weierstrass point. Use local_coordinates_at_nonweierstrass instead!

        AUTHOR:
          - Jennifer Balakrishnan (2007-12)
          - Francis Clarke (2012-08-26)
        """
        #  Ensure the input point is Weierstrass
        if not P[1].is_zero():
            raise ValueError(
                    f"P = {P} is not a finite Weierstrass point. Use local_coordinates_at_nonweierstrass instead!"
                )

        if P[2].is_zero():
            raise ValueError(f"P = {P} is the point at infinity. Use local_coordinates_at_infinity instead!")

        L = PowerSeriesRing(self.base_ring(), name)
        t = L.gen()
        pol = self.hyperelliptic_polynomials()[0]
        pol_prime = pol.derivative()
        b = P[0]
        t2 = t**2
        c = b + t2 / pol_prime(b)
        c = c.add_bigoh(prec)
        for _ in range(prec.bit_length()):
            c -= (pol(c) - t2) / pol_prime(c)
        return (c, t.add_bigoh(prec))

    def local_coordinates_at_infinity(self, prec=20, name="t"):
        """
        For the elliptic curve `y^2 = f(x)`, return
        `(x(t), y(t))` such that `(y(t))^2 = f(x(t))`, where `t = x/y` is
        the local parameter at infinity

        INPUT:

        - ``prec`` -- desired precision of the local coordinates
        - ``name`` -- generator of the power series ring (default: ``t``)

        OUTPUT:

        `(x(t),y(t))` such that `y(t)^2 = f(x(t))` and `t = x/y`
        is the local parameter at infinity

        EXAMPLES::

            sage: Q5 = pAdicField(5)
            sage: E = EllipticCurve(Q5, [1,0])
            sage: xt,yt = E.local_coordinates_at_infinity(prec=5)
            sage: xt[-2]
            1 + O(5^20)
            sage: yt[-3]
            1 + O(5^20)

        AUTHOR:

        - Jennifer Balakrishnan (2007-12)
        """
        pol = self.hyperelliptic_polynomials()[0]
        K = LaurentSeriesRing(self.base_ring(), name, default_prec=prec + 2)
        t = K.gen()
        L = PolynomialRing(K, "x")
        x = L.gen()
        w = (x / t) ** 2 - pol
        wprime = w.derivative(x)
        x = t**-2

        for _ in range((prec + 1).bit_length()):
            x = x - w(x) / wprime(x)
        y = x / t
        return x + O(t ** (prec + 2)), y + O(t ** (prec + 2))

    def local_coord(self, P, prec=20, name="t"):
        """
        Return the local coordinates of the elliptic curve at `P`.

        TODO: extend to general Weierstrass form, and move to ell_generic ?

        INPUT:

        - ``P`` -- a point on ``self``
        - ``prec`` -- desired precision of the local coordinates
        - ``name`` -- generator of the power series ring (default: ``t``)

        OUTPUT:

        `(x(t),y(t))` such that `y(t)^2 = f(x(t))`, where `t`
        is the local parameter at `P`

        EXAMPLES::

            sage: Qp = pAdicField(11)
            sage: E = EllipticCurve(Qp,[0,1,0,0,4])
            sage: f,_ = E.hyperelliptic_polynomials()
            sage: P1 = E([-2,0])
            sage: x1,y1 = E.local_coord(P1, prec=3)
            sage: y1^2 == f(x1)
            True
            sage: P2 = E([6,-16])
            sage: x2,y2 = E.local_coord(P2, prec=3)
            sage: y2^2 == f(x2)
            True
            sage: O = E.zero()
            sage: x3,y3 = E.local_coord(O, prec=3)
            sage: y3^2 == f(x3)
            True

        AUTHOR:

        - Jennifer Balakrishnan (2007-12)
        """
        if P[1].is_zero():
            return self.local_coordinates_at_weierstrass(P, prec, name)
        if P[2].is_zero():
            return self.local_coordinates_at_infinity(prec, name)
        return self.local_coordinates_at_nonweierstrass(P, prec, name)

    def monsky_washnitzer_gens(self):
        """
        Return the generators of the special hyperelliptic quotient ring.

        TODO: Should this function be moved to ell_generic and made available over
        more general base rings?

        EXAMPLES::

            sage: Q5 = pAdicField(5,10)
            sage: E = EllipticCurve(Q5,[1,0])
            sage: x,y = E.monsky_washnitzer_gens()
            sage: x^3 + x == y^2
            True
        """
        S = monsky_washnitzer.SpecialHyperellipticQuotientRing(self)
        return S.gens()

    def invariant_differential(self):
        """
        Return `dx/2y`, as an element of the Monsky-Washnitzer cohomology
        of ``self``.

        EXAMPLES::

            sage: Q5 = pAdicField(5,10)
            sage: E = EllipticCurve(Q5,[0,1])
            sage: w = E.invariant_differential(); w
            (-((4+4*5+4*5^2+4*5^3+4*5^4+4*5^5+4*5^6+4*5^7+4*5^8+4*5^9+O(5^10)))*1) dx/2y
            sage: x,y = E.monsky_washnitzer_gens()
            sage: x.diff() == 2*y*w
            True
        """
        S = monsky_washnitzer.SpecialHyperellipticQuotientRing(self)
        MW = monsky_washnitzer.MonskyWashnitzerDifferentialRing(S)
        return MW.invariant_differential()

    def local_analytic_interpolation(self, P, Q):
        """
        For points `P`, `Q` in the same residue disc,
        this constructs an interpolation from `P` to `Q`
        (in homogeneous coordinates) in a power series in
        the local parameter `t`, with precision equal to
        the `p`-adic precision of the underlying ring.

        INPUT:

        - ``P``, ``Q`` -- points on ``self`` in the same residue disc

        OUTPUT:

        Returns a point `X(t) = ( x(t) : y(t) : z(t) )` such that:

        (1) `X(0) = P` and `X(1) = Q` if `P, Q` are not in the infinite disc
        (2) `X(P[0]/P[1]) = P` and `X(Q[0]/Q[1]) = Q` if `P, Q` are in the infinite disc

        EXAMPLES::

            sage: Q7 = pAdicField(7)
            sage: E = EllipticCurve(Q7,[1,1])
            sage: P = E.lift_x(2)
            sage: Q = E.lift_x(9)
            sage: X = E.local_analytic_interpolation(P,Q)
            sage: X[0](1) == Q[0]
            True
            sage: X[0](0) == P[0]
            True

        AUTHORS:

        - Robert Bradshaw (2007-03)
        - Jennifer Balakrishnan (2010-02)
        """
        prec = self.base_ring().precision_cap()
        if not self.is_same_disc(P, Q):
            raise ValueError(f"{P} and {Q} are not in the same residue disc")
        disc = self.residue_disc(P)
        t = PowerSeriesRing(self.base_ring(), "t", prec).gen(0)
        if disc == self.change_ring(self.base_ring().residue_field())(
            0, 1, 0
        ):  # Infinite disc
            x, y = self.local_coordinates_at_infinity(2 * prec)
            return (x * t**3, y * t**3, t**3)
        if disc[1] != 0:  # non-Weierstrass disc
            x = P[0] + t * (Q[0] - P[0])
            pts = self.lift_x(x, all=True)
            if pts[0][1][0] == P[1]:
                return pts[0]
            else:
                return pts[1]
        else:  # Weierstrass disc
            S = self.find_char_zero_weier_point(P)
            x, y = self.local_coord(S, prec)
            a = P[1]
            b = Q[1] - P[1]
            y = a + b * t
            x = x.polynomial()(y).add_bigoh(x.prec())
            return (x, y, 1)

    def weierstrass_points(self):
        """
        Return the Weierstrass points of ``self`` defined over
        ``self.base_ring()``, that is, the point at infinity and those points
        in the support of the divisor of `y`.

        EXAMPLES::

            sage: Q5 = pAdicField(5,10)
            sage: E = EllipticCurve(Q5,[1,0])
            sage: E.weierstrass_points()
            [(0 : 1 + O(5^10) : 0),
             (2 + 5 + 2*5^2 + 5^3 + 3*5^4 + 4*5^5 + 2*5^6 + 3*5^7 + 3*5^9 + O(5^10) : 0 : 1 + O(5^10)),
             (3 + 3*5 + 2*5^2 + 3*5^3 + 5^4 + 2*5^6 + 5^7 + 4*5^8 + 5^9 + O(5^10) : 0 : 1 + O(5^10)),
             (O(5^10) : 0 : 1 + O(5^10))]

            sage: Q7 = pAdicField(7,10)
            sage: E = EllipticCurve(Q7,[1,0])
            sage: E.weierstrass_points()
            [(0 : 1 + O(7^10) : 0), (O(7^10) : 0 : 1 + O(7^10))]
        """
        f, h = self.hyperelliptic_polynomials()
        if not h.is_zero():
            raise NotImplementedError()
        return [self((0, 1, 0))] + [
            self((x, 0, 1)) for x in f.roots(multiplicities=False)
        ]

    def is_in_weierstrass_disc(self, P):
        """
        Check if `P` is in a Weierstrass disc.

        EXAMPLES::

            sage: K = Qp(5,8)
            sage: E = EllipticCurve(K, [-10,9])
            sage: P = E([0,3])
            sage: E.is_in_weierstrass_disc(P)
            False
            sage: Q = E([0,1,0])
            sage: E.is_in_weierstrass_disc(Q)
            True
            sage: S = E([1,0])
            sage: E.is_in_weierstrass_disc(S)
            True
            sage: T = E.lift_x(1+3*5^2); T
            (1 + 3*5^2 + O(5^8) : 3*5 + 4*5^2 + 5^4 + 3*5^5 + 5^6 + O(5^7) : 1 + O(5^8))
            sage: E.is_in_weierstrass_disc(T)
            True

        AUTHOR:

        - Jennifer Balakrishnan (2010-02)
        """
        return not (P[1].valuation() == 0 and P != self(0, 1, 0))

    def is_weierstrass(self, P):
        """
        Check if `P` is a Weierstrass point (i.e., fixed by the hyperelliptic
        involution).

        EXAMPLES::

            sage: K = Qp(5,8)
            sage: E = EllipticCurve(K, [-10,9])
            sage: P = E([0,3])
            sage: E.is_weierstrass(P)
            False
            sage: Q = E([0,1,0])
            sage: E.is_weierstrass(Q)
            True
            sage: S = E([1,0])
            sage: E.is_weierstrass(S)
            True
            sage: T = E.lift_x(1+3*5^2); T
            (1 + 3*5^2 + O(5^8) : 3*5 + 4*5^2 + 5^4 + 3*5^5 + 5^6 + O(5^7) : 1 + O(5^8))
            sage: E.is_weierstrass(T)
            False

        AUTHOR:

        - Jennifer Balakrishnan (2010-02)
        """
        return P[1].is_zero() or P[2].is_zero()

    def find_char_zero_weier_point(self, Q):
        """
        Given `Q` a point on ``self`` in a Weierstrass disc, finds the
        center of the Weierstrass disc (if defined over ``self.base_ring()``).

        EXAMPLES::

            sage: K = Qp(5,8)
            sage: E = EllipticCurve(K, [-10,9])
            sage: Q = E([0,1,0])
            sage: E.find_char_zero_weier_point(Q) == Q
            True

            sage: S = E([1,0])
            sage: T = E.lift_x(1+3*5^2)
            sage: E.find_char_zero_weier_point(T) == S
            True

            sage: P = E([0,3])
            sage: E.find_char_zero_weier_point(P)
            Traceback (most recent call last):
            ...
            ValueError: (0 : 3 + O(5^8) : 1 + O(5^8)) is not in a Weierstrass disc

        AUTHOR:

        - Jennifer Balakrishnan
        """
        if not self.is_in_weierstrass_disc(Q):
            raise ValueError(f"{Q} is not in a Weierstrass disc")
        points = self.weierstrass_points()
        for P in points:
            if self.is_same_disc(P, Q):
                return P

    def residue_disc(self, P):
        """
        Return the residue disc of `P`.

        EXAMPLES::

            sage: K = Qp(5,8)
            sage: E = EllipticCurve(K, [-10,9])
            sage: P = E.lift_x(5); P
            (5 + O(5^9) : 2 + 4*5 + 5^2 + 2*5^3 + 5^4 + 2*5^5 + 2*5^6 + 5^7 + O(5^8) : 1 + O(5^8))
            sage: E.residue_disc(P)
            (0 : 2 : 1)
            sage: E.residue_disc(P) == P.change_ring(GF(5))
            True

        Note that the residue disc can also be computed when
        the coordinates have negative valuation (in which case
        `change_ring` does not work)::

            sage: Q = E.lift_x(5^(-2))
            sage: E.residue_disc(Q)
            (0 : 1 : 0)
            sage: Q.change_ring(GF(5))
            Traceback (most recent call last):
            ...
            ValueError: element must have nonnegative valuation in order to compute residue

        AUTHOR:

        - Jennifer Balakrishnan
        """
        xPv = P[0].valuation()
        yPv = P[1].valuation()
        F = self.base_ring().residue_field()
        try:
            HF = self.change_ring(F)
        except ValueError:
            raise ValueError(
                "The base change of the elliptic curve to the residue field is not well-defined."
                )

        if P == self(0, 1, 0):
            return HF(0, 1, 0)
        elif yPv > 0:
            if xPv > 0:
                return HF(0, 0, 1)
            else:
                return HF(P[0].expansion(0), 0, 1)
        elif yPv == 0:
            if xPv > 0:
                return HF(0, P[1].expansion(0), 1)
            else:
                return HF(P[0].expansion(0), P[1].expansion(0), 1)
        else:
            return HF(0, 1, 0)

    def is_same_disc(self, P, Q):
        """
        Check if `P,Q` are in the same residue disc.

        EXAMPLES::

            sage: Q7 = pAdicField(7,6)
            sage: E = EllipticCurve(Q7,[-16,400])
            sage: P = E.lift_x(4)
            sage: Q = E.lift_x(8)
            sage: R = E.lift_x(11)
            sage: E.is_same_disc(P,Q) or E.is_same_disc(P,-Q)
            False
            sage: E.is_same_disc(P,R) or E.is_same_disc(P,-R)
            True
        """
        return self.residue_disc(P) == self.residue_disc(Q)

    def tiny_integrals(self, F, P, Q):
        r"""
        Evaluate the integrals of `f_i dx/2y` from `P` to `Q` for each `f_i` in `F`
        by formally integrating a power series in a local parameter `t`

        `P` and `Q` MUST be in the same residue disc for this result to make sense.

        INPUT:

        - ``F`` -- list of functions `f_i`
        - ``P`` -- point on ``self``
        - ``Q`` -- point on ``self`` (in the same residue disc as `P`)

        OUTPUT: the integrals `\int_P^Q f_i dx/2y`

        EXAMPLES::

            sage: K = pAdicField(17, 5)
            sage: E = EllipticCurve(K, [-31/3, -2501/108]) # 11a
            sage: P = E(K(14/3), K(11/2))
            sage: TP = E.teichmuller(P);
            sage: x,y = E.monsky_washnitzer_gens()
            sage: E.tiny_integrals([1,x], P, TP) == E.tiny_integrals_on_basis(P, TP)
            True
        """
        x, y, z = self.local_analytic_interpolation(P, Q)  # homogeneous coordinates
        x = x / z
        y = y / z
        dt = x.derivative() / (2 * y)
        integrals = []

        for f in F:
            try:
                f_dt = f(x, y) * dt
            except TypeError:  # if f is a constant, not callable
                f_dt = f * dt
            if x.valuation() != -2:
                I = sum(
                    f_dt[n] / (n + 1) for n in range(f_dt.degree() + 1)
                )  # \int_0^1 f dt
            else:
                If_dt = f_dt.integral().laurent_polynomial()
                I = If_dt(Q[0] / Q[1]) - If_dt(P[0] / P[1])
            integrals.append(I)
        return vector(integrals)

    def tiny_integrals_on_basis(self, P, Q):
        r"""
        Evaluate the integrals `\{\int_P^Q x^i dx/2y \}_{i=0}^{2g-1}`
        by formally integrating a power series in a local parameter `t`.
        `P` and `Q` MUST be in the same residue disc for this result to make sense.

        INPUT:

        - ``P`` -- point on ``self``
        - ``Q`` -- point on ``self`` (in the same residue disc as `P`)

        OUTPUT: the integrals `\{\int_P^Q x^i dx/2y \}_{i=0}^{2g-1}`

        EXAMPLES::

            sage: K = pAdicField(17, 5)
            sage: E = EllipticCurve(K, [-31/3, -2501/108]) # 11a
            sage: P = E(K(14/3), K(11/2))
            sage: TP = E.teichmuller(P);
            sage: E.tiny_integrals_on_basis(P, TP)
            (17 + 14*17^2 + 17^3 + 8*17^4 + O(17^5), 16*17 + 5*17^2 + 8*17^3 + 14*17^4 + O(17^5))
        """
        if P == Q:
            V = VectorSpace(self.base_ring(), 2)
            return V(0)
        R = PolynomialRing(self.base_ring(), ["x", "y"])
        x, _ = R.gens()
        return self.tiny_integrals([x**i for i in range(2)], P, Q)

    def teichmuller(self, P):
        r"""
        Find a Teichm\:uller point in the same residue class of `P`.

        Because this lift of frobenius acts as `x \mapsto x^p`,
        take the Teichmuller lift of `x` and then find a matching `y`
        from that.

        EXAMPLES::

            sage: K = pAdicField(7, 5)
            sage: E = EllipticCurve(K, [-31/3, -2501/108]) # 11a
            sage: P = E(K(14/3), K(11/2))
            sage: E.frobenius(P) == P
            False
            sage: TP = E.teichmuller(P); TP
            (0 : 2 + 3*7 + 3*7^2 + 3*7^4 + O(7^5) : 1 + O(7^5))
            sage: E.frobenius(TP) == TP
            True
            sage: (TP[0] - P[0]).valuation() > 0, (TP[1] - P[1]).valuation() > 0
            (True, True)
        """
        K = P[0].parent()
        x = K.teichmuller(P[0])
        pts = self.lift_x(x, all=True)
        if (pts[0][1] - P[1]).valuation() > 0:
            return pts[0]
        else:
            return pts[1]

    def coleman_integrals_on_basis(self, P, Q, algorithm=None):
        r"""
        Compute the Coleman integrals `\{\int_P^Q x^i dx/2y \}_{i=0}^{2g-1}`.

        INPUT:

        - ``P`` -- point on ``self``
        - ``Q`` -- point on ``self``
        - ``algorithm`` -- ``None`` (default, uses Frobenius) or teichmuller
          (uses Teichmuller points)

        OUTPUT: the Coleman integrals `\{\int_P^Q x^i dx/2y \}_{i=0}^{2g-1}`

        EXAMPLES::

            sage: K = Qp(5,8)
            sage: E = EllipticCurve(K, [-10,9])
            sage: S = E(1,0)
            sage: P = E(0,3)
            sage: T = E(0,1,0)
            sage: Q = E.lift_x(5^-2)
            sage: R = E.lift_x(4*5^-2)
            sage: E.coleman_integrals_on_basis(S,P)
            (2*5^2 + 5^4 + 5^5 + 3*5^6 + 3*5^7 + 2*5^8 + O(5^9), 5 + 2*5^2 + 4*5^3 + 2*5^4 + 3*5^6 + 4*5^7 + 2*5^8 + O(5^9))
            sage: E.coleman_integrals_on_basis(T,P)
            (2*5^2 + 5^4 + 5^5 + 3*5^6 + 3*5^7 + 2*5^8 + O(5^9), 5 + 2*5^2 + 4*5^3 + 2*5^4 + 3*5^6 + 4*5^7 + 2*5^8 + O(5^9))
            sage: E.coleman_integrals_on_basis(P,S) == -E.coleman_integrals_on_basis(S,P)
            True
            sage: E.coleman_integrals_on_basis(S,Q)
            (5 + O(5^4), 4*5^-1 + 4 + 4*5 + 4*5^2 + O(5^3))
            sage: E.coleman_integrals_on_basis(Q,R)
            (5 + 2*5^2 + 2*5^3 + 2*5^4 + 3*5^5 + 3*5^6 + 3*5^7 + 5^8 + O(5^9), 3*5^-1 + 2*5^4 + 5^5 + 2*5^6 + O(5^7))
            sage: E.coleman_integrals_on_basis(S,R) == E.coleman_integrals_on_basis(S,Q) + E.coleman_integrals_on_basis(Q,R)
            True
            sage: E.coleman_integrals_on_basis(T,T)
            (0, 0)
            sage: E.coleman_integrals_on_basis(S,T)
            (0, 0)

        AUTHORS:

        - Robert Bradshaw (2007-03): non-Weierstrass points
        - Jennifer Balakrishnan and Robert Bradshaw (2010-02): Weierstrass points
        """
        K = self.base_ring()
        p = K.prime()
        prec = K.precision_cap()
        dim = 2
        V = VectorSpace(K, dim)
        # if P or Q is Weierstrass, use the Frobenius algorithm
        if self.is_weierstrass(P):
            if self.is_weierstrass(Q):
                return V(0)
            else:
                PP = None
                QQ = Q
                TP = None
                TQ = self.frobenius(Q)
        elif self.is_weierstrass(Q):
            PP = P
            QQ = None
            TQ = None
            TP = self.frobenius(P)
        elif self.is_same_disc(P, Q):
            return self.tiny_integrals_on_basis(P, Q)
        elif algorithm == "teichmuller":
            PP = TP = self.teichmuller(P)
            QQ = TQ = self.teichmuller(Q)
        else:
            TP = self.frobenius(P)
            TQ = self.frobenius(Q)
            PP, QQ = P, Q
        if TP is None:
            P_to_TP = V(0)
        else:
            if TP is not None:
                TPv = (TP[0] / TP[1]).valuation()
                xTPv = TP[0].valuation()
            else:
                xTPv = TPv = +Infinity
            if TQ is not None:
                TQv = (TQ[0] / TQ[1]).valuation()
                xTQv = TQ[0].valuation()
            else:
                xTQv = TQv = +Infinity
            offset = (2 - 1) * max(TPv, TQv)
            if offset == +Infinity:
                offset = (2 - 1) * min(TPv, TQv)
            if (
                offset > prec
                and (xTPv < 0 or xTQv < 0)
                and (
                    self.residue_disc(P) == self.change_ring(GF(p))(0, 1, 0)
                    or self.residue_disc(Q) == self.change_ring(GF(p))(0, 1, 0)
                )
            ):
                newprec = offset + prec
                K = pAdicField(p, newprec)
                A = PolynomialRing(RationalField(), "x")
                f = A(self.hyperelliptic_polynomials()[0])
                self = EllipticCurve(f).change_ring(K)
                xP = P[0]
                xPv = xP.valuation()
                xPnew = K(sum(c * p ** (xPv + i) for i, c in enumerate(xP.expansion())))
                PP = P = self.lift_x(xPnew)
                TP = self.frobenius(P)
                xQ = Q[0]
                xQv = xQ.valuation()
                xQnew = K(sum(c * p ** (xQv + i) for i, c in enumerate(xQ.expansion())))
                QQ = Q = self.lift_x(xQnew)
                TQ = self.frobenius(Q)
                V = VectorSpace(K, dim)
            P_to_TP = V(self.tiny_integrals_on_basis(P, TP))
        if TQ is None:
            TQ_to_Q = V(0)
        else:
            TQ_to_Q = V(self.tiny_integrals_on_basis(TQ, Q))
        try:
            M_frob, forms = self._frob_calc
        except AttributeError:
            M_frob, forms = self._frob_calc = (
                monsky_washnitzer.matrix_of_frobenius_hyperelliptic(self)
            )
        R = forms[0].base_ring()
        try:
            if PP is None:
                L = [-ff(R(QQ[0]), R(QQ[1])) for ff in forms]  ##changed
            elif QQ is None:
                L = [ff(R(PP[0]), R(PP[1])) for ff in forms]
            else:
                L = [ff(R(PP[0]), R(PP[1])) - ff(R(QQ[0]), R(QQ[1])) for ff in forms]
        except ValueError:
            forms = [ff.change_ring(self.base_ring()) for ff in forms]
            if PP is None:
                L = [-ff(QQ[0], QQ[1]) for ff in forms]  ##changed
            elif QQ is None:
                L = [ff(PP[0], PP[1]) for ff in forms]
            else:
                L = [ff(PP[0], PP[1]) - ff(QQ[0], QQ[1]) for ff in forms]
        b = V(L)
        if PP is None:
            b -= TQ_to_Q
        elif QQ is None:
            b -= P_to_TP
        elif algorithm != "teichmuller":
            b -= P_to_TP + TQ_to_Q
        M_sys = matrix(K, M_frob).transpose() - 1
        TP_to_TQ = M_sys ** (-1) * b
        if algorithm == "teichmuller":
            return P_to_TP + TP_to_TQ + TQ_to_Q
        else:
            return TP_to_TQ

    def coleman_integral(self, w, P, Q, algorithm="None"):
        r"""
        Return the Coleman integral `\int_P^Q w`.

        INPUT:

        - ``w`` -- differential (if one of P,Q is Weierstrass, w must be odd)
        - ``P`` -- point on ``self``
        - ``Q`` -- point on ``self``
        - ``algorithm`` -- ``None`` (default, uses Frobenius) or teichmuller
          (uses Teichmuller points)

        OUTPUT: the Coleman integral `\int_P^Q w`

        EXAMPLES:

        A simple example, integrating dx::

            sage: K = Qp(5,10)
            sage: E = EllipticCurve(K, [-4,4])
            sage: P = E(2, 2)
            sage: Q = E.teichmuller(P)
            sage: x, y = E.monsky_washnitzer_gens()
            sage: E.coleman_integral(x.diff(), P, Q)
            5 + 2*5^2 + 5^3 + 3*5^4 + 4*5^5 + 2*5^6 + 3*5^7 + 3*5^9 + O(5^10)
            sage: Q[0] - P[0]
            5 + 2*5^2 + 5^3 + 3*5^4 + 4*5^5 + 2*5^6 + 3*5^7 + 3*5^9 + O(5^10)

        Another example::

            sage: K = Qp(7,10)
            sage: E = EllipticCurve(K, [0,8,0,-9,0])
            sage: _, forms = monsky_washnitzer.matrix_of_frobenius_hyperelliptic(E)
            sage: w = E.invariant_differential()
            sage: x,y = E.monsky_washnitzer_gens()
            sage: f = forms[0]
            sage: S = E(9,36)
            sage: Q = E.teichmuller(S)
            sage: P = E(-1,4)
            sage: b = x*w*w._coeff.parent()(f)
            sage: E.coleman_integral(b,P,Q)
            7 + 7^2 + 4*7^3 + 5*7^4 + 3*7^5 + 7^6 + 5*7^7 + 3*7^8 + 4*7^9 + 4*7^10 + O(7^11)

        ::

            sage: K = Qp(5,8)
            sage: E = EllipticCurve(K, [0,1])
            sage: w = E.invariant_differential()
            sage: P = E(0,1)
            sage: Q = E(5, 1 + 3*5^3 + 2*5^4 + 2*5^5 + 3*5^7)
            sage: x,y = E.monsky_washnitzer_gens()
            sage: (2*y*w).coleman_integral(P,Q)
            5 + O(5^9)
            sage: xloc,yloc,zloc = E.local_analytic_interpolation(P,Q)
            sage: I2 = (xloc.derivative()/(2*yloc)).integral()
            sage: I2.polynomial()(1) - I2(0)
            3*5 + 2*5^2 + 2*5^3 + 5^4 + 4*5^6 + 5^7 + O(5^9)
            sage: E.coleman_integral(w,P,Q)
            3*5 + 2*5^2 + 2*5^3 + 5^4 + 4*5^6 + 5^7 + O(5^9)

        Integrals involving Weierstrass points::

            sage: K = Qp(5,8)
            sage: E = EllipticCurve(K, [-10,9])
            sage: S = E(1,0)
            sage: P = E(0,3)
            sage: negP = E(0,-3)
            sage: T = E(0,1,0)
            sage: w = E.invariant_differential()
            sage: x,y = E.monsky_washnitzer_gens()
            sage: E.coleman_integral(w*x^3,S,T)
            0
            sage: E.coleman_integral(w*x^3,T,S)
            0
            sage: E.coleman_integral(w,S,P)
            2*5^2 + 5^4 + 5^5 + 3*5^6 + 3*5^7 + 2*5^8 + O(5^9)
            sage: E.coleman_integral(w,T,P)
            2*5^2 + 5^4 + 5^5 + 3*5^6 + 3*5^7 + 2*5^8 + O(5^9)
            sage: E.coleman_integral(w*x^3,T,P)
            5^2 + 2*5^3 + 3*5^6 + 3*5^7 + O(5^8)
            sage: E.coleman_integral(w*x^3,S,P)
            5^2 + 2*5^3 + 3*5^6 + 3*5^7 + O(5^8)
            sage: E.coleman_integral(w, P, negP, algorithm='teichmuller')
            5^2 + 4*5^3 + 2*5^4 + 2*5^5 + 3*5^6 + 2*5^7 + 4*5^8 + O(5^9)
            sage: E.coleman_integral(w, P, negP)
            5^2 + 4*5^3 + 2*5^4 + 2*5^5 + 3*5^6 + 2*5^7 + 4*5^8 + O(5^9)

        AUTHORS:

        - Robert Bradshaw (2007-03)
        - Kiran Kedlaya (2008-05)
        - Jennifer Balakrishnan (2010-02)
        """
        K = self.base_ring()
        prec = K.precision_cap()
        S = monsky_washnitzer.SpecialHyperellipticQuotientRing(self, K)
        MW = monsky_washnitzer.MonskyWashnitzerDifferentialRing(S)
        w = MW(w)
        f, vec = w.reduce_fast()
        basis_values = self.coleman_integrals_on_basis(P, Q, algorithm)
        dim = len(basis_values)
        x, y = self.local_coordinates_at_infinity(2 * prec)
        if self.is_weierstrass(P):
            if self.is_weierstrass(Q):
                return 0
            elif f == 0:
                return sum([vec[i] * basis_values[i] for i in range(dim)])
            elif (
                w._coeff(x, -y) * x.derivative() / (-2 * y)
                + w._coeff(x, y) * x.derivative() / (2 * y)
                == 0
            ):
                return (
                    self.coleman_integral(
                        w, self(Q[0], -Q[1]), self(Q[0], Q[1]), algorithm
                    )
                    / 2
                )
            else:
                raise ValueError(
                    "The differential is not odd: use coleman_integral_from_weierstrass_via_boundary"
                )

        elif self.is_weierstrass(Q):
            if f == 0:
                return sum([vec[i] * basis_values[i] for i in range(dim)])
            elif (
                w._coeff(x, -y) * x.derivative() / (-2 * y)
                + w._coeff(x, y) * x.derivative() / (2 * y)
                == 0
            ):
                return (
                    -self.coleman_integral(
                        w, self(P[0], -P[1]), self(P[0], P[1]), algorithm
                    )
                    / 2
                )
            else:
                raise ValueError(
                    "The differential is not odd: use coleman_integral_from_weierstrass_via_boundary"
                )
        else:
            return (
                f(Q[0], Q[1])
                - f(P[0], P[1])
                + sum([vec[i] * basis_values[i] for i in range(dim)])
            )  # this is just a dot product...

    def curve_over_ram_extn(self, deg):
        r"""
        Return ``self`` over `\QQ_p(p^(1/deg))`.

        INPUT:

        - ``deg`` -- the degree of the ramified extension

        OUTPUT: ``self`` over `\QQ_p(p^(1/deg))`

        EXAMPLES::

            sage: K = Qp(5,8)
            sage: E = EllipticCurve(K, [0,1])
            sage: E.curve_over_ram_extn(2)
            Elliptic Curve defined by y^2 = x^3 + (1+O(a^16)) over 5-adic Eisenstein Extension Field in a defined by x^2 - 5

        AUTHOR:

        - Jennifer Balakrishnan
        """
        K = self.base_ring()
        p = K.prime()
        A = PolynomialRing(QQ, "x")
        x = A.gen()
        J = K.extension(x**deg - p, names="a")
        EJ = self.change_ring(J)
        self._curve_over_ram_extn = EJ
        self._curve_over_ram_extn._curve_over_Qp = self
        return EJ

    def get_boundary_point(self, curve_over_extn, P):
        """
        Given ``self`` over an extension field, find a point in the disc of `P`
        near the boundary.

        INPUT:

        - ``curve_over_extn`` -- ``self`` over a totally ramified extension
        - ``P`` -- Weierstrass point

        OUTPUT: a point in the disc of `P` near the boundary

        EXAMPLES::

            sage: K = Qp(3,6)
            sage: E = EllipticCurve(K,[-10,9])
            sage: P = E(1,0)
            sage: J.<a> = K.extension(x^30-3)
            sage: EJ  = E.change_ring(J)
            sage: S = E.get_boundary_point(EJ,P)
            sage: S
            (1 + 2*a^2 + 2*a^6 + 2*a^18 + a^32 + a^34 + a^36 + 2*a^38 + 2*a^40 + a^42 + 2*a^44 + a^48 + 2*a^50 + 2*a^52 + a^54 + a^56 + 2*a^60 + 2*a^62 + a^70 + 2*a^72 + a^76 + 2*a^78 + a^82 + a^88 + a^96 + 2*a^98 + 2*a^102 + a^104 + 2*a^106 + a^108 + 2*a^110 + a^112 + 2*a^116 + a^126 + 2*a^130 + 2*a^132 + a^144 + 2*a^148 + 2*a^150 + a^152 + 2*a^154 + a^162 + a^164 + a^166 + a^168 + a^170 + a^176 + a^178 + O(a^180) : a + O(a^180) : 1 + O(a^180))

        AUTHOR:

        - Jennifer Balakrishnan
        """
        J = curve_over_extn.base_ring()
        a = J.gen()
        prec2 = J.precision_cap()
        x, y = self.local_coord(P, prec2)
        return curve_over_extn(x(a), y(a))

    def P_to_S(self, P, S):
        r"""
        Given a finite Weierstrass point `P` and a point `S` in the same disc,
        compute the Coleman integrals `\{\int_P^S x^i dx/2y \}_{i=0}^{2g-1}`.

        INPUT:

        - ``P`` -- finite Weierstrass point
        - ``S`` -- point in disc of `P`

        OUTPUT: Coleman integrals `\{\int_P^S x^i dx/2y \}_{i=0}^{2g-1}`

        EXAMPLES::

            sage: K = Qp(5,4)
            sage: E = EllipticCurve(K, [-10,9])
            sage: P = E(1,0)
            sage: EJ = E.curve_over_ram_extn(10)
            sage: S = E.get_boundary_point(EJ,P)
            sage: E.P_to_S(P, S)
            (2*a + 4*a^3 + 2*a^11 + 4*a^13 + 2*a^17 + 2*a^19 + a^21 + 4*a^23 + a^25 + 2*a^27 + 2*a^29 + 3*a^31 + 4*a^33 + O(a^35), a^-5 + 2*a + 2*a^3 + a^7 + 3*a^11 + a^13 + 3*a^15 + 3*a^17 + 2*a^19 + 4*a^21 + 4*a^23 + 4*a^25 + 2*a^27 + a^29 + a^31 + O(a^33))

        AUTHOR:

        - Jennifer Balakrishnan
        """
        prec = self.base_ring().precision_cap()
        deg = (S[0]).parent().defining_polynomial().degree()
        prec2 = prec * deg
        x, y = self.local_coord(P, prec2)
        integrals = [((x**k * x.derivative() / (2 * y)).integral()) for k in range(2)]
        val = [I(S[1]) for I in integrals]
        return vector(val)

    def coleman_integral_P_to_S(self, w, P, S):
        r"""
        Given a finite Weierstrass point `P` and a point `S`
        in the same disc, compute the Coleman integral `\int_P^S w`.

        INPUT:

        - ``w`` -- differential
        - ``P`` -- Weierstrass point
        - ``S`` -- point in the same disc of `P` (S is defined over an extension
          of `\QQ_p`; coordinates of S are given in terms of uniformizer `a`)

        OUTPUT: Coleman integral `\int_P^S w` in terms of `a`

        EXAMPLES::

            sage: K = Qp(5,4)
            sage: E = EllipticCurve(K, [-10,9])
            sage: P = E(1,0)
            sage: J.<a> = K.extension(x^10-5)
            sage: EJ  = E.change_ring(J)
            sage: S = E.get_boundary_point(EJ,P)
            sage: x,y = E.monsky_washnitzer_gens()
            sage: S[0]-P[0] == E.coleman_integral_P_to_S(x.diff(),P,S)
            True
            sage: E.coleman_integral_P_to_S(E.invariant_differential(),P,S) == E.P_to_S(P,S)[0]
            True

        AUTHOR:

        - Jennifer Balakrishnan
        """
        prec = self.base_ring().precision_cap()
        deg = S[0].parent().defining_polynomial().degree()
        prec2 = prec * deg
        x, y = self.local_coord(P, prec2)
        int_sing = (w.coeff()(x, y) * x.derivative() / (2 * y)).integral()
        return int_sing(S[1])

    def S_to_Q(self, S, Q):
        r"""
        Given `S` a point on ``self`` over an extension field, compute the
        Coleman integrals `\{\int_S^Q x^i dx/2y \}_{i=0}^{2g-1}`.

        **one should be able to feed `S,Q` into coleman_integral,
        but currently that segfaults**

        INPUT:

        - ``S`` -- a point with coordinates in an extension of `\QQ_p` (with unif. a)
        - ``Q`` -- a non-Weierstrass point defined over `\QQ_p`

        OUTPUT:

        The Coleman integrals `\{\int_S^Q x^i dx/2y \}_{i=0}^{2g-1}` in terms of `a`.

        EXAMPLES::

            sage: K = Qp(5,6)
            sage: E = EllipticCurve(K, [-10,9])
            sage: J.<a> = K.extension(x^20-5)
            sage: EJ  = E.change_ring(J)
            sage: w = E.invariant_differential()
            sage: x,y = E.monsky_washnitzer_gens()
            sage: P = E(1,0)
            sage: Q = E(0,3)
            sage: S = E.get_boundary_point(EJ,P)
            sage: P_to_S = E.P_to_S(P,S)
            sage: S_to_Q = EJ.S_to_Q(S,Q)
            sage: P_to_S + S_to_Q
            (2*a^40 + a^80 + a^100 + O(a^105), a^20 + 2*a^40 + 4*a^60 + 2*a^80 + O(a^103))
            sage: E.coleman_integrals_on_basis(P,Q)
            (2*5^2 + 5^4 + 5^5 + 3*5^6 + O(5^7), 5 + 2*5^2 + 4*5^3 + 2*5^4 + 5^6 + O(5^7))

        AUTHOR:

        - Jennifer Balakrishnan
        """
        FS = self.frobenius(S)
        FS = (FS[0], FS[1])
        FQ = self.frobenius(Q)

        try:
            M_frob, forms = self._frob_calc
        except AttributeError:
            M_frob, forms = self._frob_calc = (
                monsky_washnitzer.matrix_of_frobenius_hyperelliptic(self)
            )
        try:
            HJ = self._curve_over_ram_extn
            K = HJ.base_ring()
        except AttributeError:
            HJ = S.scheme()
            K = self.base_ring()

        prec2 = K.precision_cap()
        p = K.prime()
        dim = 2
        V = VectorSpace(K, dim)
        if S == FS:
            S_to_FS = V(dim * [0])
        else:
            P = self(ZZ(FS[0].expansion(0)), ZZ(FS[1].expansion(0)))
            x, y = self.local_coord(P, prec2)
            integrals = [
                (x**i * x.derivative() / (2 * y)).integral() for i in range(dim)
            ]
            S_to_FS = vector(
                [I.polynomial()(FS[1]) - I.polynomial()(S[1]) for I in integrals]
            )
        if HJ(Q[0], Q[1]) == HJ(FQ):
            FQ_to_Q = V(dim * [0])
        else:
            FQ_to_Q = V(self.tiny_integrals_on_basis(FQ, Q))
        try:
            L = [f(K(S[0]), K(S[1])) - f(K(Q[0]), K(Q[1])) for f in forms]
        except ValueError:
            forms = [f.change_ring(K) for f in forms]
            L = [f(S[0], S[1]) - f(Q[0], Q[1]) for f in forms]
        b = V(L)
        M_sys = matrix(K, M_frob).transpose() - 1
        B = ~M_sys
        vv = min(c.valuation() for c in B.list())
        B = (p ** (-vv) * B).change_ring(K)
        B = p ** (vv) * B
        return B * (b - S_to_FS - FQ_to_Q)

    def coleman_integral_S_to_Q(self, w, S, Q):
        r"""
        Compute the Coleman integral `\int_S^Q w`.

        **one should be able to feed `S,Q` into coleman_integral,
        but currently that segfaults**

        INPUT:

        - ``w`` -- a differential
        - ``S`` -- a point with coordinates in an extension of `\QQ_p`
        - ``Q`` -- a non-Weierstrass point defined over `\QQ_p`

        OUTPUT: the Coleman integral `\int_S^Q w`

        EXAMPLES::

            sage: K = Qp(5,6)
            sage: E = EllipticCurve(K, [-10,9])
            sage: J.<a> = K.extension(x^20-5)
            sage: EJ = E.change_ring(J)
            sage: x,y = E.monsky_washnitzer_gens()
            sage: P = E(1,0)
            sage: Q = E(0,3)
            sage: S = E.get_boundary_point(EJ,P)
            sage: P_to_S = E.coleman_integral_P_to_S(y.diff(),P,S)
            sage: S_to_Q = EJ.coleman_integral_S_to_Q(y.diff(),S,Q)
            sage: P_to_S + S_to_Q
            3 + O(a^119)
            sage: E.coleman_integral(y.diff(),P,Q)
            3 + O(5^6)

        AUTHOR:

        - Jennifer Balakrishnan
        """
        K = self.base_ring()
        R = monsky_washnitzer.SpecialHyperellipticQuotientRing(self, K)
        MW = monsky_washnitzer.MonskyWashnitzerDifferentialRing(R)
        w = MW(w)
        f, vec = w.reduce_fast()

        const = f(Q[0], Q[1]) - f(S[0], S[1])
        if vec == vector([0, 0]):
            return const
        else:
            basis_values = self.S_to_Q(S, Q)
            dim = len(basis_values)
            dot = sum([vec[i] * basis_values[i] for i in range(dim)])
            return const + dot

    def coleman_integral_from_weierstrass_via_boundary(self, w, P, Q, d):
        r"""
        Compute the Coleman integral `\int_P^Q w` via a boundary point
        in the disc of `P`, defined over a degree `d` extension

        INPUT:

        - ``w`` -- a differential
        - ``P`` -- a Weierstrass point
        - ``Q`` -- a non-Weierstrass point
        - ``d`` -- degree of extension where coordinates of boundary point lie

        OUTPUT:

        the Coleman integral `\int_P^Q w`, written in terms of the uniformizer
        `a` of the degree `d` extension

        EXAMPLES::

            sage: K = Qp(5,6)
            sage: E = EllipticCurve(K, [-10,9])
            sage: P = E(1,0)
            sage: Q = E(0,3)
            sage: x,y = E.monsky_washnitzer_gens()
            sage: E.coleman_integral_from_weierstrass_via_boundary(y.diff(),P,Q,20)
            3 + O(a^119)
            sage: E.coleman_integral(y.diff(),P,Q)
            3 + O(5^6)
            sage: w = E.invariant_differential()
            sage: E.coleman_integral_from_weierstrass_via_boundary(w,P,Q,20)
            2*a^40 + a^80 + a^100 + O(a^105)
            sage: E.coleman_integral(w,P,Q)
            2*5^2 + 5^4 + 5^5 + 3*5^6 + O(5^7)

        AUTHOR:

        - Jennifer Balakrishnan
        """
        HJ = self.curve_over_ram_extn(d)
        S = self.get_boundary_point(HJ, P)
        P_to_S = self.coleman_integral_P_to_S(w, P, S)
        S_to_Q = HJ.coleman_integral_S_to_Q(w, S, Q)
        return P_to_S + S_to_Q
