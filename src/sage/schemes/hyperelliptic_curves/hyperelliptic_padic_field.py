# sage.doctest: needs sage.rings.padics
"""
Hyperelliptic curves over a `p`-adic field
"""
# ****************************************************************************
#       Copyright (C) 2007 Robert Bradshaw <robertwb@math.washington.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.matrix.constructor import matrix
from sage.misc.lazy_import import lazy_import
from sage.modules.free_module import VectorSpace
from sage.modules.free_module_element import vector
from sage.rings.finite_rings.finite_field_constructor import FiniteField as GF
from sage.rings.infinity import Infinity
from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.power_series_ring import PowerSeriesRing
from sage.rings.rational_field import QQ
from sage.rings.rational_field import RationalField
from sage.rings.real_mpfr import RR
from sage.schemes.curves.projective_curve import ProjectivePlaneCurve_field

lazy_import("sage.functions.log", "log")
lazy_import("sage.rings.padics.factory", "Qp", as_="pAdicField")

from . import hyperelliptic_generic


class HyperellipticCurve_padic_field(hyperelliptic_generic.HyperellipticCurve_generic,
                                     ProjectivePlaneCurve_field):

    # The functions below were prototyped at the 2007 Arizona Winter School by
    # Robert Bradshaw and Ralf Gerkmann, working with Miljan Brakovevic and
    # Kiran Kedlaya
    # All of the below is with respect to the Monsky Washnitzer cohomology.

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
        (2) `X(P[0]^g/P[1]) = P` and `X(Q[0]^g/Q[1]) = Q` if `P, Q` are in the infinite disc

        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: H = HyperellipticCurve(x^3-10*x+9)
            sage: K = Qp(5,8)
            sage: HK = H.change_ring(K)

        A non-Weierstrass disc::

            sage: P = HK(0,3)
            sage: Q = HK(5, 3 + 3*5^2 + 2*5^3 + 3*5^4 + 2*5^5 + 2*5^6 + 3*5^7 + O(5^8))
            sage: x,y,z, = HK.local_analytic_interpolation(P,Q)
            sage: x(0) == P[0], x(1) == Q[0], y(0) == P[1], y.polynomial()(1) == Q[1]
            (True, True, True, True)

        A finite Weierstrass disc::

            sage: P = HK.lift_x(1 + 2*5^2)
            sage: Q = HK.lift_x(1 + 3*5^2)
            sage: x,y,z = HK.local_analytic_interpolation(P,Q)
            sage: x(0) == P[0], x.polynomial()(1) == Q[0], y(0) == P[1], y(1) == Q[1]
            (True, True, True, True)

        The infinite disc::

            sage: P = HK.lift_x(5^-2)
            sage: Q = HK.lift_x(4*5^-2)
            sage: x,y,z = HK.local_analytic_interpolation(P,Q)
            sage: x = x/z
            sage: y = y/z
            sage: x(P[0]/P[1]) == P[0]
            True
            sage: x(Q[0]/Q[1]) == Q[0]
            True
            sage: y(P[0]/P[1]) == P[1]
            True
            sage: y(Q[0]/Q[1]) == Q[1]
            True

        An error if points are not in the same disc::

            sage: x,y,z = HK.local_analytic_interpolation(P,HK(1,0))
            Traceback (most recent call last):
            ...
            ValueError: (5^-2 + O(5^6) : 4*5^-3 + 4*5^-2 + 4*5^-1 + 4 + 4*5 + 3*5^3 + 5^4 + O(5^5) : 1 + O(5^8)) and (1 + O(5^8) : 0 : 1 + O(5^8)) are not in the same residue disc

        TESTS:

        Check that :issue:`26005` is fixed::

            sage: L = Qp(5, 100)
            sage: HL = H.change_ring(L)
            sage: P = HL.lift_x(1 + 2*5^2)
            sage: Q = HL.lift_x(1 + 3*5^2)
            sage: x,y,z=HL.local_analytic_interpolation(P, Q)
            sage: x.polynomial().degree()
            98

        AUTHORS:

        - Robert Bradshaw (2007-03)
        - Jennifer Balakrishnan (2010-02)
        """
        prec = self.base_ring().precision_cap()
        if not self.is_same_disc(P,Q):
            raise ValueError("%s and %s are not in the same residue disc" % (P,Q))
        disc = self.residue_disc(P)
        t = PowerSeriesRing(self.base_ring(), 't', prec).gen(0)
        if disc == self.change_ring(self.base_ring().residue_field())(0,1,0): # Infinite disc
            x,y = self.local_coordinates_at_infinity(2*prec)
            g = self.genus()
            return (x*t**(2*g+1),y*t**(2*g+1),t**(2*g+1))
        if disc[1] != 0: # non-Weierstrass disc
            x = P[0]+t*(Q[0]-P[0])
            pts = self.lift_x(x, all=True)
            if pts[0][1][0] == P[1]:
                return pts[0]
            else:
                return pts[1]
        else: # Weierstrass disc
            S = self.find_char_zero_weier_point(P)
            x,y = self.local_coord(S, prec)
            a = P[1]
            b = Q[1] - P[1]
            y = a + b*t
            x = x.polynomial()(y).add_bigoh(x.prec())
            return (x, y, 1)

    def weierstrass_points(self):
        """
        Return the Weierstrass points of ``self`` defined over
        ``self.base_ring()``, that is, the point at infinity and those points
        in the support of the divisor of `y`.

        EXAMPLES::

            sage: K = pAdicField(11, 5)
            sage: x = polygen(K)
            sage: C = HyperellipticCurve(x^5 + 33/16*x^4 + 3/4*x^3 + 3/8*x^2 - 1/4*x + 1/16)
            sage: C.weierstrass_points()
            [(0 : 1 + O(11^5) : 0), (7 + 10*11 + 4*11^3 + O(11^5) : 0 : 1 + O(11^5))]
        """
        f, h = self.hyperelliptic_polynomials()
        if h != 0:
            raise NotImplementedError()
        return [self((0,1,0))] + [self((x, 0, 1)) for x in f.roots(multiplicities=False)]

    def is_in_weierstrass_disc(self, P):
        """
        Check if `P` is in a Weierstrass disc.

        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: H = HyperellipticCurve(x^3-10*x+9)
            sage: K = Qp(5,8)
            sage: HK = H.change_ring(K)
            sage: P = HK(0,3)
            sage: HK.is_in_weierstrass_disc(P)
            False
            sage: Q = HK(0,1,0)
            sage: HK.is_in_weierstrass_disc(Q)
            True
            sage: S = HK(1,0)
            sage: HK.is_in_weierstrass_disc(S)
            True
            sage: T = HK.lift_x(1+3*5^2); T
            (1 + 3*5^2 + O(5^8) : 3*5 + 4*5^2 + 5^4 + 3*5^5 + 5^6 + O(5^7) : 1 + O(5^8))
            sage: HK.is_in_weierstrass_disc(T)
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

            sage: R.<x> = QQ['x']
            sage: H = HyperellipticCurve(x^3-10*x+9)
            sage: K = Qp(5,8)
            sage: HK = H.change_ring(K)
            sage: P = HK(0,3)
            sage: HK.is_weierstrass(P)
            False
            sage: Q = HK(0,1,0)
            sage: HK.is_weierstrass(Q)
            True
            sage: S = HK(1,0)
            sage: HK.is_weierstrass(S)
            True
            sage: T = HK.lift_x(1+3*5^2); T
            (1 + 3*5^2 + O(5^8) : 3*5 + 4*5^2 + 5^4 + 3*5^5 + 5^6 + O(5^7) : 1 + O(5^8))
            sage: HK.is_weierstrass(T)
            False

        AUTHOR:

        - Jennifer Balakrishnan (2010-02)
        """
        return (P[1] == 0 or P[2] == 0)

    def find_char_zero_weier_point(self, Q):
        """
        Given `Q` a point on ``self`` in a Weierstrass disc, finds the
        center of the Weierstrass disc (if defined over ``self.base_ring()``).

        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: H = HyperellipticCurve(x^3-10*x+9)
            sage: K = Qp(5,8)
            sage: HK = H.change_ring(K)
            sage: P = HK.lift_x(1 + 2*5^2)
            sage: Q = HK.lift_x(5^-2)
            sage: S = HK(1,0)
            sage: T = HK(0,1,0)
            sage: HK.find_char_zero_weier_point(P)
            (1 + O(5^8) : 0 : 1 + O(5^8))
            sage: HK.find_char_zero_weier_point(Q)
            (0 : 1 + O(5^8) : 0)
            sage: HK.find_char_zero_weier_point(S)
            (1 + O(5^8) : 0 : 1 + O(5^8))
            sage: HK.find_char_zero_weier_point(T)
            (0 : 1 + O(5^8) : 0)

        AUTHOR:

        - Jennifer Balakrishnan
        """
        if not self.is_in_weierstrass_disc(Q):
            raise ValueError("%s is not in a Weierstrass disc" % Q)
        points = self.weierstrass_points()
        for P in points:
            if self.is_same_disc(P,Q):
                return P

    def residue_disc(self, P):
        """
        Give the residue disc of `P`.

        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: H = HyperellipticCurve(x^3-10*x+9)
            sage: K = Qp(5,8)
            sage: HK = H.change_ring(K)
            sage: P = HK.lift_x(1 + 2*5^2)
            sage: HK.residue_disc(P)
            (1 : 0 : 1)
            sage: Q = HK(0,3)
            sage: HK.residue_disc(Q)
            (0 : 3 : 1)
            sage: S = HK.lift_x(5^-2)
            sage: HK.residue_disc(S)
            (0 : 1 : 0)
            sage: T = HK(0,1,0)
            sage: HK.residue_disc(T)
            (0 : 1 : 0)

        AUTHOR:

        - Jennifer Balakrishnan
        """
        xPv = P[0].valuation()
        yPv = P[1].valuation()
        F = self.base_ring().residue_field()
        HF = self.change_ring(F)
        if P == self(0,1,0):
            return HF(0,1,0)
        elif yPv > 0:
            if xPv > 0:
                return HF(0,0,1)
            if xPv == 0:
                return HF(P[0].expansion(0), 0,1)
        elif yPv == 0:
            if xPv > 0:
                return HF(0, P[1].expansion(0),1)
            if xPv == 0:
                return HF(P[0].expansion(0), P[1].expansion(0),1)
        else:
            return HF(0,1,0)

    def is_same_disc(self, P, Q):
        """
        Check if `P,Q` are in same residue disc.

        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: H = HyperellipticCurve(x^3-10*x+9)
            sage: K = Qp(5,8)
            sage: HK = H.change_ring(K)
            sage: P = HK.lift_x(1 + 2*5^2)
            sage: Q = HK.lift_x(5^-2)
            sage: S = HK(1,0)
            sage: HK.is_same_disc(P,Q)
            False
            sage: HK.is_same_disc(P,S)
            True
            sage: HK.is_same_disc(Q,S)
            False
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
            sage: E.tiny_integrals([1,x],P, TP) == E.tiny_integrals_on_basis(P,TP)
            True

        ::

            sage: K = pAdicField(11, 5)
            sage: x = polygen(K)
            sage: C = HyperellipticCurve(x^5 + 33/16*x^4 + 3/4*x^3 + 3/8*x^2 - 1/4*x + 1/16)
            sage: P = C.lift_x(11^(-2))
            sage: Q = C.lift_x(3*11^(-2))
            sage: C.tiny_integrals([1],P,Q)
            (5*11^3 + 7*11^4 + 2*11^5 + 6*11^6 + 11^7 + O(11^8))

        Note that this fails if the points are not in the same residue disc::

            sage: S = C(0,1/4)
            sage: C.tiny_integrals([1,x,x^2,x^3],P,S)
            Traceback (most recent call last):
            ...
            ValueError: (11^-2 + O(11^3) : 11^-5 + 8*11^-2 + O(11^0) : 1 + O(11^5)) and (0 : 3 + 8*11 + 2*11^2 + 8*11^3 + 2*11^4 + O(11^5) : 1 + O(11^5)) are not in the same residue disc
        """
        x, y, z = self.local_analytic_interpolation(P, Q)  #homogeneous coordinates
        x = x/z
        y = y/z
        dt = x.derivative() / (2*y)
        integrals = []
        g = self.genus()
        for f in F:
            try:
                f_dt = f(x,y)*dt
            except TypeError:   #if f is a constant, not callable
                f_dt = f*dt
            if x.valuation() != -2:
                I = sum(f_dt[n]/(n+1) for n in range(f_dt.degree() + 1)) # \int_0^1 f dt
            else:
                If_dt = f_dt.integral().laurent_polynomial()
                I = If_dt(Q[0]**g/Q[1]) - If_dt(P[0]**g/P[1])
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

        ::

            sage: K = pAdicField(11, 5)
            sage: x = polygen(K)
            sage: C = HyperellipticCurve(x^5 + 33/16*x^4 + 3/4*x^3 + 3/8*x^2 - 1/4*x + 1/16)
            sage: P = C.lift_x(11^(-2))
            sage: Q = C.lift_x(3*11^(-2))
            sage: C.tiny_integrals_on_basis(P,Q)
            (5*11^3 + 7*11^4 + 2*11^5 + 6*11^6 + 11^7 + O(11^8), 10*11 + 2*11^3 + 3*11^4 + 5*11^5 + O(11^6), 5*11^-1 + 8 + 4*11 + 10*11^2 + 7*11^3 + O(11^4), 2*11^-3 + 11^-2 + 11^-1 + 10 + 8*11 + O(11^2))


        Note that this fails if the points are not in the same residue disc::

            sage: S = C(0,1/4)
            sage: C.tiny_integrals_on_basis(P,S)
            Traceback (most recent call last):
            ...
            ValueError: (11^-2 + O(11^3) : 11^-5 + 8*11^-2 + O(11^0) : 1 + O(11^5)) and (0 : 3 + 8*11 + 2*11^2 + 8*11^3 + 2*11^4 + O(11^5) : 1 + O(11^5)) are not in the same residue disc
        """
        if P == Q:
            V = VectorSpace(self.base_ring(), 2*self.genus())
            return V(0)
        R = PolynomialRing(self.base_ring(), ['x', 'y'])
        x, y = R.gens()
        return self.tiny_integrals([x**i for i in range(2*self.genus())], P, Q)

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

            sage: K = pAdicField(11, 5)
            sage: x = polygen(K)
            sage: C = HyperellipticCurve(x^5 + 33/16*x^4 + 3/4*x^3 + 3/8*x^2 - 1/4*x + 1/16)
            sage: P = C.lift_x(2)
            sage: Q = C.lift_x(3)
            sage: C.coleman_integrals_on_basis(P, Q)
            (9*11^2 + 7*11^3 + 5*11^4 + O(11^5), 11 + 3*11^2 + 7*11^3 + 11^4 + O(11^5), 10*11 + 11^2 + 5*11^3 + 5*11^4 + O(11^5), 3 + 9*11^2 + 6*11^3 + 11^4 + O(11^5))
            sage: C.coleman_integrals_on_basis(P, Q, algorithm='teichmuller')
            (9*11^2 + 7*11^3 + 5*11^4 + O(11^5), 11 + 3*11^2 + 7*11^3 + 11^4 + O(11^5), 10*11 + 11^2 + 5*11^3 + 5*11^4 + O(11^5), 3 + 9*11^2 + 6*11^3 + 11^4 + O(11^5))

        ::

            sage: K = pAdicField(11,5)
            sage: x = polygen(K)
            sage: C = HyperellipticCurve(x^5 + 33/16*x^4 + 3/4*x^3 + 3/8*x^2 - 1/4*x + 1/16)
            sage: P = C(11^-2, 10*11^-5 + 10*11^-4 + 10*11^-3 + 2*11^-2 + 10*11^-1)
            sage: Q = C(3*11^-2, 11^-5 + 11^-3 + 10*11^-2 + 7*11^-1)
            sage: C.coleman_integrals_on_basis(P, Q)
            (6*11^3 + 3*11^4 + 8*11^5 + 4*11^6 + 9*11^7 + O(11^8), 11 + 10*11^2 + 8*11^3 + 7*11^4 + 5*11^5 + O(11^6), 6*11^-1 + 2 + 6*11 + 3*11^3 + O(11^4), 9*11^-3 + 9*11^-2 + 9*11^-1 + 2*11 + O(11^2))

        ::

            sage: R = C(0,1/4)
            sage: a = C.coleman_integrals_on_basis(P,R)  # long time (7s on sage.math, 2011)
            sage: b = C.coleman_integrals_on_basis(R,Q)  # long time (9s on sage.math, 2011)
            sage: c = C.coleman_integrals_on_basis(P,Q)  # long time
            sage: a+b == c  # long time
            True

        ::

            sage: R.<x> = QQ['x']
            sage: H = HyperellipticCurve(x^3-10*x+9)
            sage: K = Qp(5,8)
            sage: HK = H.change_ring(K)
            sage: S = HK(1,0)
            sage: P = HK(0,3)
            sage: T = HK(0,1,0)
            sage: Q = HK.lift_x(5^-2)
            sage: R = HK.lift_x(4*5^-2)
            sage: HK.coleman_integrals_on_basis(S,P)
            (2*5^2 + 5^4 + 5^5 + 3*5^6 + 3*5^7 + 2*5^8 + O(5^9), 5 + 2*5^2 + 4*5^3 + 2*5^4 + 3*5^6 + 4*5^7 + 2*5^8 + O(5^9))
            sage: HK.coleman_integrals_on_basis(T,P)
            (2*5^2 + 5^4 + 5^5 + 3*5^6 + 3*5^7 + 2*5^8 + O(5^9), 5 + 2*5^2 + 4*5^3 + 2*5^4 + 3*5^6 + 4*5^7 + 2*5^8 + O(5^9))
            sage: HK.coleman_integrals_on_basis(P,S) == -HK.coleman_integrals_on_basis(S,P)
            True
            sage: HK.coleman_integrals_on_basis(S,Q)
            (5 + O(5^4), 4*5^-1 + 4 + 4*5 + 4*5^2 + O(5^3))
            sage: HK.coleman_integrals_on_basis(Q,R)
            (5 + 2*5^2 + 2*5^3 + 2*5^4 + 3*5^5 + 3*5^6 + 3*5^7 + 5^8 + O(5^9), 3*5^-1 + 2*5^4 + 5^5 + 2*5^6 + O(5^7))
            sage: HK.coleman_integrals_on_basis(S,R) == HK.coleman_integrals_on_basis(S,Q) + HK.coleman_integrals_on_basis(Q,R)
            True
            sage: HK.coleman_integrals_on_basis(T,T)
            (0, 0)
            sage: HK.coleman_integrals_on_basis(S,T)
            (0, 0)

        AUTHORS:

        - Robert Bradshaw (2007-03): non-Weierstrass points
        - Jennifer Balakrishnan and Robert Bradshaw (2010-02): Weierstrass points
        """
        import sage.schemes.hyperelliptic_curves.monsky_washnitzer as monsky_washnitzer
        from sage.misc.profiler import Profiler
        prof = Profiler()
        prof("setup")
        K = self.base_ring()
        p = K.prime()
        prec = K.precision_cap()
        g = self.genus()
        dim = 2*g
        V = VectorSpace(K, dim)
        #if P or Q is Weierstrass, use the Frobenius algorithm
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
        elif self.is_same_disc(P,Q):
            return self.tiny_integrals_on_basis(P,Q)
        elif algorithm == 'teichmuller':
            prof("teichmuller")
            PP = TP = self.teichmuller(P)
            QQ = TQ = self.teichmuller(Q)
        else:
            prof("frobPQ")
            TP = self.frobenius(P)
            TQ = self.frobenius(Q)
            PP, QQ = P, Q
        prof("tiny integrals")
        if TP is None:
            P_to_TP = V(0)
        else:
            if TP is not None:
                TPv = (TP[0]**g/TP[1]).valuation()
                xTPv = TP[0].valuation()
            else:
                xTPv = TPv = +Infinity
            if TQ is not None:
                TQv = (TQ[0]**g/TQ[1]).valuation()
                xTQv = TQ[0].valuation()
            else:
                xTQv = TQv = +Infinity
            offset = (2*g-1)*max(TPv, TQv)
            if offset == +Infinity:
                offset = (2*g-1)*min(TPv,TQv)
            if (offset > prec and (xTPv < 0 or xTQv < 0) and (self.residue_disc(P) == self.change_ring(GF(p))(0,1,0) or self.residue_disc(Q) == self.change_ring(GF(p))(0,1,0))):
                newprec = offset + prec
                K = pAdicField(p,newprec)
                A = PolynomialRing(RationalField(),'x')
                f = A(self.hyperelliptic_polynomials()[0])
                from sage.schemes.hyperelliptic_curves.constructor import HyperellipticCurve
                self = HyperellipticCurve(f).change_ring(K)
                xP = P[0]
                xPv = xP.valuation()
                xPnew = K(sum(c * p**(xPv + i) for i, c in enumerate(xP.expansion())))
                PP = P = self.lift_x(xPnew)
                TP = self.frobenius(P)
                xQ = Q[0]
                xQv = xQ.valuation()
                xQnew = K(sum(c * p**(xQv + i) for i, c in enumerate(xQ.expansion())))
                QQ = Q = self.lift_x(xQnew)
                TQ = self.frobenius(Q)
                V = VectorSpace(K,dim)
            P_to_TP = V(self.tiny_integrals_on_basis(P, TP))
        if TQ is None:
            TQ_to_Q = V(0)
        else:
            TQ_to_Q = V(self.tiny_integrals_on_basis(TQ, Q))
        prof("mw calc")
        try:
            M_frob, forms = self._frob_calc
        except AttributeError:
            M_frob, forms = self._frob_calc = monsky_washnitzer.matrix_of_frobenius_hyperelliptic(self)
        prof("eval f")
        R = forms[0].base_ring()
        try:
            prof("eval f %s" % R)
            if PP is None:
                L = [-ff(R(QQ[0]), R(QQ[1])) for ff in forms]  ##changed
            elif QQ is None:
                L = [ff(R(PP[0]), R(PP[1])) for ff in forms]
            else:
                L = [ff(R(PP[0]), R(PP[1])) - ff(R(QQ[0]), R(QQ[1]))
                     for ff in forms]
        except ValueError:
            prof("changing rings")
            forms = [ff.change_ring(self.base_ring()) for ff in forms]
            prof("eval f %s" % self.base_ring())
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
        elif algorithm != 'teichmuller':
            b -= P_to_TP + TQ_to_Q
        prof("lin alg")
        M_sys = matrix(K, M_frob).transpose() - 1
        TP_to_TQ = M_sys**(-1) * b
        prof("done")
        if algorithm == 'teichmuller':
            return P_to_TP + TP_to_TQ + TQ_to_Q
        else:
            return TP_to_TQ

    coleman_integrals_on_basis_hyperelliptic = coleman_integrals_on_basis

#    def invariant_differential(self):
#        """
#        Returns the invariant differential `dx/2y` on self
#
#        EXAMPLES::
#
#            sage: R.<x> = QQ['x']
#            sage: H = HyperellipticCurve(x^3+1)
#            sage: K = Qp(5,8)
#            sage: HK = H.change_ring(K)
#            sage: w = HK.invariant_differential(); w
#            (((1+O(5^8)))*1) dx/2y
#
#        ::
#
#            sage: K = pAdicField(11, 6)
#            sage: x = polygen(K)
#            sage: C = HyperellipticCurve(x^5 + 33/16*x^4 + 3/4*x^3 + 3/8*x^2 - 1/4*x + 1/16)
#            sage: C.invariant_differential()
#            (((1+O(11^6)))*1) dx/2y
#
#        """
#        import sage.schemes.hyperelliptic_curves.monsky_washnitzer as monsky_washnitzer
#        S = monsky_washnitzer.SpecialHyperellipticQuotientRing(self)
#        MW = monsky_washnitzer.MonskyWashnitzerDifferentialRing(S)
#        return MW.invariant_differential()

    def coleman_integral(self, w, P, Q, algorithm='None'):
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

        Example of Leprevost from Kiran Kedlaya
        The first two should be zero as `(P-Q) = 30(P-Q)` in the Jacobian
        and `dx/2y` and `x dx/2y` are holomorphic. ::

            sage: K = pAdicField(11, 6)
            sage: x = polygen(K)
            sage: C = HyperellipticCurve(x^5 + 33/16*x^4 + 3/4*x^3 + 3/8*x^2 - 1/4*x + 1/16)
            sage: P = C(-1, 1); P1 = C(-1, -1)
            sage: Q = C(0, 1/4); Q1 = C(0, -1/4)
            sage: x, y = C.monsky_washnitzer_gens()
            sage: w = C.invariant_differential()
            sage: w.coleman_integral(P, Q)
            O(11^6)
            sage: C.coleman_integral(x*w, P, Q)
            O(11^6)
            sage: C.coleman_integral(x^2*w, P, Q)
            7*11 + 6*11^2 + 3*11^3 + 11^4 + 5*11^5 + O(11^6)

        ::

            sage: p = 71; m = 4
            sage: K = pAdicField(p, m)
            sage: x = polygen(K)
            sage: C = HyperellipticCurve(x^5 + 33/16*x^4 + 3/4*x^3 + 3/8*x^2 - 1/4*x + 1/16)
            sage: P = C(-1, 1); P1 = C(-1, -1)
            sage: Q = C(0, 1/4); Q1 = C(0, -1/4)
            sage: x, y = C.monsky_washnitzer_gens()
            sage: w = C.invariant_differential()
            sage: w.integrate(P, Q), (x*w).integrate(P, Q)
            (O(71^4), O(71^4))
            sage: R, R1 = C.lift_x(4, all=True)
            sage: w.integrate(P, R)
            50*71 + 3*71^2 + 43*71^3 + O(71^4)
            sage: w.integrate(P, R) + w.integrate(P1, R1)
            O(71^4)

        A simple example, integrating dx::

            sage: R.<x> = QQ['x']
            sage: E= HyperellipticCurve(x^3-4*x+4)
            sage: K = Qp(5,10)
            sage: EK = E.change_ring(K)
            sage: P = EK(2, 2)
            sage: Q = EK.teichmuller(P)
            sage: x, y = EK.monsky_washnitzer_gens()
            sage: EK.coleman_integral(x.diff(), P, Q)
            5 + 2*5^2 + 5^3 + 3*5^4 + 4*5^5 + 2*5^6 + 3*5^7 + 3*5^9 + O(5^10)
            sage: Q[0] - P[0]
            5 + 2*5^2 + 5^3 + 3*5^4 + 4*5^5 + 2*5^6 + 3*5^7 + 3*5^9 + O(5^10)

        Yet another example::

            sage: R.<x> = QQ['x']
            sage: H = HyperellipticCurve(x*(x-1)*(x+9))
            sage: K = Qp(7,10)
            sage: HK = H.change_ring(K)
            sage: import sage.schemes.hyperelliptic_curves.monsky_washnitzer as mw
            sage: M_frob, forms = mw.matrix_of_frobenius_hyperelliptic(HK)
            sage: w = HK.invariant_differential()
            sage: x,y = HK.monsky_washnitzer_gens()
            sage: f = forms[0]
            sage: S = HK(9,36)
            sage: Q = HK.teichmuller(S)
            sage: P = HK(-1,4)
            sage: b = x*w*w._coeff.parent()(f)
            sage: HK.coleman_integral(b,P,Q)
            7 + 7^2 + 4*7^3 + 5*7^4 + 3*7^5 + 7^6 + 5*7^7 + 3*7^8 + 4*7^9 + 4*7^10 + O(7^11)

        ::

            sage: R.<x> = QQ['x']
            sage: H = HyperellipticCurve(x^3+1)
            sage: K = Qp(5,8)
            sage: HK = H.change_ring(K)
            sage: w = HK.invariant_differential()
            sage: P = HK(0,1)
            sage: Q = HK(5, 1 + 3*5^3 + 2*5^4 + 2*5^5 + 3*5^7)
            sage: x,y = HK.monsky_washnitzer_gens()
            sage: (2*y*w).coleman_integral(P,Q)
            5 + O(5^9)
            sage: xloc,yloc,zloc = HK.local_analytic_interpolation(P,Q)
            sage: I2 = (xloc.derivative()/(2*yloc)).integral()
            sage: I2.polynomial()(1) - I2(0)
            3*5 + 2*5^2 + 2*5^3 + 5^4 + 4*5^6 + 5^7 + O(5^9)
            sage: HK.coleman_integral(w,P,Q)
            3*5 + 2*5^2 + 2*5^3 + 5^4 + 4*5^6 + 5^7 + O(5^9)

        Integrals involving Weierstrass points::

            sage: R.<x> = QQ['x']
            sage: H = HyperellipticCurve(x^3-10*x+9)
            sage: K = Qp(5,8)
            sage: HK = H.change_ring(K)
            sage: S = HK(1,0)
            sage: P = HK(0,3)
            sage: negP = HK(0,-3)
            sage: T = HK(0,1,0)
            sage: w = HK.invariant_differential()
            sage: x,y = HK.monsky_washnitzer_gens()
            sage: HK.coleman_integral(w*x^3,S,T)
            0
            sage: HK.coleman_integral(w*x^3,T,S)
            0
            sage: HK.coleman_integral(w,S,P)
            2*5^2 + 5^4 + 5^5 + 3*5^6 + 3*5^7 + 2*5^8 + O(5^9)
            sage: HK.coleman_integral(w,T,P)
            2*5^2 + 5^4 + 5^5 + 3*5^6 + 3*5^7 + 2*5^8 + O(5^9)
            sage: HK.coleman_integral(w*x^3,T,P)
            5^2 + 2*5^3 + 3*5^6 + 3*5^7 + O(5^8)
            sage: HK.coleman_integral(w*x^3,S,P)
            5^2 + 2*5^3 + 3*5^6 + 3*5^7 + O(5^8)
            sage: HK.coleman_integral(w, P, negP, algorithm='teichmuller')
            5^2 + 4*5^3 + 2*5^4 + 2*5^5 + 3*5^6 + 2*5^7 + 4*5^8 + O(5^9)
            sage: HK.coleman_integral(w, P, negP)
            5^2 + 4*5^3 + 2*5^4 + 2*5^5 + 3*5^6 + 2*5^7 + 4*5^8 + O(5^9)

        AUTHORS:

        - Robert Bradshaw (2007-03)
        - Kiran Kedlaya (2008-05)
        - Jennifer Balakrishnan (2010-02)
        """
        # TODO: implement Jacobians and show the relationship directly
        import sage.schemes.hyperelliptic_curves.monsky_washnitzer as monsky_washnitzer
        K = self.base_ring()
        prec = K.precision_cap()
        S = monsky_washnitzer.SpecialHyperellipticQuotientRing(self, K)
        MW = monsky_washnitzer.MonskyWashnitzerDifferentialRing(S)
        w = MW(w)
        f, vec = w.reduce_fast()
        basis_values = self.coleman_integrals_on_basis(P, Q, algorithm)
        dim = len(basis_values)
        x,y = self.local_coordinates_at_infinity(2*prec)
        if self.is_weierstrass(P):
            if self.is_weierstrass(Q):
                return 0
            elif f == 0:
                return sum([vec[i] * basis_values[i] for i in range(dim)])
            elif w._coeff(x,-y)*x.derivative()/(-2*y)+w._coeff(x,y)*x.derivative()/(2*y) == 0:
                return self.coleman_integral(w,self(Q[0],-Q[1]), self(Q[0],Q[1]), algorithm)/2
            else:
                raise ValueError("The differential is not odd: use coleman_integral_from_weierstrass_via_boundary")

        elif self.is_weierstrass(Q):
            if f == 0:
                return sum([vec[i] * basis_values[i] for i in range(dim)])
            elif w._coeff(x,-y)*x.derivative()/(-2*y)+w._coeff(x,y)*x.derivative()/(2*y) == 0:
                return -self.coleman_integral(w,self(P[0],-P[1]), self(P[0],P[1]), algorithm)/2
            else:
                raise ValueError("The differential is not odd: use coleman_integral_from_weierstrass_via_boundary")
        else:
            return f(Q[0], Q[1]) - f(P[0], P[1]) + sum([vec[i] * basis_values[i] for i in range(dim)]) # this is just a dot product...

    def frobenius(self, P=None):
        """
        Return the `p`-th power lift of Frobenius of `P`.

        EXAMPLES::

            sage: K = Qp(11, 5)
            sage: R.<x> = K[]
            sage: E = HyperellipticCurve(x^5 - 21*x - 20)
            sage: P = E.lift_x(2)
            sage: E.frobenius(P)
            (2 + 10*11 + 5*11^2 + 11^3 + O(11^5) : 6 + 11 + 8*11^2 + 8*11^3 + 10*11^4 + O(11^5) : 1 + O(11^5))
            sage: Q = E.teichmuller(P); Q
            (2 + 10*11 + 4*11^2 + 9*11^3 + 11^4 + O(11^5) : 6 + 11 + 4*11^2 + 9*11^3 + 4*11^4 + O(11^5) : 1 + O(11^5))
            sage: E.frobenius(Q)
            (2 + 10*11 + 4*11^2 + 9*11^3 + 11^4 + O(11^5) : 6 + 11 + 4*11^2 + 9*11^3 + 4*11^4 + O(11^5) : 1 + O(11^5))

        ::

            sage: R.<x> = QQ[]
            sage: H = HyperellipticCurve(x^5-23*x^3+18*x^2+40*x)
            sage: Q = H(0,0)
            sage: u,v = H.local_coord(Q,prec=100)
            sage: K = Qp(11,5)
            sage: L.<a> = K.extension(x^20-11)
            sage: HL = H.change_ring(L)
            sage: S = HL(u(a),v(a))
            sage: HL.frobenius(S)
            (8*a^22 + 10*a^42 + 4*a^44 + 2*a^46 + 9*a^48 + 8*a^50 + a^52 + 7*a^54 +
            7*a^56 + 5*a^58 + 9*a^62 + 5*a^64 + a^66 + 6*a^68 + a^70 + 6*a^74 +
            2*a^76 + 2*a^78 + 4*a^82 + 5*a^84 + 2*a^86 + 7*a^88 + a^90 + 6*a^92 +
            a^96 + 5*a^98 + 2*a^102 + 2*a^106 + 6*a^108 + 8*a^110 + 3*a^112 +
            a^114 + 8*a^116 + 10*a^118 + 3*a^120 + O(a^122) :
            a^11 + 7*a^33 + 7*a^35 + 4*a^37 + 6*a^39 + 9*a^41 + 8*a^43 + 8*a^45 +
            a^47 + 7*a^51 + 4*a^53 + 5*a^55 + a^57 + 7*a^59 + 5*a^61 + 9*a^63 +
            4*a^65 + 10*a^69 + 3*a^71 + 2*a^73 + 9*a^75 + 10*a^77 + 6*a^79 +
            10*a^81 + 7*a^85 + a^87 + 4*a^89 + 8*a^91 + a^93 + 8*a^95 + 2*a^97 +
            7*a^99 + a^101 + 3*a^103 + 6*a^105 + 7*a^107 + 4*a^109 + O(a^111) :
            1 + O(a^100))

        AUTHORS:

        - Robert Bradshaw and Jennifer Balakrishnan (2010-02)
        """
        try:
            _frob = self._frob
        except AttributeError:
            K = self.base_ring()
            p = K.prime()
            x = K['x'].gen(0)

            f, f2 = self.hyperelliptic_polynomials()
            if f2 != 0:
                raise NotImplementedError("Curve must be in weierstrass normal form.")
            h = (f(x**p) - f**p)

            def _frob(P):
                if P == self(0,1,0):
                    return P
                x0 = P[0]
                y0 = P[1]
                try:
                    uN = (1 + h(x0)/y0**(2*p)).sqrt()
                    yres = y0**p * uN
                    xres = x0**p
                    if (yres-y0).valuation() == 0:
                        yres = -yres
                    return self.point([xres,yres, K(1)])
                except (TypeError, NotImplementedError):
                    uN2 = 1 + h(x0)/y0**(2*p)
                    #yfrob2 = f(x)
                    c = uN2.expansion(0)
                    v = uN2.valuation()
                    a = uN2.parent().gen()
                    uN = self.newton_sqrt(uN2,c.sqrt()*a**(v//2),K.precision_cap())
                    yres = y0**p * uN
                    xres = x0**p
                    if (yres - y0).valuation() == 0:
                        yres = -yres
                    try:
                        return self(xres,yres)
                    except ValueError:
                        return self._curve_over_ram_extn(xres,yres)

            self._frob = _frob

        if P is None:
            return _frob
        else:
            return _frob(P)

    def newton_sqrt(self, f, x0, prec):
        r"""
        Take the square root of the power series `f` by Newton's method.

        NOTE:

        this function should eventually be moved to `p`-adic power series ring

        INPUT:

        - ``f`` -- power series with coefficients in `\QQ_p` or an extension
        - ``x0`` -- seeds the Newton iteration
        - ``prec`` -- precision

        OUTPUT: the square root of `f`

        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: H = HyperellipticCurve(x^5-23*x^3+18*x^2+40*x)
            sage: Q = H(0,0)
            sage: u,v = H.local_coord(Q,prec=100)
            sage: K = Qp(11,5)
            sage: HK = H.change_ring(K)
            sage: L.<a> = K.extension(x^20-11)
            sage: HL = H.change_ring(L)
            sage: S = HL(u(a),v(a))
            sage: f = H.hyperelliptic_polynomials()[0]
            sage: y = HK.newton_sqrt( f(u(a)^11), a^11,5)
            sage: y^2 - f(u(a)^11)
            O(a^122)

        AUTHOR:

        - Jennifer Balakrishnan
        """
        z = x0
        loop_prec = log(RR(prec), 2).ceil()
        for i in range(loop_prec):
            z = (z + f / z) / 2
        return z

    def curve_over_ram_extn(self, deg):
        r"""
        Return ``self`` over `\QQ_p(p^(1/deg))`.

        INPUT:

        - ``deg`` -- the degree of the ramified extension

        OUTPUT: ``self`` over `\QQ_p(p^(1/deg))`

        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: H = HyperellipticCurve(x^5-23*x^3+18*x^2+40*x)
            sage: K = Qp(11,5)
            sage: HK = H.change_ring(K)
            sage: HL = HK.curve_over_ram_extn(2)
            sage: HL
            Hyperelliptic Curve over 11-adic Eisenstein Extension Field in a defined by x^2 - 11 defined by (1 + O(a^10))*y^2 = (1 + O(a^10))*x^5 + (10 + 8*a^2 + 10*a^4 + 10*a^6 + 10*a^8 + O(a^10))*x^3 + (7 + a^2 + O(a^10))*x^2 + (7 + 3*a^2 + O(a^10))*x

        AUTHOR:

        - Jennifer Balakrishnan
        """
        from sage.schemes.hyperelliptic_curves.constructor import HyperellipticCurve
        K = self.base_ring()
        p = K.prime()
        A = PolynomialRing(QQ,'x')
        x = A.gen()
        J = K.extension(x**deg-p,names='a')
        pol = self.hyperelliptic_polynomials()[0]
        H = HyperellipticCurve(A(pol))
        HJ = H.change_ring(J)
        self._curve_over_ram_extn = HJ
        self._curve_over_ram_extn._curve_over_Qp = self
        return HJ

    def get_boundary_point(self, curve_over_extn, P):
        """
        Given ``self`` over an extension field, find a point in the disc of `P`
        near the boundary.

        INPUT:

        - ``curve_over_extn`` -- ``self`` over a totally ramified extension
        - ``P`` -- Weierstrass point

        OUTPUT: a point in the disc of `P` near the boundary

        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: H = HyperellipticCurve(x^3-10*x+9)
            sage: K = Qp(3,6)
            sage: HK = H.change_ring(K)
            sage: P = HK(1,0)
            sage: J.<a> = K.extension(x^30-3)
            sage: HJ  = H.change_ring(J)
            sage: S = HK.get_boundary_point(HJ,P)
            sage: S
            (1 + 2*a^2 + 2*a^6 + 2*a^18 + a^32 + a^34 + a^36 + 2*a^38 + 2*a^40 + a^42 + 2*a^44 + a^48 + 2*a^50 + 2*a^52 + a^54 + a^56 + 2*a^60 + 2*a^62 + a^70 + 2*a^72 + a^76 + 2*a^78 + a^82 + a^88 + a^96 + 2*a^98 + 2*a^102 + a^104 + 2*a^106 + a^108 + 2*a^110 + a^112 + 2*a^116 + a^126 + 2*a^130 + 2*a^132 + a^144 + 2*a^148 + 2*a^150 + a^152 + 2*a^154 + a^162 + a^164 + a^166 + a^168 + a^170 + a^176 + a^178 + O(a^180) : a + O(a^180) : 1 + O(a^180))

        AUTHOR:

        - Jennifer Balakrishnan
        """
        J = curve_over_extn.base_ring()
        a = J.gen()
        prec2 = J.precision_cap()
        x,y = self.local_coord(P,prec2)
        return curve_over_extn(x(a),y(a))

    def P_to_S(self, P, S):
        r"""
        Given a finite Weierstrass point `P` and a point `S` in the same disc,
        compute the Coleman integrals `\{\int_P^S x^i dx/2y \}_{i=0}^{2g-1}`.

        INPUT:

        - ``P`` -- finite Weierstrass point
        - ``S`` -- point in disc of `P`

        OUTPUT: Coleman integrals `\{\int_P^S x^i dx/2y \}_{i=0}^{2g-1}`

        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: H = HyperellipticCurve(x^3-10*x+9)
            sage: K = Qp(5,4)
            sage: HK = H.change_ring(K)
            sage: P = HK(1,0)
            sage: HJ = HK.curve_over_ram_extn(10)
            sage: S = HK.get_boundary_point(HJ,P)
            sage: HK.P_to_S(P, S)
            (2*a + 4*a^3 + 2*a^11 + 4*a^13 + 2*a^17 + 2*a^19 + a^21 + 4*a^23 + a^25 + 2*a^27 + 2*a^29 + 3*a^31 + 4*a^33 + O(a^35), a^-5 + 2*a + 2*a^3 + a^7 + 3*a^11 + a^13 + 3*a^15 + 3*a^17 + 2*a^19 + 4*a^21 + 4*a^23 + 4*a^25 + 2*a^27 + a^29 + a^31 + O(a^33))

        AUTHOR:

        - Jennifer Balakrishnan
        """
        prec = self.base_ring().precision_cap()
        deg = (S[0]).parent().defining_polynomial().degree()
        prec2 = prec*deg
        x,y = self.local_coord(P,prec2)
        g = self.genus()
        integrals = [((x**k*x.derivative()/(2*y)).integral()) for k in range(2*g)]
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

            sage: R.<x> = QQ['x']
            sage: H = HyperellipticCurve(x^3-10*x+9)
            sage: K = Qp(5,4)
            sage: HK = H.change_ring(K)
            sage: P = HK(1,0)
            sage: J.<a> = K.extension(x^10-5)
            sage: HJ  = H.change_ring(J)
            sage: S = HK.get_boundary_point(HJ,P)
            sage: x,y = HK.monsky_washnitzer_gens()
            sage: S[0]-P[0] == HK.coleman_integral_P_to_S(x.diff(),P,S)
            True
            sage: HK.coleman_integral_P_to_S(HK.invariant_differential(),P,S) == HK.P_to_S(P,S)[0]
            True

        AUTHOR:

        - Jennifer Balakrishnan
        """
        prec = self.base_ring().precision_cap()
        deg = S[0].parent().defining_polynomial().degree()
        prec2 = prec * deg
        x,y = self.local_coord(P,prec2)
        int_sing = (w.coeff()(x,y)*x.derivative()/(2*y)).integral()
        int_sing_a = int_sing(S[1])
        return int_sing_a

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

            sage: R.<x> = QQ['x']
            sage: H = HyperellipticCurve(x^3-10*x+9)
            sage: K = Qp(5,6)
            sage: HK = H.change_ring(K)
            sage: J.<a> = K.extension(x^20-5)
            sage: HJ  = H.change_ring(J)
            sage: w = HK.invariant_differential()
            sage: x,y = HK.monsky_washnitzer_gens()
            sage: P = HK(1,0)
            sage: Q = HK(0,3)
            sage: S = HK.get_boundary_point(HJ,P)
            sage: P_to_S = HK.P_to_S(P,S)
            sage: S_to_Q = HJ.S_to_Q(S,Q)
            sage: P_to_S + S_to_Q
            (2*a^40 + a^80 + a^100 + O(a^105), a^20 + 2*a^40 + 4*a^60 + 2*a^80 + O(a^103))
            sage: HK.coleman_integrals_on_basis(P,Q)
            (2*5^2 + 5^4 + 5^5 + 3*5^6 + O(5^7), 5 + 2*5^2 + 4*5^3 + 2*5^4 + 5^6 + O(5^7))

        AUTHOR:

        - Jennifer Balakrishnan
        """
        FS = self.frobenius(S)
        FS = (FS[0],FS[1])
        FQ = self.frobenius(Q)
        import sage.schemes.hyperelliptic_curves.monsky_washnitzer as monsky_washnitzer
        try:
            M_frob, forms = self._frob_calc
        except AttributeError:
            M_frob, forms = self._frob_calc = monsky_washnitzer.matrix_of_frobenius_hyperelliptic(self)
        try:
            HJ = self._curve_over_ram_extn
            K = HJ.base_ring()
        except AttributeError:
            HJ = S.scheme()
            K = self.base_ring()
        g = self.genus()
        prec2 = K.precision_cap()
        p = K.prime()
        dim = 2*g
        V = VectorSpace(K,dim)
        if S == FS:
            S_to_FS = V(dim*[0])
        else:
            P = self(ZZ(FS[0].expansion(0)),ZZ(FS[1].expansion(0)))
            x,y = self.local_coord(P,prec2)
            integrals = [(x**i*x.derivative()/(2*y)).integral() for i in range(dim)]
            S_to_FS = vector([I.polynomial()(FS[1]) - I.polynomial()(S[1]) for I in integrals])
        if HJ(Q[0],Q[1]) == HJ(FQ):
            FQ_to_Q = V(dim*[0])
        else:
            FQ_to_Q = V(self.tiny_integrals_on_basis(FQ, Q))
        try:
            L = [f(K(S[0]), K(S[1])) - f(K(Q[0]), K(Q[1])) for f in forms]
        except ValueError:
            forms = [f.change_ring(K) for f in forms]
            L = [f(S[0], S[1]) - f(Q[0], Q[1]) for f in forms]
        b = V(L)
        M_sys = matrix(K, M_frob).transpose() - 1
        B = (~M_sys)
        vv = min(c.valuation() for c in B.list())
        B = (p**(-vv)*B).change_ring(K)
        B = p**(vv)*B
        return B*(b-S_to_FS-FQ_to_Q)

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

            sage: R.<x> = QQ['x']
            sage: H = HyperellipticCurve(x^3-10*x+9)
            sage: K = Qp(5,6)
            sage: HK = H.change_ring(K)
            sage: J.<a> = K.extension(x^20-5)
            sage: HJ  = H.change_ring(J)
            sage: x,y = HK.monsky_washnitzer_gens()
            sage: P = HK(1,0)
            sage: Q = HK(0,3)
            sage: S = HK.get_boundary_point(HJ,P)
            sage: P_to_S = HK.coleman_integral_P_to_S(y.diff(),P,S)
            sage: S_to_Q = HJ.coleman_integral_S_to_Q(y.diff(),S,Q)
            sage: P_to_S  + S_to_Q
            3 + O(a^119)
            sage: HK.coleman_integral(y.diff(),P,Q)
            3 + O(5^6)

        AUTHOR:

        - Jennifer Balakrishnan
        """
        import sage.schemes.hyperelliptic_curves.monsky_washnitzer as monsky_washnitzer
        K = self.base_ring()
        R = monsky_washnitzer.SpecialHyperellipticQuotientRing(self, K)
        MW = monsky_washnitzer.MonskyWashnitzerDifferentialRing(R)
        w = MW(w)
        f, vec = w.reduce_fast()
        g = self.genus()
        const = f(Q[0],Q[1])-f(S[0],S[1])
        if vec == vector(2*g*[0]):
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

            sage: R.<x> = QQ['x']
            sage: H = HyperellipticCurve(x^3-10*x+9)
            sage: K = Qp(5,6)
            sage: HK = H.change_ring(K)
            sage: P = HK(1,0)
            sage: Q = HK(0,3)
            sage: x,y = HK.monsky_washnitzer_gens()
            sage: HK.coleman_integral_from_weierstrass_via_boundary(y.diff(),P,Q,20)
            3 + O(a^119)
            sage: HK.coleman_integral(y.diff(),P,Q)
            3 + O(5^6)
            sage: w = HK.invariant_differential()
            sage: HK.coleman_integral_from_weierstrass_via_boundary(w,P,Q,20)
            2*a^40 + a^80 + a^100 + O(a^105)
            sage: HK.coleman_integral(w,P,Q)
            2*5^2 + 5^4 + 5^5 + 3*5^6 + O(5^7)

        AUTHOR:

        - Jennifer Balakrishnan
        """
        HJ = self.curve_over_ram_extn(d)
        S = self.get_boundary_point(HJ, P)
        P_to_S = self.coleman_integral_P_to_S(w, P, S)
        S_to_Q = HJ.coleman_integral_S_to_Q(w, S, Q)
        return P_to_S + S_to_Q
