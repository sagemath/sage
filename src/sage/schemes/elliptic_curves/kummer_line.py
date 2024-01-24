r"""
Kummer Line for short Weierstrass elliptic curves

Sage defines Kummer lines derived from an elliptic curve in short Weierstrass or Montgomery form, known as `x`-only arithmetic.
Also defines a point on a Kummer line as an element of the projective line as well as isogenies between Kummer lines.

EXAMPLES:

We construct a Kummer line over an already defined elliptic curve in short Weierstrass form::

    sage: E = EllipticCurve(GF(101), [2, 3])
    sage: KummerLine(E)
    Kummer line of the elliptic curve y^2 = x^3 + 2*x + 3 over Finite Field of size 101

We can also give the curve coefficients directly with a base ring::

    sage: KummerLine(QQ, [4, 5/6])
    Kummer line of the elliptic curve y^2 = x^3 + 4*x + 5/6 over Rational Field

TODO KummerPoint examples, as well as many more usage (see integer for instance)
TODO <Lots and lots of examples>

AUTHORS:

- Giacomo Pope (2023): original code

- Elif Özbay Gürler, Nicolas Sarkis (2024-01): ported code to short Weierstrass and documentation
"""

# ****************************************************************************
#       Copyright (C) 2023 Giacomo Pope <giacomopope@gmail.com>
#       Copyright (C) 2024 Elif Özbay Gürler <TODO>
#       Copyright (C) 2024 Nicolas Sarkis <nicolas.sarkis@math.u-bordeaux.fr>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

# from sage.all import cached_method, EllipticCurve
from sage.all import cached_method


from sage.schemes.elliptic_curves.ell_generic import EllipticCurve_generic


class KummerLine:
    r"""
    Kummer line class for a short Weierstrass elliptic curve.

    EXAMPLES::
        sage: E = EllipticCurve(GF(101), [2, 3])
        sage: KummerLine(E)
        Kummer line of the elliptic curve y^2 = x^3 + 2*x + 3 over Finite Field of size 101
    """

    def __init__(self, curve):
        r"""
        Constructor for a Kummer line from a short Weierstrass elliptic curve.

        INPUT::

            - ``curve`` -- an elliptic curve in short Weierstrass form `y^2 = x^3 + a*x + b`

        EXAMPLES::

            sage: E = EllipticCurve(GF(101), [2, 3])
            sage: K = KummerLine(E); K
            Kummer line of the elliptic curve y^2 = x^3 + 2*x + 3 over Finite Field of size 101
        """

        if not isinstance(curve, EllipticCurve_generic):
            raise TypeError(
                "input must be an elliptic curve in short Weierstrass form."
            )

        ainvs = curve.a_invariants()
        a, b = ainvs[3], ainvs[4]  # a = a4, b = a6
        if ainvs != (0, 0, 0, a, b):
            raise ValueError(
                "input must be an elliptic curve in short Weierstrass form."
            )
        # TODO could deal with Montgomery model using an isomorphism

        self._curve = curve
        self._base_ring = curve.base_ring()

        # Initialize variables
        self._a = self._base_ring(a)
        self._b = self._base_ring(b)

        # TODO
        # Not necessary since this check is already done while initializing the curve
        # # Make sure the curve is not singular
        # if self.discriminant() == 0:
        #     raise ValueError(
        #         f"Constants {curve_constants} do not define an elliptic curve in short \
        #         Weierstrass form."
        #     )

    def __eq__(self, other):
        r"""
        Test equality of Kummer lines by checking if the above curves are the same.

        EXAMPLES::

        Base ring must be the same::

            sage: E1 = EllipticCurve(GF(101), [2, 3])
            sage: E2 = EllipticCurve(QQ, [2, 3])
            sage: K1 = KummerLine(E1)
            sage: K2 = KummerLine(E2)
            sage: K1 == K2
            False

        Kummer lines are not equals even if the curves are isomorphic, which
        is also the case for elliptic curves::

            sage: E1 = EllipticCurve(GF(101), [2, 3])
            sage: E2 = EllipticCurve(GF(101), [2, 98])
            sage: K1 = KummerLine(E1)
            sage: K2 = KummerLine(E2)
            sage: K1 == K2
            False
            sage: E1.is_isomorphic(E2)
            True
        """
        return self.curve() == other.curve()

    def __repr__(self):
        r"""
        String representation of a Kummer line.

        EXAMPLES::

            sage: E = EllipticCurve(GF(101), [2, 3])
            sage: K = KummerLine(E); K.__repr__()
            'Kummer line of the elliptic curve y^2 = x^3 + 2*x + 3 over Finite Field of size 101'
        """
        return f"Kummer line of the elliptic curve y^2 = x^3 + {self._a}*x + {self._b} over {self.base_ring()}"

    # TODO
    # def __call__(self, coords):
    #     r"""
    #     Create a point from the coordinates.
    #
    #     INPUT::
    #
    #         - ``coords`` - either a point P on EllipticCurve or a list or tuple (X, Z) where P = (X : * : Z); Z is optional
    #
    #     EXAMPLES::
    #
    #         sage: E = EllipticCurve(GF(101), [2, 3])
    #         sage: P = E(95, 49)
    #         sage: K = KummerLine(E)
    #
    #     Point as a parameter::
    #         sage: K(P)
    #         Kummer Point [95 : 1] on Kummer line of the elliptic curve y^2 = x^3 + 2*x + 3 over Finite Field of size 101
    #
    #     Same output if we just give the XZ-coordinates::
    #         sage: K([95, 1]) == K(P)
    #         True
    #
    #     Z is optional::
    #         sage: K(95) == K(P)
    #         True
    #
    #         sage: K([95]) == K(P)
    #         True
    #     """
    #     return KummerPoint(self, coords)

    def base_ring(self):
        r"""
        Return the base ring of the Kummer Line.

        EXAMPLES::

            sage: E = EllipticCurve(GF(101), [2, 3])
            sage: K = KummerLine(E); K.base_ring()
            Finite Field of size 101
        """
        return self._base_ring

    def extract_constants(self):
        r"""
        Return the short Weierstrass coefficients a, b as a tuple.

        EXAMPLES::

            sage: E = EllipticCurve(GF(101), [2, 3])
            sage: K = KummerLine(E); K.extract_constants()
            (2, 3)
        """
        return self._a, self._b

    def zero(self):
        r"""
        Return the identity point on the Kummer Line.

        EXAMPLES::

            sage: E = EllipticCurve(GF(101), [2, 3])
            sage: K = KummerLine(E); K.zero()
            (1 : 0)
        """
        return self(None)

    def curve(self):
        r"""
        Returnn the associated elliptic curve.

        EXAMPLES::

            sage: E = EllipticCurve(GF(101), [2, 3])
            sage: K = KummerLine(E); K.curve()
            Elliptic Curve defined by y^2 = x^3 + 2*x + 3 over Finite Field of size 101
        """
        return self._curve

    # TODO will be useful for Montgomery but not at the moment
    # @cached_method
    # def short_weierstrass_curve(self):
    #     r"""
    #     Return the short Weierstrass curve associated with the Kummer line.
    #
    #     EXAMPLES::
    #
    #         sage: E = EllipticCurve(GF(101), [2, 3])
    #         sage: K = KummerLine(E)
    #         sage: K.short_weierstrass_curve() == E
    #         True
    #
    #         sage: K = KummerLine(QQ, [3, 4/5])
    #         sage: K.short_weierstreass_curve()
    #         Elliptic Curve defined by y^2 = x^3 + 3*x + 4/5 over Rational Field
    #     """
    #
    #     F = self.base_ring()
    #     a, b = self.extract_constants()
    #     return EllipticCurve(F, [a, b])

    @cached_method
    def j_invariant(self):
        r"""
        Compute the j-invariant of the Kummer line using the formula

        ..MATH::

            j(E) = 1728\frac{4a^3}{4a^3 + 27b^2}.

        EXAMPLES::

            sage: E = EllipticCurve(j=42)
            sage: K = KummerLine(E)
            sage: K.j_invariant() == 42
            True

        ::

            sage: K = KummerLine(GF(101), [2, 3])
            sage: K.j_invariant()
            74
        """
        a, b = self.extract_constants()

        j_num = 4 * a**3
        j_den = j_num + 27 * b**2
        j_num = 1728 * j_num
        return j_num / j_den

    @cached_method
    def discriminant(self):
        r"""
        Compute the discriminant of the Kummer line using the formula

        ..MATH::

            \Delta(E) = -16(4a^3 + 27b^2).

        EXAMPLES::

            sage: K = KummerLine(GF(101), [2, 3])
            sage: K.discriminant()
            44
        """
        a, b = self.extract_constants()
        return -16 * (4 * a**3 + 27 * b**2)

    # TODO code cleanup (new class) + fixing 2P // P.double()
    # TODO isogeny call instantiate a KummerLineIsogeny object
    # def isogeny(self, S, P):
    #     # TODO
    #     XP, _ZP = P.XZ()
    #     a4, a6 = self.extract_constants()
    #     if self.zero() in S:
    #         S.remove(self.zero())
    #     v, w = 0, 0
    #     alpha = XP
    #     for Q in S:
    #         XQ, _ZQ = Q.XZ()
    #         gQx = 3 * XQ**2 + a4
    #         if Q.double() == self.zero():
    #             # if 2 * Q == self.zero(): # broken
    #             # # print((2 * Q).x(), (Q.double().x()))
    #             # print(type(2 * Q), 2 * Q)
    #             # DEBUG
    #             # print(Q)
    #             # print(Q.double(), Q.double().is_zero())
    #             # print(2 * Q, (2 * Q).is_zero())
    #             vQ = gQx
    #         else:
    #             vQ = 2 * gQx
    #         uQ = 4 * (XQ**3 + a4 * XQ + a6)
    #         v += vQ
    #         w += uQ + XQ * vQ
    #         alpha += vQ / (XP - XQ) + uQ / (XP - XQ) ** 2
    #     A4 = a4 - 5 * v
    #     A6 = a6 - 7 * w
    #     K2 = KummerLine(self.base_ring(), [A4, A6])
    #     return K2, K2(alpha)


# class KummerPoint:
#     r"""
#     Kummer line point class
#
#     EXAMPLES::
#
#         sage: E = EllipticCurve(GF(101), [2, 3])
#         sage: K = KummerLine(E)
#         sage: P = E(95, 52)
#
#     A KummerPoint can be constructed from coordinates::
#         sage: xP = KummerPoint(K, [95, 1]); xP
#         Kummer Point [95 : 1] on Kummer line of the elliptic curve y^2 = x^3 + 2*x + 3 over Finite Field of size 101
#
#     Z-coordinate is optional, assume it is 1 then::
#         sage: xP = KummerPoint(K, 95); xP
#         Kummer Point [95 : 1] on Kummer line of the elliptic curve y^2 = x^3 + 2*x + 3 over Finite Field of size 101
#
#     Or, it can be constructed from a point on the elliptic curve::
#         sage: xP = KummerPoint(K, P); xP
#         Kummer Point [95 : 1] on Kummer line of the elliptic curve y^2 = x^3 + 2*x + 3 over Finite Field of size 101
#
#     The point can also be created by calling the Kummer line::
#         sage: xP = K(P); xP
#         Kummer Point [95 : 1] on Kummer line of the elliptic curve y^2 = x^3 + 2*x + 3 over Finite Field of size 101
#     """
#
#     def __init__(self, parent, coords):
#         r"""
#         Create a point from the coordinates.
#
#         INPUT::
#
#             - ``parent`` - A Kummer line
#
#             - ``coords`` - either a point P on EllipticCurve or a list or tuple (X, Z) where P = (X : * : Z); Z is optional
#
#         EXAMPLES::
#
#             sage: E = EllipticCurve(GF(101), [2, 3])
#             sage: P = E(95, 49)
#             sage: K = KummerLine(E)
#
#         Point as a parameter::
#             sage: xP = KummerPoint(K, [95, 1]); xP
#             Kummer Point [95 : 1] on Kummer line of the elliptic curve y^2 = x^3 + 2*x + 3 over Finite Field of size 101
#
#         Same output if we just give the XZ-coordinates::
#             sage: KummerPoint(K, [95, 1]) == KummerPoint(K, P)
#             True
#
#         Z is optional::
#             sage: KummerPoint(K, 95) == KummerPoint(K, P)
#             True
#
#             sage: KummerPoint(K, [95]) == KummerPoint(K, P)
#             True
#         """
#         # Ensure the parent is the right type
#         if not isinstance(parent, KummerLine):
#             raise TypeError("not a Weierstrass Kummer line")
#
#         R = parent.base_ring()
#
#         # Point at infinity
#         if coords is None:
#             coords = (Integer(1), Integer(0))
#         # Construct point from P on an elliptic curve in Montgomery form
#         elif isinstance(coords, EllipticCurvePoint_field):
#             # Make sure point's parent curve matches with Kummer Line
#             a, b = parent._a, parent._b
#             assert coords.curve().a_invariants() == (0, 0, 0, a, b)
#             coords = coords[0], coords[2]
#         # Construct from X coordinate only
#         elif isinstance(coords, RingElement):
#             coords = (coords,)
#         # Construct from a tuple (X : Z)
#         else:
#             coords = tuple(coords)
#
#         # Sanitise the input coordinates
#         if len(coords) == 1:
#             coords += (R.one(),)
#         if len(coords) != 2:
#             raise ValueError("not a point on ℙ¹")
#         coords = tuple(map(R, coords))
#
#         self._base_ring = R
#         self._parent = parent
#         self._X, self._Z = coords
#
#         # TODO does not consider a potential twist of the curve
#         # TODO broken
#         # if not self.is_x_coord():
#         #     raise ValueError("Not a valid x-coordinate")
#
#     def __repr__(self):
#         r"""
#         String representation of a point on a Kummer line.
#
#         EXAMPLES::
#
#             sage: E = EllipticCurve(GF(101), [2, 3])
#             sage: P = E(95, 49)
#             sage: K = KummerLine(E)
#             sage: xP = K(P); xP.__repr__()
#             TODO
#         """
#         return f"Kummer Point [{self._X} : {self._Z}] on {self.parent()}"
#
#     def __bool__(self):
#         r"""
#         Boolean value for a Kummer line point.
#
#         It represents False if it is the point at infinity and True otherwise
#
#         EXAMPLES::
#
#             sage: E = EllipticCurve(GF(101), [2, 3])
#             sage: P = E(95, 49)
#             sage: K = KummerLine(E)
#             sage: xO = K(None)
#             sage: xO == False
#             True
#
#             sage: xP = K(P)
#             sage: xP == True
#             True
#         """
#
#         return bool(self._Z)
#
#     def __eq__(self, other):
#         r"""
#         Equality of two Kummer points.
#
#         They are equal if the above line match and if the projective points are the same.
#
#         INPUT::
#
#             - ``other`` - A Kummer point
#
#         EXAMPLES::
#
#             sage: E = EllipticCurve(GF(101), [2, 3])
#             sage: K = KummerLine(E)
#
#         Two points can have different coordinates and still be the same projectively::
#             sage: P = E(95, 49)
#             sage: xP1 = K(P)
#             sage: xP2 = K([29, 12])
#             sage: xP1 == xP2  # 95*12 % 101 = 29
#             True
#
#         They must represent the same x-coordinate::
#             sage: P = E(95, 49)
#             sage: Q = E(, ) TODO
#             sage: xP, xQ = K(P), K(Q)
#             sage: xP == xQ
#             False
#
#         Base Kummer line must be the same even if the x-coordinates are equal::
#             sage: E2 = EllipticCurve(GF(101), [3, 4])
#             sage: P2 = E2(95, ) TODO
#             sage: K2 = KummerLine(E2)
#             sage: xP1 = K(P)
#             sage: xP2 = K2(P2)
#             sage: xP == xP2
#             False
#         """
#
#         if not isinstance(other, KummerPoint):
#             raise ValueError("Can only compare equality between to Kummer Points")
#         if self._parent != other._parent:
#             return False
#         return self._X * other._Z == other._X * self._Z
#
#     def is_zero(self):
#         r"""
#         A Kummer Point is zero if it corresponds to the point at infinity
#         on the parent curve.
#
#         EXAMPLES::
#
#             sage: E = EllipticCurve(GF(101), [2, 3])
#             sage: P = E(95, 49)
#             sage: K = KummerLine(E)
#             sage: xP = K(P); xP.is_zero()
#
#             sage: xO = K(None); xO.is_zero()
#             True
#         """
#
#         return self._Z == 0
#
#     def is_x_coord(self):
#         r"""
#         Check if the x coordinate of a Kummer Point is a valid
#         one on the parent curve.
#
#         False should never be encountered as an error is raised on initialisation if so.
#
#         EXAMPLES::
#
#             sage: E = EllipticCurve(GF(101), [2, 3])
#             sage: P = E(95, 49)
#             sage: K = KummerLine(E)
#             sage: xP = K(95); xP.is_x_coord()
#             True
#         """
#
#         if self.is_zero():
#             return True
#         return self._parent.curve().is_x_coord(self.x())
#
#     def base_ring(self):
#         r"""
#         Return the base ring of the Kummer Point.
#
#         EXAMPLES::
#
#             sage: E = EllipticCurve(GF(101), [2, 3])
#             sage: P = E(95, 49)
#             sage: K = KummerLine(E)
#             sage: xP = K(P)
#             sage: xP.base_ring()
#             Finite Field of size 101
#         """
#
#         return self._base_ring
#
#     def parent(self):
#         r"""
#         Return the Kummer line on which the point is constructed.
#
#         EXAMPLES::
#
#             sage: E = EllipticCurve(GF(101), [2, 3])
#             sage: P = E(95, 49)
#             sage: K = KummerLine(E)
#             sage: xP = K(P)
#             sage: xP.parent()
#             TODO
#         """
#
#         return self._parent
#
#     def XZ(self):
#         r"""
#         Return the projective coordinates (X : Z) of the point.
#
#         EXAMPLES::
#
#             sage: E = EllipticCurve(GF(101), [2, 3])
#             sage: P = E(95, 49)
#             sage: K = KummerLine(E)
#             sage: xP = K(P)
#             sage: xP.XZ()
#             (95, 1)
#         """
#
#         return self._X, self._Z
#
#     def x(self):
#         r"""
#         Return the affine x-coordinate of the point, if well-defined.
#
#         Raise an error on the point at infinity.
#
#         EXAMPLES::
#
#             sage: E = EllipticCurve(GF(101), [2, 3])
#             sage: P = E(95, 49)
#             sage: K = KummerLine(E)
#             sage: xP = K(P)
#             sage: xP.x()
#             95
#
#         The point at infinity has no valid x-coordinate::
#             sage: xO = K(None)
#             sage: xO.x()
#             TODO
#         """
#
#         if self.is_zero():
#             raise ValueError("The identity point has no valid x-coordinate")
#         if self._Z == 1:
#             return self._base_ring(self._X)
#         return self._base_ring(self._X / self._Z)
#
#     @cached_method
#     def curve_point(self):
#         r"""
#         Return a lift of the point on the parent curve.
#         y is chose as the smallest square root.
#
#         EXAMPLES::
#
#             sage: E = EllipticCurve(GF(101), [2, 3])
#             sage: P = E(95, 49)
#             sage: K = KummerLine(E)
#             sage: xP = K(P)
#             sage: xP.curve_point()
#             TODO
#         """
#         # Get the Montgomery curve and constant A
#         L = self.parent()
#         E = L.curve()
#         a, b = L._a, L._b
#
#         # TODO change to lift_x
#         # Compute y2, assume x is a valid coordinate
#         x = self.x()
#         y2 = x**3 + a * x + b
#         y = y2.sqrt()
#         return E(x, y)
#
#     @staticmethod
#     def xDBL(X, Z, a, b2, b4):
#         r"""
#         Doubling formula for x-only coordinate on short Weierstrass curve.
#
#         Cost: 2M + 5S + 3m0 + 8a
#
#         INPUT::
#
#             - ``X``, ``Z`` -- the coordinates of the Kummer point `P = (X : Z)`
#
#             - ``a``, ``b2``, ``b4`` -- curve constants where `b2 = 2*b` and `b4 = 4*b`
#
#         OUTPUT:: The coordinates of `[2]P = (X2 : Z2)`
#
#         EXAMPLES::
#
#             TODO
#         """
#
#         XX = X * X
#         ZZ = Z * Z
#         t0 = X + Z
#         t1 = t0 * t0 - XX - ZZ
#         t1 = t1 + t1
#         t0 = a * ZZ
#         X2 = XX - t0
#         X2 = X2 * X2 - b2 * t1 * ZZ
#         Z2 = XX + t0
#         ZZ = ZZ * ZZ
#         Z2 = t1 * Z2 + b4 * ZZ
#
#         return X2, Z2
#
#     @staticmethod
#     def xADD(XP, ZP, XQ, ZQ, xPQ, zPQ, a, b):
#         r"""
#         Differential addition formula formula for x-only coordinate on short Weierstrass curve.
#
#         Cost: 7M + 2S + 2m0 + 6a
#
#         INPUT::
#
#             - ``XP``, ``ZP`` -- the coordinates of the Kummer point `x(P) = (XP : ZP)`
#
#             - ``XQ``, ``ZQ`` -- the coordinates of the Kummer point `x(Q) = (XQ : ZQ)`
#
#             - ``xPQ``, ``zPQ`` -- the coordinates of the Kummer point `x(P-Q) = (xPQ : zPQ)`
#
#             - ``a``, ``b`` -- curve constants
#
#         OUTPUT:: The coordinates of `x(P+Q) = (XQP : ZQP)`
#
#         EXAMPLES::
#
#             TODO
#         """
#         T1 = XP * XQ
#         T2 = ZP * ZQ
#         T3 = XP * ZQ
#         T4 = ZP * XQ
#         T5 = a * T2
#         T6 = T1 - T5
#         T7 = T6 * T6
#         T8 = b * T2
#         T9 = T8 + T8
#         T9 = T9 + T9
#         T10 = T3 + T4
#         T11 = T9 * T10
#         T12 = T7 - T11
#         XQP = zPQ * T12
#         T13 = T3 - T4
#         T13 = T13 * T13
#         ZQP = xPQ * T13
#
#         return XQP, ZQP
#
#     @staticmethod
#     def xDBLADD(XP, ZP, XQ, ZQ, xPQ, zPQ, a, b4):
#         r"""
#         Used in Montgomery ladder to do simultaneously doubling and differential addition.
#
#         Cost: TODO
#
#         INPUT::
#
#             - ``XP``, ``ZP`` -- the coordinates of the Kummer point `x(P) = (XP : ZP)`
#
#             - ``XQ``, ``ZQ`` -- the coordinates of the Kummer point `x(Q) = (XQ : ZQ)`
#
#             - ``xPQ``, ``zPQ`` -- the coordinates of the Kummer point `x(P-Q) = (xPQ : zPQ)`
#
#             - ``a``, ``b4`` -- curve constants where `b4 = 4*b`
#
#         OUTPUT:: The coordinates of `x(2P) = (X2P : Z2P)` and `x(P+Q) = (XQP : ZQP)`
#
#         EXAMPLES::
#
#             TODO
#         """
#
#         # TODO there is a weird edge case with these formulas if P = (0 : 1) and Q = (1 : 0)
#         # XX = XP**2
#         # ZZ = ZP**2
#         # aZZ = a * ZZ
#         # t0 = XP + ZP
#         # t1 = t0**2
#         # t2 = t1 - XX
#         # E = t2 - ZZ
#         # t3 = XX - aZZ
#         # t4 = t3**2
#         # t5 = E * ZZ
#         # t6 = b4 * t5
#         # X2P = t4 - t6
#         # t7 = XX + aZZ
#         # t8 = ZZ**2
#         # t9 = b4 * t8
#         # t10 = E * t7
#         # t11 = t10 + t10
#         # Z2P = t11 + t9
#         # A = XP * XQ
#         # B = ZP * ZQ
#         # C = XP * ZQ
#         # D = XQ * ZP
#         # t12 = a * B
#         # t13 = A - t12
#         # t14 = C + D
#         # t15 = t13**2
#         # t16 = B * t14
#         # t17 = b4 * t16
#         # t18 = t15 - t17
#         # XQP = zPQ * t18
#         # t19 = C - D
#         # t20 = t19**2
#         # ZQP = xPQ * t20
#         # print(XQP, ZQP)
#
#         X1, Z1 = xPQ, zPQ
#         X2, Z2 = XP, ZP
#         X3, Z3 = XQ, ZQ
#
#         XX = X2**2
#         ZZ = Z2**2
#         aZZ = a * ZZ
#         t0 = X2 + Z2
#         t1 = t0**2
#         t2 = t1 - XX
#         E = t2 - ZZ
#         t3 = XX - aZZ
#         t4 = t3**2
#         t5 = E * ZZ
#         t6 = b4 * t5
#         X4 = t4 - t6
#         t7 = XX + aZZ
#         t8 = ZZ**2
#         t9 = b4 * t8
#         t10 = E * t7
#         t11 = 2 * t10
#         Z4 = t11 + t9
#         A = X2 * X3
#         B = Z2 * Z3
#         C = X2 * Z3
#         D = X3 * Z2
#         t12 = a * B
#         t13 = C + D
#         t14 = A + t12
#         t15 = B**2
#         t16 = b4 * t15
#         t17 = t13 * t14
#         t18 = 2 * t17
#         R = t18 + t16
#         t19 = C - D
#         S = t19**2
#         t20 = S * X1
#         t21 = R * Z1
#         X5 = t21 - t20
#         Z5 = S * Z1
#
#         X2P, Z2P = X4, Z4
#         XQP, ZQP = X5, Z5
#
#         return X2P, Z2P, XQP, ZQP
#
#     def _double(self):
#         r""" """
#         X, Z = self.XZ()
#         a, b = self._parent.extract_constants()
#         b2 = b + b
#         b4 = b2 + b2
#         X2, Z2 = self.xDBL(X, Z, a, b2, b4)
#         return self._parent((X2, Z2))
#
#     def _double_iter(self, n):
#         r""" """
#         X, Z = self.XZ()
#         a, b = self._parent.extract_constants()
#         b2 = b + b
#         b4 = b2 + b2
#         for _ in range(n):
#             X, Z = self.xDBL(X, Z, a, b2, b4)
#         return self._parent((X, Z))
#
#     def double(self):
#         """
#         Wrapper function which deals with the doubling of
#         the identity
#
#         Returns [2] * self
#         """
#         # Deal with identity
#         if not self._Z:
#             return self
#         return self._double()
#
#     def double_iter(self, n):
#         """
#         Wrapper function which deals with the repeated
#         doubling
#
#         Returns [2^n] * self
#
#         This avoids the ADD part of xDBLADD, and so is
#         faster when we know our scalar is a power of two
#         """
#         # Deal with identity
#         if not self._Z:
#             return self
#         return self._double_iter(n)
#
#     def _add(self, Q, PQ):
#         """
#         Performs differential addition assuming
#         P, Q and PQ are all not the point at
#         infinity
#         """
#         XP, ZP = self.XZ()
#         XQ, ZQ = Q.XZ()
#         XPQ, ZPQ = PQ.XZ()
#
#         a, b = self._parent.extract_constants()
#
#         X_new, Z_new = self.xADD(XP, ZP, XQ, ZQ, XPQ, ZPQ, a, b)
#         return self._parent((X_new, Z_new))
#
#     def add(self, Q, PQ):
#         """
#         Function to perform differential addition and
#         handle the cases when P, Q or PQ are the points
#         at infinity
#         """
#         # Adding O + Q = Q
#         if not self._Z:
#             return Q
#
#         # Adding P + O = P
#         if not Q._Z:
#             return self
#
#         # Difference is the identity
#         # so P = Q and P+Q = [2]P
#         if not PQ._Z:
#             return self._double()
#
#         return self._add(Q, PQ)
#
#     def __mul__(self, m):
#         """
#         Montgomery-ladder to compute [m]P
#
#         Input: coordinates of P=(XP:ZP)
#                scalar factor m, curve constants (A:C)
#         Output: KummerPoint [m]P=(X0:Z0)
#         """
#         if not isinstance(m, (int, Integer)):
#             try:
#                 m = Integer(m)
#             except TypeError:
#                 raise TypeError(f"Cannot coerce input scalar {m = } to an integer")
#
#         # If m is zero, return identity
#         if not m:
#             return self.parent().zero()
#
#         # [m]P = [-m]P for x-only
#         m = abs(m)
#
#         # Extract base field and coefficients
#         R = self.base_ring()
#         XP, ZP = self.XZ()
#
#         # Initialise for loop
#         X0, Z0 = R(1), R(0)
#         X1, Z1 = XP, ZP
#
#         # Converting parameters for projective DBLADD -> (A24:C24)=(A+2C:4C)
#         a, b = self.parent().extract_constants()
#         b2 = b + b
#         b4 = b2 + b2
#         # A24 = C + C
#         # C24 = pari(A24 + A24)
#         # A24 = pari(A24 + A)
#
#         # Montgomery-ladder
#         for bit in bin(m)[2:]:
#             if bit == "0":
#                 # X0, Z0, X1, Z1 = self.xDBLADD(X0, Z0, X1, Z1, XP, ZP, a, b4)
#                 X1, Z1 = self.xADD(X0, Z0, X1, Z1, XP, ZP, a, b)
#                 X0, Z0 = self.xDBL(X0, Z0, a, b2, b4)
#             else:
#                 # X1, Z1, X0, Z0 = self.xDBLADD(X1, Z1, X0, Z0, XP, ZP, a, b4)
#                 X0, Z0 = self.xADD(X0, Z0, X1, Z1, XP, ZP, a, b)
#                 X1, Z1 = self.xDBL(X1, Z1, a, b2, b4)
#
#         return self._parent((X0, Z0))
#
#     def __rmul__(self, m):
#         return self * m
#
#     def __imul__(self, m):
#         self = self * m
#         return self
#
#     def ladder_3_pt(self, xP, xPQ, m):
#         """
#         Function to compute x(P + [m]Q) using x-only
#         arithmetic. Very similar to the Montgomery ladder above
#
#         Note: self = xQ
#         """
#         if not isinstance(m, (int, Integer)):
#             try:
#                 m = Integer(m)
#             except TypeError:
#                 raise TypeError(f"Cannot coerce input scalar {m = } to an integer")
#
#         # If m is zero, return xP
#         if not m:
#             return xP
#
#         # [m]P = [-m]P for x-only
#         m = abs(m)
#
#         # Converting parameters for projective DBLADD -> (A24:C24)=(A+2C:4C)
#         a, b = self.parent().extract_constants()
#         b4 = b + b
#         b4 = b4 + b4
#
#         # Extract out coordinates
#         XQ, ZQ = self.XZ()
#         XP, ZP = xP.XZ()
#         XPQ, ZPQ = xPQ.XZ()
#
#         # Montgomery-ladder
#         for bit in bin(m)[:1:-1]:
#             if bit == "1":
#                 XQ, ZQ, XP, ZP = self.xDBLADD(XQ, ZQ, XP, ZP, XPQ, ZPQ, a, b4)
#             else:
#                 XQ, ZQ, XPQ, ZPQ = self.xDBLADD(XQ, ZQ, XPQ, ZPQ, XP, ZP, a, b4)
#         return self._parent((XP, ZP))
#
#     def multiples(self):
#         """
#         A generator of points [l]P for self = P
#         Stops when it has generated the full subgroup generated by P
#         (without the identity point).
#
#         NOTE: this is implemented to make Vélu-like computations easy
#         """
#         yield self
#         R = self.double()
#         # Order 2 case
#         if not R:
#             return
#
#         # Odd order case
#         Q = self
#         while R:
#             yield R
#             S = R.add(self, Q)
#             Q, R = R, S
#
#         return
