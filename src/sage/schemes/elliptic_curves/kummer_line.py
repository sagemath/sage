r"""
Kummer Line for short Weierstrass elliptic curves

Sage defines Kummer lines derived from an elliptic curve in short Weierstrass or Montgomery form, known as `x`-only arithmetic.
Also defines a point on a Kummer line as an element of the projective line as well as isogenies between Kummer lines.

EXAMPLES:

We construct a Kummer line over an already defined elliptic curve in short Weierstrass form::

    sage: E = EllipticCurve(GF(101), [2, 3])
    sage: KummerLine(E)
    Kummer line of Elliptic curve defined by y^2 = x^3 + 2*x + 3 over Finite Field of size 101

We can also give the curve coefficients directly with a base ring::

    sage: KummerLine(QQ, [4, 5/6])
    Kummer line of Elliptic curve defined by y^2 = x^3 + 4*x + 5/6 over Rational Field

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

#
# ::
#
# TODO KummerLinePoint examples, as well as many more usage (see integer for instance)
# TODO <Lots and lots of examples>


from sage.all import cached_method

from sage.rings.integer import Integer
from sage.schemes.elliptic_curves.ell_generic import EllipticCurve_generic
from sage.schemes.elliptic_curves.constructor import EllipticCurve
from sage.schemes.elliptic_curves.ell_point import EllipticCurvePoint_field
from sage.structure.element import RingElement
from sage.structure.sage_object import SageObject


class KummerLine(SageObject):
    r"""
    Kummer line class for a short Weierstrass elliptic curve.

    EXAMPLES::

        sage: E = EllipticCurve(GF(101), [2, 3])
        sage: K = KummerLine(E); K
        Kummer line of Elliptic curve defined by y^2 = x^3 + 2*x + 3 over Finite Field of size 101
    """

    def __init__(self, curve):
        r"""
        Constructor for a Kummer line from a short Weierstrass elliptic curve.

        INPUT::

            - ``curve`` -- an elliptic curve in short Weierstrass form `y^2 = x^3 + a*x + b`

        EXAMPLES::

            sage: E = EllipticCurve(GF(101), [2, 3])
            sage: K = KummerLine(E); K
            Kummer line of Elliptic curve defined by y^2 = x^3 + 2*x + 3 over Finite Field of size 101
        """

        if not isinstance(curve, EllipticCurve_generic):
            raise TypeError("input must be an elliptic curve in short Weierstrass form")

        ainvs = curve.a_invariants()
        a, b = ainvs[3], ainvs[4]  # a = a4, b = a6
        if ainvs != (0, 0, 0, a, b):
            raise ValueError(
                "input must be an elliptic curve in short Weierstrass form"
            )
        # TODO could deal with Montgomery model using an isomorphism

        self._curve = curve
        self._base_ring = curve.base_ring()

        # Initialize variables
        self._a = self._base_ring(a)
        self._b = self._base_ring(b)

        self._montgomery = False

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
            'Kummer line of Elliptic curve defined by y^2 = x^3 + 2*x + 3 over Finite Field of size 101'
        """
        return f"Kummer line of {self._curve}"

    def __call__(self, coords):
        r"""
        Create a Kummer point from the coordinates.

        INPUT::

            - ``coords`` - either a point `P` on the elliptic curve or a valid tuple `(X, Z)` where `P = (X : * : Z)`; `Z` is optional

        OUTPUT: A Kummer point on the Kummer line

        EXAMPLES::

            sage: E = EllipticCurve(GF(101), [2, 3])
            sage: P = E(95, 49)
            sage: K = KummerLine(E)

        Point as a parameter::

            sage: xP = K(P); xP
            (95 : 1)

        The output is the same if we just give the `XZ`-coordinates::

            sage: K([95, 1]) == K(P)
            True

        `Z` is optional::

            sage: K(95) == K(P)
            True
            sage: K([95]) == K(P)
            True
        """
        return KummerLinePoint(self, coords)

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
        Return the short Weierstrass coefficients `(a, b)` as a tuple.

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
        Return the associated elliptic curve.

        EXAMPLES::

            sage: E = EllipticCurve(GF(101), [2, 3])
            sage: K = KummerLine(E); K.curve()
            Elliptic Curve defined by y^2 = x^3 + 2*x + 3 over Finite Field of size 101
        """
        return self._curve

    def is_montgomery(self):
        r"""
        Return if the associated curve is a Montgomery one or not.

        EXAMPLES::

            sage: E = EllipticCurve(GF(101), [2, 3])
            sage: K = KummerLine(E); K.is_montgomery()
            False
        """
        return self._montgomery

    @cached_method
    def j_invariant(self):
        r"""
        Compute the j-invariant of the Kummer line using the formula `j(E) = 1728\frac{4a^3}{4a^3 + 27b^2}`.

        Cached method.

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
        Compute the discriminant of the Kummer line using the formula `\Delta(E) = -16(4a^3 + 27b^2)`.

        Cached method.

        EXAMPLES::

            sage: K = KummerLine(GF(101), [2, 3])
            sage: K.discriminant()
            44
        """
        a, b = self.extract_constants()
        return -16 * (4 * a**3 + 27 * b**2)

    def isogeny(self, xP):
        r"""
        Compute the isogeny with kernel ``xP``.

        EXAMPLES::

            sage: E = EllipticCurve(GF(101), [2, 3]); K = KummerLine(E)
            sage: P = E(76, 65); xP = K(P)
            sage: Q = E(49, 61); xQ = K(Q)
            sage: KummerLineIsogeny(xP)
            Isogeny between Kummer lines of degree 12 with kernel (76 : 1)
            Domain: Kummer line of Elliptic curve defined by y^2 = x^3 + 2*x + 3 over Finite Field of size 101
            Codomain: Kummer line of Elliptic curve defined by y^2 = x^3 + 23*x + 73 over Finite Field of size 101
        """

        if xP._parent != self:
            raise ValueError("point is not on the Kummer line")
        return KummerLineIsogeny(xP)


class KummerLinePoint(SageObject):
    r"""
    Kummer point on a Kummer line.

    EXAMPLES::

        sage: E = EllipticCurve(GF(101), [2, 3])
        sage: K = KummerLine(E)
        sage: P = E(95, 52)

    A Kummer point can be constructed from `XZ`-coordinates::

        sage: xP = KummerLinePoint(K, [95, 1]); xP
        (95 : 1)

    Or, it can be constructed from a point on the elliptic curve::

        sage: xP = KummerLinePoint(K, P); xP
        (95 : 1)

    The point can also be created by calling the Kummer line::

        sage: xP = K(P); xP
        (95 : 1)
    """

    def __init__(self, parent, coords):
        r"""
        Create a point from the coordinates.

        INPUT::

            - ``parent`` - A Kummer line

            - ``coords`` - either a point P on EllipticCurve or a tuple (X, Z) where `P = (X : * : Z)`; Z is optional

        EXAMPLES::

            sage: E = EllipticCurve(GF(101), [2, 3])
            sage: P = E(95, 52)
            sage: K = KummerLine(E)

        A Kummer point can be constructed from `XZ`-coordinates::

            sage: xP = KummerLinePoint(K, [95, 1]); xP
            (95 : 1)

        Z-coordinate is optional, assumed to be 1 then::

            sage: KummerLinePoint(K, 95) == KummerLinePoint(K, [95, 1])
            True
            sage: KummerLinePoint(K, 95) == KummerLinePoint(K, [95, 2])
            False
            sage: KummerLinePoint(K, 95) == KummerLinePoint(K, [95])
            True

        Or, it can be constructed from a point on the elliptic curve::

            sage: xP = KummerLinePoint(K, P); xP
            (95 : 1)
        """
        # Ensure the parent is the right type
        if not isinstance(parent, KummerLine):
            raise TypeError("not a Kummer line")

        R = parent.base_ring()

        # Point at infinity
        if coords is None:
            coords = (R(1), R(0))

        # Construct point from P on an elliptic curve in Montgomery form
        elif isinstance(coords, EllipticCurvePoint_field):
            # Make sure point's parent curve matches with Kummer Line
            if coords.is_zero():
                coords = (R(1), R(0))
            else:
                a, b = parent.extract_constants()
                # assert coords.curve().a_invariants() == (0, 0, 0, a, b)
                # Done this way to handle Montgomery curves
                assert coords.curve().is_isomorphic(EllipticCurve(R, [a, b]))
                coords = coords[0], coords[2]

        # Construct from X coordinate only
        elif isinstance(coords, RingElement):
            coords = (coords,)

        # Construct from a tuple (X : Z)
        else:
            coords = tuple(coords)

        # Sanitise the input coordinates
        if len(coords) == 1:
            coords += (R(1),)
        if len(coords) != 2:
            raise ValueError("not a point on the projective line P^1")
        coords = tuple(map(R, coords))

        self._base_ring = R
        self._parent = parent
        self._X, self._Z = coords

        if not self.is_x_coord():
            raise ValueError("not a valid x-coordinate")

    def __repr__(self):
        r"""
        String representation of a Kummer point.

        EXAMPLES::

            sage: E = EllipticCurve(GF(101), [2, 3])
            sage: P = E(95, 49)
            sage: K = KummerLine(E)
            sage: xP = K(P); xP.__repr__()
        """
        return f"({self._X} : {self._Z})"

    def __bool__(self):
        r"""
        Boolean value for a Kummer line point.

        OUTPUT: ``False`` if this is the point at infinity, otherwise ``True``

        .. SEEALSO::

            meth:`is_zero`

        EXAMPLES::

            sage: E = EllipticCurve(GF(101), [2, 3])
            sage: P = E(95, 49)
            sage: K = KummerLine(E)
            sage: xO = K(None); bool(xO)
            False
            sage: xP = K(P); bool(xP)
            False
        """

        return bool(self._Z)

    def __eq__(self, other):
        r"""
        Equality of two Kummer points as `\PP^1` elements.

        EXAMPLES::

            sage: E1 = EllipticCurve(GF(101), [2, 3])
            sage: K1 = KummerLine(E1)
            sage: P = E1(95, 49)

        Two points can have different coordinates and still be the same projectively::
            sage: xP = K1(P1)
            sage: xQ = K1([29, 12])
            sage: xP == xQ  # 95*12 % 101 = 29
            True

        They must represent the same `x`-coordinate::
            sage: Q = E1(50, 41)
            sage: xP, xQ = K1(P), K1(Q)
            sage: xP == xQ
            False

        Base Kummer line must be the same even if the `x`-coordinates are equal::
            sage: E2 = EllipticCurve(GF(101), [3, 7])
            sage: Q = E2(95, 50)
            sage: K2 = KummerLine(E2)
            sage: xP, xQ = K1(P), K2(Q)
            sage: xP == xQ
            False
            sage: xP.x() == xQ.x()
            True
            sage: xP.parent() == xQ.parent()
            False
        """
        if not isinstance(other, KummerLinePoint):
            raise ValueError("can only compare equality between two Kummer points")

        if self.parent() != other.parent():
            return False

        XP, ZP = self.XZ()
        XQ, ZQ = other.XZ()
        return XP * ZQ == XQ * ZP

    def is_zero(self):
        r"""
        A Kummer Point is zero if it corresponds to the point at infinity
        on the parent curve.

        OUTPUT: ``False`` if this is the point at infinity, otherwise ``True``

        EXAMPLES::

            sage: E = EllipticCurve(GF(101), [2, 3])
            sage: P = E(95, 49)
            sage: K = KummerLine(E)
            sage: xO = K(None); xO.is_zero()
            False
            sage: xP = K(P); xP.is_zero()
            False
        """

        return not bool(self)

    def is_x_coord(self):
        r"""
        Check if the `x`-coordinate of a Kummer point is a valid
        one on the parent curve.

        Raise an error in initialisation if ``False``.

        .. TODO::

            Dealing with a twist of the curve eventually

        EXAMPLES::

            sage: E = EllipticCurve(GF(101), [2, 3])
            sage: P = E(95, 49)
            sage: K = KummerLine(E)
            sage: xP = K(95); xP.is_x_coord()
            True
        """
        if self.is_zero():
            return True
        return self.parent().curve().is_x_coord(self.x())

    def base_ring(self):
        r"""
        Return the base ring of the Kummer Point.

        EXAMPLES::

            sage: E = EllipticCurve(GF(101), [2, 3])
            sage: P = E(95, 49)
            sage: K = KummerLine(E)
            sage: xP = K(P); xP.base_ring()
            Finite Field of size 101
        """
        return self._base_ring

    def parent(self):
        r"""
        Return the Kummer line on which the point is constructed.

        EXAMPLES::

            sage: E = EllipticCurve(GF(101), [2, 3])
            sage: P = E(95, 49)
            sage: K = KummerLine(E)
            sage: xP = K(P); xP.parent()
            Kummer line of Elliptic curve defined by y^2 = x^3 + 2*x + 3 over Finite Field of size 101
        """
        return self._parent

    def XZ(self):
        r"""
        Return the projective coordinates `(X : Z)` of the point.

        EXAMPLES::

            sage: E = EllipticCurve(GF(101), [2, 3])
            sage: P = E(95, 49)
            sage: K = KummerLine(E)
            sage: xP = K(P)
            sage: xP.XZ()
            (95, 1)
        """
        return self._X, self._Z

    def x(self):
        r"""
        Return the affine `x`-coordinate of the point, if well-defined.

        Raise an error on the point at infinity.

        EXAMPLES::

            sage: E = EllipticCurve(GF(101), [2, 3])
            sage: P = E(95, 49)
            sage: K = KummerLine(E)
            sage: xP = K(P); xP.x()
            95

        The point at infinity has no valid `x`-coordinate::

            sage: xO = K(None); xO.x()
            ValueError: the point at infinity has no valid x-coordinate
        """
        if self.is_zero():
            raise ValueError("the point at infinity has no valid x-coordinate")

        X, Z = self.XZ()

        if Z == 1:
            return self._base_ring(X)
        return self._base_ring(X / Z)

    @cached_method
    def curve_point(self):
        r"""
        Return a lift of the point on the parent curve.

        `y`-coordinate is chosen as the smallest square root.
        Cached method.

        EXAMPLES::

            sage: E = EllipticCurve(GF(101), [2, 3])
            sage: P = E(95, 49)
            sage: K = KummerLine(E)
            sage: xP = K(P); xP.curve_point()
            (95 : 49 : 1)
        """
        K = self.parent()
        E = K.curve()

        if self.is_zero():
            return E(0)

        x = self.x()
        return E.lift_x(x)

    @staticmethod
    def xDBL(X, Z, a, b2, b4):
        r"""
        Doubling formula for `x`-only coordinate on short Weierstrass curve.

        Cost: ``2M + 5S + 3m0 + 8a``

        INPUT:

            - ``X``, ``Z`` -- the coordinates of the Kummer point `x(P)`

            - ``a``, ``b2``, ``b4`` -- curve constants where ``b2 = 2*b`` and ``b4 = 4*b``

        OUTPUT: ``(X2, Z2)``, the coordinates of `x(2P)`

        EXAMPLES::

            sage: E = EllipticCurve(GF(101), [2, 3])
            sage: K = KummerLine(E)
            sage: P = E(95, 49); xP = K(P)
            sage: a, b = K.extract_constants()
            sage: b2, b4 = 2*b, 4*b
            sage: X, Z = xP.XZ()
            sage: X2, Z2 = xP.xDBL(X, Z, a, b2, b4)
            sage: if Z2:
            sage:     (2*P)[0] == X2/Z2
            sage: else:
            sage:     P == -P
            True
        """
        XX = X * X
        ZZ = Z * Z
        t1 = X + Z
        t1 = t1 * t1
        t1 = t1 - XX
        t1 = t1 - ZZ
        A = t1 + t1
        aZZ = a * ZZ
        t1 = XX - aZZ
        t1 = t1 * t1
        t2 = A * ZZ
        t2 = b2 * t2
        X2 = t1 - t2
        t1 = XX + aZZ
        t2 = ZZ * ZZ
        t1 = A * t1
        t2 = b4 * t2
        Z2 = t1 + t2

        return X2, Z2

    @staticmethod
    def xADD(XP, ZP, XQ, ZQ, XPQ, ZPQ, a, b):
        r"""
        Differential addition formula formula for `x`-only coordinate on short Weierstrass curve.

        Cost: ``7M + 2S + 2m0 + 6a``

        INPUT:

        - ``XP``, ``ZP`` -- the coordinates of the Kummer point `x(P)`

        - ``XQ``, ``ZQ`` -- the coordinates of the Kummer point `x(Q)`

        - ``XPQ``, ``ZPQ`` -- the coordinates of the Kummer point `x(P-Q)`

        - ``a``, ``b`` -- curve constants

        OUTPUT: ``(XQP, ZQP)``, the coordinates of `x(P+Q)`

        .. NOTE::

            The formula does not behave well when `x(P-Q) = (1 : 0)` (should use `xDBL` instead)

        EXAMPLES::

            sage: E = EllipticCurve(GF(101), [2, 3])
            sage: K = KummerLine(E)
            sage: P, Q = E(95, 49), E(30, 46)
            sage: xP, xQ, xPQ = K(P), K(Q), K(P-Q)
            sage: a, b = K.extract_constants()
            sage: (XP, ZP), (XQ, ZQ), (XPQ, ZPQ) = xP.XZ(), xQ.XZ(), xPQ.XZ()
            sage: XQP, ZQP = xP.xADD(XP, ZP, XQ, ZQ, XPQ, ZPQ, a, b)
            sage: if ZQP:
            sage:     (P+Q)[0] == XQP/ZQP
            sage: else:
            sage:     P == -Q
            True
        """

        t1 = XP * XQ
        t2 = ZP * ZQ
        t3 = XP * ZQ
        t4 = ZP * XQ
        t5 = a * t2
        t1 = t1 - t5
        t1 = t1 * t1
        t2 = b * t2
        t2 = t2 + t2
        t2 = t2 + t2
        t5 = t3 + t4
        t5 = t2 * t5
        t1 = t1 - t5
        XQP = ZPQ * t1
        t1 = t3 - t4
        t1 = t1 * t1
        ZQP = XPQ * t1

        return XQP, ZQP

    @staticmethod
    def xDBLADD(XP, ZP, XQ, ZQ, XPQ, ZPQ, a, b4):
        r"""
        Differential addition and doubling for `x`-only coordinate on short Weierstrass curve.

        Cost: ``9M + 7S + 5m0 + 12a``

        INPUT:

        - ``XP``, ``ZP`` -- the coordinates of the Kummer point `x(P)`

        - ``XQ``, ``ZQ`` -- the coordinates of the Kummer point `x(Q)`

        - ``XPQ``, ``ZPQ`` -- the coordinates of the Kummer point `x(P-Q)`

        - ``a``, ``b4`` -- curve constants where ``b4 = 4*b``

        OUTPUT: ``(X2P, Z2P)`` and ``(XQP, ZQP)``, the coordinates of `x(2P)` and `x(P+Q)` respectively

        .. NOTE::

            The formula does not behave well when `x(P-Q) = (1 : 0)` (should use `xDBL` instead)

        EXAMPLES::

            sage: E = EllipticCurve(GF(101), [2, 3])
            sage: K = KummerLine(E)
            sage: P, Q = E(95, 49), E(30, 46)
            sage: xP, xQ, xPQ = K(P), K(Q), K(P-Q)
            sage: a, b = K.extract_constants()
            sage: b4 = 4*b
            sage: (XP, ZP), (XQ, ZQ), (XPQ, ZPQ) = xP.XZ(), xQ.XZ(), xPQ.XZ()
            sage: X2P, Z2P, XQP, ZQP = xP.xDBLADD(XP, ZP, XQ, ZQ, XPQ, ZPQ, a, b4)
            sage: if ZQP:
            sage:     (P+Q)[0] == XQP/ZQP
            sage: else:
            sage:     P == -Q
            True
            sage: if Z2P:
            sage:     (2*P)[0] == X2P/Z2P
            sage: else:
            sage:     P == -P
            True
        """
        XX = XP * XP
        ZZ = ZP * ZP
        aZZ = a * ZZ
        t1 = XP + ZP
        t1 = t1 * t1
        t1 = t1 - XX
        E = t1 - ZZ
        t1 = XX - aZZ
        t1 = t1 * t1
        t2 = E * ZZ
        t2 = b4 * t2
        X2P = t1 - t2
        t1 = XX + aZZ
        t1 = E * t1
        t1 = t1 + t1
        t2 = ZZ * ZZ
        t2 = b4 * t2
        Z2P = t1 + t2
        A = XP * XQ
        B = ZP * ZQ
        C = XP * ZQ
        D = XQ * ZP
        t1 = a * B
        t1 = A - t1
        t1 = t1 * t1
        t2 = C + D
        t2 = B * t2
        t2 = b4 * t2
        t1 = t1 - t2
        XQP = ZPQ * t1
        t1 = C - D
        t1 = t1 * t1
        ZQP = XPQ * t1

        return X2P, Z2P, XQP, ZQP

    def _double(self):
        r"""
        Return the doubling of the point.

        .. SEEALSO::

            meth:``double``

        EXAMPLES::

            sage: E = EllipticCurve(GF(101), [2, 3])
            sage: K = KummerLine(E)
            sage: P = E(95, 49); xP = K(P)
            sage: xP._double()
            (88 : 9)
        """
        X, Z = self.XZ()
        a, b = self._parent.extract_constants()
        b2 = b + b
        b4 = b2 + b2
        X, Z = self.xDBL(X, Z, a, b2, b4)
        return self._parent((X, Z))

    def _double_iter(self, n):
        r"""
        Return the repeated doubling of the point `n` times.

        .. SEEALSO::

            meth:``double_iter``

        EXAMPLES::

            sage: E = EllipticCurve(GF(101), [2, 3])
            sage: K = KummerLine(E)
            sage: P = E(95, 49); xP = K(P)
            sage: xP._double_iter(42)
            (2 : 52)
        """
        X, Z = self.XZ()
        a, b = self._parent.extract_constants()
        b2 = b + b
        b4 = b2 + b2
        for _ in range(n):
            X, Z = self.xDBL(X, Z, a, b2, b4)
        return self._parent((X, Z))

    def double(self):
        r"""
        Return the doubling of the point.

        EXAMPLES::

            sage: E = EllipticCurve(GF(101), [2, 3])
            sage: K = KummerLine(E)
            sage: P = E(95, 49); xP = K(P)
            sage: xP.double()
            (88 : 9)
        """

        # Deal with identity
        if self.is_zero():
            return self
        return self._double()

    def double_iter(self, n):
        r"""
        Return the repeated doubling of the point `n` times.

        This avoids the ``ADD`` part of ``xDBLADD``, hence is faster
        when the scalar is known to be a power of `2`.

        EXAMPLES::

            sage: E = EllipticCurve(GF(101), [2, 3])
            sage: K = KummerLine(E)
            sage: P = E(95, 49); xP = K(P)
            sage: xP.double_iter(42)
            (2 : 52)
        """
        if self.is_zero():
            return self
        return self._double_iter(n)

    def _add(self, Q, PQ):
        r"""
        Return the differential addition of `x(P)`, `x(Q)` and `x(P-Q)` assuming none is the point at infinity.

        INPUT:

        - ``self`` -- the Kummer point `x(P)`

        - ``Q`` -- the Kummer point `x(Q)`

        - ``PQ`` -- the Kummer point `x(P-Q)`

        .. SEEALSO::

            meth:``add``

        EXAMPLES::

            sage: E = EllipticCurve(GF(101), [2, 3])
            sage: K = KummerLine(E)
            sage: P, Q = E(95, 49), E(30, 46)
            sage: xP, xQ, xPQ = K(P), K(Q), K(P-Q)
            sage: xP._add(xQ, xPQ)
            (11 : 40)
        """
        XP, ZP = self.XZ()
        XQ, ZQ = Q.XZ()
        XPQ, ZPQ = PQ.XZ()

        a, b = self._parent.extract_constants()

        XQP, ZQP = self.xADD(XP, ZP, XQ, ZQ, XPQ, ZPQ, a, b)
        return self._parent((XQP, ZQP))

    def add(self, Q, PQ):
        r"""
        Return the differential addition of `x(P)`, `x(Q)` and `x(P-Q)`.

        INPUT:

        - ``self`` -- the Kummer point `x(P)`

        - ``Q`` -- the Kummer point `x(Q)`

        - ``PQ`` -- the Kummer point `x(P-Q)`

        EXAMPLES::

            sage: E = EllipticCurve(GF(101), [2, 3])
            sage: K = KummerLine(E)
            sage: P, Q = E(95, 49), E(30, 46)
            sage: xP, xQ, xPQ = K(P), K(Q), K(P-Q)
            sage: xP.add(xQ, xPQ)
            (11 : 40)
        """
        # Adding O + Q = Q
        if self.is_zero():
            return Q

        # Adding P + O = P
        if Q.is_zero():
            return self

        # P-Q = O so P+Q = 2P
        if PQ.is_zero():
            return self._double()

        return self._add(Q, PQ)

    def __mul__(self, m):
        r"""
        Montgomery-ladder to compute `x(mP)`.

        INPUT: `m` -- scalar to multiply `P` by.

        EXAMPLES::

            sage: E = EllipticCurve(GF(101), [2, 3])
            sage: K = KummerLine(E)
            sage: P = E(95, 49); xP = K(P)
            sage: m = 3820495
            sage: (m*P)[0] == (m*xP).x()
            True
        """
        if not isinstance(m, (int, Integer)):
            try:
                m = Integer(m)
            except TypeError:
                raise TypeError(f"cannot coerce input scalar {m = } to an integer")

        # If m is zero, return identity
        if not m:
            return self.parent().zero()

        # [m]P = [-m]P for x-only
        m = abs(m)

        # Compute necessary constants
        a, b = self._parent.extract_constants()
        b2 = b + b
        b4 = b2 + b2

        # Initialise for loop
        XP, ZP = self.XZ()
        R = self._parent.base_ring()
        X0, Z0 = R(1), R(0)
        X1, Z1 = XP, ZP

        # Montgomery ladder
        for bit in bin(m)[2:]:
            if bit == "0":
                X0, Z0, X1, Z1 = self.xDBLADD(X0, Z0, X1, Z1, XP, ZP, a, b4)
            else:
                X1, Z1, X0, Z0 = self.xDBLADD(X1, Z1, X0, Z0, XP, ZP, a, b4)

        return self._parent((X0, Z0))

    def __rmul__(self, m):
        r"""
        Multiplication by a scalar on the right side.

        .. SEEALSO::

            meth:``__mul__``

        EXAMPLES::

            sage: E = EllipticCurve(GF(101), [2, 3])
            sage: K = KummerLine(E)
            sage: P = E(95, 49); xP = K(P)
            sage: m = 3820495
            sage: xP*m == m*xP
            True
        """
        return self * m

    def __imul__(self, m):
        r"""
        Inplace multiplication by a scalar.

        .. SEEALSO::

            meth:``__mul__``

        EXAMPLES::

            sage: E = EllipticCurve(GF(101), [2, 3])
            sage: K = KummerLine(E)
            sage: P = E(95, 49); xP = K(P)
            sage: m = 3820495
            sage: xP *= m; xP
            (69 : 33)
        """
        self = self * m
        return self

    def ladder_3_pt(self, xQ, xPQ, m):
        """
        Compute `x(P + mQ)` using `x`-only arithmetic.

        Very similar to the Montgomery ladder.

        INPUT:

        - ``self`` -- the Kummer point `x(P)`

        - ``Q`` -- the Kummer point `x(Q)`

        - ``PQ`` -- the Kummer point `x(P-Q)`

        - ``m`` -- integer

        EXAMPLES::

            sage: E = EllipticCurve(GF(101), [2, 3])
            sage: K = KummerLine(E)
            sage: P, Q = E(95, 49), E(30, 46)
            sage: xP, xQ, xPQ = K(P), K(Q), K(P-Q)
            sage: m = 29432
            sage: xR = xP.ladder_3_pt(xQ, xPQ, m)
            sage: xR.x() == (P + m*Q)[0]
            True
        """
        if not isinstance(m, (int, Integer)):
            try:
                m = Integer(m)
            except TypeError:
                raise TypeError(f"cannot coerce input scalar {m = } to an integer")

        # If m is zero, return xP
        if not m:
            return self

        # [m]P = [-m]P for x-only
        m = abs(m)

        # Compute necessary constants
        a, b = self._parent.extract_constants()
        b4 = b + b
        b4 = b4 + b4

        # Extract out coordinates
        XP, ZP = self.XZ()
        XQ, ZQ = xQ.XZ()
        XPQ, ZPQ = xPQ.XZ()

        # Montgomery ladder
        for bit in bin(m)[:1:-1]:
            if bit == "1":
                XQ, ZQ, XP, ZP = self.xDBLADD(XQ, ZQ, XP, ZP, XPQ, ZPQ, a, b4)
            else:
                XQ, ZQ, XPQ, ZPQ = self.xDBLADD(XQ, ZQ, XPQ, ZPQ, XP, ZP, a, b4)
        return self._parent((XP, ZP))

    # TODO ditch it or improve it, there could be less points in this
    def multiples(self):
        r"""
        Generator for points `x(mP)` where ``P = self``.

        Stops when the full subgroup is generated, minus the point at infinity

        .. NOTE::

            Usage is intended for Vélu-like computations

        EXAMPLES::

            sage: E = EllipticCurve(GF(101), [2, 3])
            sage: K = KummerLine(E)
            sage: P = E(95, 49); xP = K(P)
            sage: L = list(xP.multiples()); len(L)
            47 # P has order 48 on E
            sage: L[36] == 37*xP
            True
        """
        yield self
        R = self.double()
        # Order 2 case
        if R.is_zero():
            return

        # Odd order case
        Q = self
        while R:
            yield R
            S = R.add(self, Q)
            Q, R = R, S

        return


class KummerLineIsogeny(SageObject):
    r"""
    Isogeny between Kummer lines using Vélu formulas.

    .. TODO::

        Optimize formulas: VéluSqrt, use projective coordinates, ...

    EXAMPLES::

        sage: E = EllipticCurve(GF(101), [2, 3]); K = KummerLine(E)
        sage: P = E(76, 65); xP = K(P)
        sage: Q = E(49, 61); xQ = K(Q)
        sage: KummerLineIsogeny(xP)
        Isogeny between Kummer lines of degree 12 with kernel (76 : 1)
        Domain: Kummer line of Elliptic curve defined by y^2 = x^3 + 2*x + 3 over Finite Field of size 101
        Codomain: Kummer line of Elliptic curve defined by y^2 = x^3 + 23*x + 73 over Finite Field of size 101

    It can also be constructed with the ``isogeny`` method from a Kummer line::

        sage: K.isogeny(xP)
        Isogeny between Kummer lines of degree 12 with kernel (76 : 1)
        Domain: Kummer line of Elliptic curve defined by y^2 = x^3 + 2*x + 3 over Finite Field of size 101
        Codomain: Kummer line of Elliptic curve defined by y^2 = x^3 + 23*x + 73 over Finite Field of size 101
    """

    def __init__(self, kernel):
        r"""
        Create an isogeny with a given kernel on a Kummer line.

        EXAMPLES::

            sage: E = EllipticCurve(GF(101), [2, 3]); K = KummerLine(E)
            sage: P = E(76, 65); xP = K(P)
            sage: Q = E(49, 61); xQ = K(Q)
            sage: KummerLineIsogeny(xP)
            Isogeny between Kummer lines of degree 12 with kernel (76 : 1)
            Domain: Kummer line of Elliptic curve defined by y^2 = x^3 + 2*x + 3 over Finite Field of size 101
            Codomain: Kummer line of Elliptic curve defined by y^2 = x^3 + 23*x + 73 over Finite Field of size 101
        """

        # Ensure the kernel is the right type
        if not isinstance(kernel, KummerLinePoint):
            raise TypeError("kernel is not a Kummer point")

        self._kernel = kernel
        self._domain = kernel.parent()
        # Corresponds to the set "S" used in Vélu's formulas
        self._kernel_subset, self._degree = self._compute_kernel_subset_and_degree()
        self._kernel_x = [xQ.x() for xQ in self._kernel_subset]

        (
            self._velu_constants_t,
            self._velu_constants_u,
            self._codomain,
        ) = self._compute_velu_constants_and_codomain()

    def degree(self):
        r"""
        Return the degree of the isogeny.

        EXAMPLES::

            sage: E = EllipticCurve(GF(101), [2, 3]); K = KummerLine(E)
            sage: P = E(76, 65); xP = K(P)
            sage: xf = K.isogeny(xP); xf.degree()
            12
        """
        return self._degree

    def domain(self):
        r"""
        Return the domain of the isogeny.

        EXAMPLES::

            sage: E = EllipticCurve(GF(101), [2, 3]); K = KummerLine(E)
            sage: P = E(76, 65); xP = K(P)
            sage: xf = K.isogeny(xP); xf.domain()
            Kummer line of Elliptic curve defined by y^2 = x^3 + 2*x + 3 over Finite Field of size 101
        """
        return self._domain

    def codomain(self):
        r"""
        Return the codomain of the isogeny.

        EXAMPLES::

            sage: E = EllipticCurve(GF(101), [2, 3]); K = KummerLine(E)
            sage: P = E(76, 65); xP = K(P)
            sage: xf = K.isogeny(xP); xf.codomain()
            Kummer line of Elliptic curve defined by y^2 = x^3 + 23*x + 73 over Finite Field of size 101
        """
        return self._codomain

    def kernel(self):
        r"""
        Return the kernel generator of the isogeny.

        EXAMPLES::

            sage: E = EllipticCurve(GF(101), [2, 3]); K = KummerLine(E)
            sage: P = E(76, 65); xP = K(P)
            sage: xf = K.isogeny(xP); xf.kernel()
            (76 : 1)
        """
        return self._kernel

    def __repr__(self):
        r"""
        String representation of an isogeny between Kummer lines.

        EXAMPLES::

            sage: E = EllipticCurve(GF(101), [2, 3]); K = KummerLine(E)
            sage: P = E(76, 65); xP = K(P)
            sage: xf = K.isogeny(xP); xf.__repr__()
            Isogeny between Kummer lines of degree 12 with kernel (76 : 1)
            Domain: Kummer line of Elliptic curve defined by y^2 = x^3 + 2*x + 3 over Finite Field of size 101
            Codomain: Kummer line of Elliptic curve defined by y^2 = x^3 + 23*x + 73 over Finite Field of size 101
        """
        return f"Isogeny between Kummer lines of degree {self._degree} with kernel {self._kernel}\nDomain: {self._domain}\nCodomain: {self._codomain}"

    def __call__(self, xP):
        r"""
        Evaluate the isogeny at ``xP``.

        .. SEEALSO::

            meth:`evaluate`

        EXAMPLES::

            sage: E = EllipticCurve(GF(101), [2, 3]); K = KummerLine(E)
            sage: P = E(76, 65); xP = K(P)
            sage: Q = E(49, 61); xQ = K(Q)
            sage: f = E.isogeny(P)
            sage: xf = K.isogeny(xP) # Isogeny on K with kernel xP
            sage: f.codomain() == xf.codomain().curve()
            True
            sage: f(Q)[0] == xf(Q).x()
            True
        """
        return self.evaluate(xP)

    @cached_method
    def _compute_kernel_subset_and_degree(self):
        r"""
        Compute the subset of multiples of the kernel needed as well as the degree.

        This first compute all the multiples of the generator until reaching
        point at infinity, then returns half of them.
        Called during initialisation.
        Cached method.
        """
        L = list(self._kernel.multiples())
        deg = len(L) + 1
        return L[: (len(L) + 1) // 2], deg

    @cached_method
    def _compute_velu_constants_and_codomain(self):
        r"""
        Compute the needed constants in Vélu formulas as well as the codomain.

        Called during initialisation.
        Cached method.
        """
        R = self._domain.base_ring()
        a, b = self._domain.extract_constants()
        t, w = R(0), R(0)
        cst_t, cst_u = [], []
        S = self._kernel_subset
        Sx = self._kernel_x

        for xQx in Sx[:-1]:
            tQ = 2 * (3 * xQx**2 + a)
            uQ = 4 * (xQx**3 + a * xQx + b)
            cst_t.append(tQ)
            cst_u.append(uQ)
            t += tQ
            w += uQ + xQx * tQ

        xQ = S[-1]
        xQx = Sx[-1]
        gQx = 3 * xQx**2 + a
        uQ = 4 * (xQx**3 + a * xQx + b)
        if xQ.double().is_zero():
            tQ = gQx
        else:
            tQ = 2 * gQx
        cst_t.append(tQ)
        cst_u.append(uQ)
        t += tQ
        w += uQ + xQx * tQ

        A = a - 5 * t
        B = b - 7 * w

        E2 = EllipticCurve(R, [A, B])
        K2 = KummerLine(E2)

        return cst_t, cst_u, K2

    def evaluate(self, xP):
        r"""
        Evaluate the isogeny at ``xP``.

        EXAMPLES::

            sage: E = EllipticCurve(GF(101), [2, 3]); K = KummerLine(E)
            sage: P = E(76, 65); xP = K(P)
            sage: Q = E(49, 61); xQ = K(Q)
            sage: f = E.isogeny(P)
            sage: xf = K.isogeny(xP) # Isogeny on K with kernel xP
            sage: f.codomain() == xf.codomain().curve()
            True
            sage: f(Q)[0] == xf.evaluate(Q).x()
            True
        """
        if xP in self._kernel_subset:
            return self._codomain.zero()

        Sx = self._kernel_x
        xPx = xP.x()
        cst_t, cst_u = self._velu_constants_t, self._velu_constants_u

        fxP = xPx
        for i in range(len(Sx)):
            xQx, tQ, uQ = Sx[i], cst_t[i], cst_u[i]
            alpha = xPx - xQx
            fxP += (tQ * alpha + uQ) / alpha**2

        return self._codomain(fxP)
