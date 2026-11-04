r"""
Twisted Edwards curves

This module provides a native additive-group parent for a twisted Edwards
curve over a field.  A twisted Edwards curve is given by

.. MATH::

    a x^2 + y^2 = 1 + d x^2 y^2.

The implementation is deliberately separate from
:class:`~sage.schemes.elliptic_curves.ell_generic.EllipticCurve_generic`.
Sage's ``EllipticCurve`` constructor represents elliptic curves by a
Weierstrass model, whereas the Edwards model is particularly useful when the
coordinate system and its complete addition law are part of the application
(for example, Ed25519).

Only the complete twisted Edwards setting is accepted here: the base ring is
a field of characteristic different from 2, ``a`` and ``d`` are nonzero and
distinct, ``a`` is a square, and ``d`` is a nonsquare.  Under these
conditions the affine addition formulas are complete.  The points exposed by
this class are affine pairs ``(x, y)``.

This class models the curve group only.  The Ed25519 protocol layer, including
SHA-512 hashing, scalar pruning, little-endian encoding, the cofactor, and the
particular base point, is available in
:mod:`sage.crypto.ed25519`.

In the usual Ed25519 signature notation, ``r`` is the nonce scalar,
``R=[r]B`` is the resulting group point (normally transmitted in encoded
form), ``h`` is the challenge hash reduced modulo the subgroup order, and
``s = r + h a (mod L)`` is the response.  Thus ``R`` and ``r`` are not two
spellings for the same object.  A fault mask or a fault location in an
implementation is not part of the Ed25519 standard; it must be inferred from
the implementation or from the supplied faulty outputs.

EXAMPLES::

    sage: from sage.schemes.elliptic_curves.ell_edwards import TwistedEdwardsCurve
    sage: F = GF(19)
    sage: C = TwistedEdwardsCurve(F, 1, 2)
    sage: C
    Twisted Edwards curve over Finite Field of size 19 (a=1, d=2)
    sage: P = C(1, 0)
    sage: O = C(0, 1)
    sage: P + O == P
    True
    sage: 4*P == O
    True

The point operations use Sage's ordinary additive-group operators.  No
separate scalar-multiplication routine is required::

    sage: Q = C(0, -1)
    sage: P + P == Q
    True
    sage: -(P + Q) == -P - Q
    True
    sage: points = [C(x, y) for x in F for y in F if C.is_on_curve(x, y)]
    sage: all((P + Q) + R == P + (Q + R)
    ....:     for P in points for Q in points for R in points)
    True

The curve can be converted to the Weierstrass curve used by Sage's existing
elliptic-curve algorithms.  The conversion is birational, so the identity
and the exceptional 2-torsion point are handled explicitly::

    sage: E = C.to_elliptic_curve()
    sage: C.from_elliptic_curve(P.to_elliptic_curve()) == P
    True
    sage: C.from_elliptic_curve(O.to_elliptic_curve()) == O
    True
    sage: C.from_elliptic_curve(Q.to_elliptic_curve()) == Q
    True
    sage: all(C.from_elliptic_curve(P.to_elliptic_curve()) == P for P in points)
    True

REFERENCES:

- [Edwards2007]_
- [Hisil2008]_

.. [Edwards2007] H. M. Edwards, *A normal form for elliptic curves*,
   Bulletin of the American Mathematical Society 44 (2007), 393--422.

.. [Hisil2008] H. Hisil, K. K.-H. Wong, G. Carter, and E. Dawson,
   *Twisted Edwards Curves Revisited*, ASIACRYPT 2008, Lecture Notes in
   Computer Science 5350, 326--343.

AUTHORS:

- SageMath developers (2026): initial native twisted Edwards group support
"""

# ****************************************************************************
#       Copyright (C) 2026 SageMath developers
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.categories.commutative_additive_groups import CommutativeAdditiveGroups
from sage.categories.fields import Fields
from sage.misc.cachefunc import cached_method
from sage.schemes.elliptic_curves.constructor import EllipticCurve
from sage.structure.element import AdditiveGroupElement
from sage.structure.parent import Parent
from sage.structure.richcmp import op_EQ, op_NE
from sage.structure.unique_representation import UniqueRepresentation


class TwistedEdwardsPoint(AdditiveGroupElement):
    r"""A point on a :class:`TwistedEdwardsCurve`.

    The point is stored in the affine chart as a pair ``(x, y)``.  The
    parent validates the equation when ``check=True`` (the default).
    """

    def __init__(self, parent, coordinates, check=True):
        AdditiveGroupElement.__init__(self, parent)
        try:
            x, y = coordinates
        except (TypeError, ValueError):
            raise TypeError("a twisted Edwards point must have two coordinates")

        R = parent.base_ring()
        self._x = R(x)
        self._y = R(y)
        if check and not parent.is_on_curve(self._x, self._y):
            raise ValueError("the coordinates do not define a point on the curve")

    def _richcmp_(self, other, op):
        if op == op_EQ:
            return self._x == other._x and self._y == other._y
        if op == op_NE:
            return self._x != other._x or self._y != other._y
        return NotImplemented

    def __hash__(self):
        return hash((self.parent(), self._x, self._y))

    def __bool__(self):
        """Return whether this point is different from the identity."""
        return self._x != 0 or self._y != 1

    def __iter__(self):
        return iter((self._x, self._y))

    def __getitem__(self, index):
        return (self._x, self._y)[index]

    def _repr_(self):
        return "(%s : %s)" % (self._x, self._y)

    def coordinates(self):
        r"""Return the affine coordinates of this point."""
        return self._x, self._y

    def xy(self):
        r"""Return the affine coordinates of this point."""
        return self.coordinates()

    def x(self):
        r"""Return the first affine coordinate."""
        return self._x

    def y(self):
        r"""Return the second affine coordinate."""
        return self._y

    def curve(self):
        r"""Return the twisted Edwards curve containing this point."""
        return self.parent()

    def _add_(self, other):
        r"""Add two points using the complete affine formula.

        For ``P_i=(x_i,y_i)``, set

        .. MATH::

            \begin{aligned}
            x_3 &= \frac{x_1y_2+y_1x_2}
                         {1+d x_1x_2y_1y_2},\\
            y_3 &= \frac{y_1y_2-a x_1x_2}
                         {1-d x_1x_2y_1y_2}.
            \end{aligned}

        The parameter restrictions checked by the parent ensure that these
        denominators are nonzero for every pair of affine points.
        """
        C = self.parent()
        x1, y1 = self.coordinates()
        x2, y2 = other.coordinates()
        t = C.d() * x1 * x2 * y1 * y2
        dx = 1 + t
        dy = 1 - t
        if dx == 0 or dy == 0:
            raise ZeroDivisionError("the Edwards addition denominator is zero")
        x3 = (x1 * y2 + y1 * x2) / dx
        y3 = (y1 * y2 - C.a() * x1 * x2) / dy
        return C((x3, y3), check=False)

    def _neg_(self):
        return self.parent()((-self._x, self._y), check=False)

    def to_elliptic_curve(self):
        r"""Map this point to Sage's Weierstrass elliptic-curve model.

        The map is obtained through the Montgomery model

        .. MATH::

            Bv^2=u^3+Au^2+u,

        where ``u=(1+y)/(1-y)``, ``v=u/x``,
        ``A=2(a+d)/(a-d)``, and ``B=4/(a-d)``.  The corresponding
        Weierstrass coordinates are ``(X,Y)=(Bu,B^2v)``.

        The identity maps to the point at infinity.  The point ``(0,-1)``
        maps to ``(0,0)``.
        """
        C = self.parent()
        E = C.to_elliptic_curve()
        if not self:
            return E(0)
        if self._x == 0 and self._y == -1:
            return E(0, 0)

        B = C._montgomery_B()
        u = (1 + self._y) / (1 - self._y)
        v = u / self._x
        return E(B * u, B**2 * v)


class TwistedEdwardsCurve(UniqueRepresentation, Parent):
    r"""A complete twisted Edwards curve over a field.

    INPUT:

    - ``base_ring`` -- a field of characteristic different from 2

    - ``a``, ``d`` -- nonzero, distinct elements of ``base_ring`` such that
      ``a`` is a square and ``d`` is a nonsquare

    OUTPUT:

    The additive group of affine points satisfying

    .. MATH::

        a x^2 + y^2 = 1 + d x^2 y^2.

    The parameter restrictions are exactly those used by the complete
    twisted Edwards addition law.  For Ed25519 over
    ``GF(2^255 - 19)``, use ``a=-1`` and
    ``d=-121665/121666``.  The protocol-level Ed25519 encoding and signing
    operations are provided by :mod:`sage.crypto.ed25519`.

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.ell_edwards import TwistedEdwardsCurve
        sage: C = TwistedEdwardsCurve(GF(19), 1, 2)
        sage: C.a(), C.d()
        (1, 2)
        sage: C(0, 1).is_zero()
        True
        sage: C((0, -1))
        (0 : 18)

    Invalid models and points are rejected::

        sage: TwistedEdwardsCurve(ZZ, 1, 2)
        Traceback (most recent call last):
        ...
        TypeError: the base ring must be a field
        sage: TwistedEdwardsCurve(GF(19), 1, 1)
        Traceback (most recent call last):
        ...
        ValueError: ``a`` and ``d`` must be distinct
        sage: C(1, 1)
        Traceback (most recent call last):
        ...
        ValueError: the coordinates do not define a point on the curve
    """

    Element = TwistedEdwardsPoint

    @staticmethod
    def __classcall__(cls, base_ring, a, d):
        a = base_ring(a)
        d = base_ring(d)
        return super().__classcall__(cls, base_ring, a, d)

    def __init__(self, base_ring, a, d):
        if base_ring not in Fields():
            raise TypeError("the base ring must be a field")
        if base_ring.characteristic() == 2:
            raise ValueError("the base field must have characteristic different from 2")
        if a == 0 or d == 0:
            raise ValueError("``a`` and ``d`` must be nonzero")
        if a == d:
            raise ValueError("``a`` and ``d`` must be distinct")
        if not bool(a.is_square()):
            raise ValueError("``a`` must be a square in the base field")
        if bool(d.is_square()):
            raise ValueError("``d`` must be a nonsquare in the base field")

        self._a = a
        self._d = d
        Parent.__init__(self, base=base_ring, category=CommutativeAdditiveGroups())

    def _repr_(self):
        return "Twisted Edwards curve over %s (a=%s, d=%s)" % (
            self.base_ring(), self._a, self._d
        )

    def a(self):
        r"""Return the coefficient ``a``."""
        return self._a

    def d(self):
        r"""Return the coefficient ``d``."""
        return self._d

    def is_on_curve(self, x, y=None):
        r"""Return whether the given affine coordinates satisfy the equation."""
        if y is None:
            try:
                x, y = x
            except (TypeError, ValueError):
                return False
        try:
            x = self.base_ring()(x)
            y = self.base_ring()(y)
        except (TypeError, ValueError):
            return False
        return self._a * x**2 + y**2 == 1 + self._d * x**2 * y**2

    def _element_constructor_(self, x=0, y=None, check=True):
        if y is None:
            if isinstance(x, TwistedEdwardsPoint):
                if x.parent() is self:
                    return x
                x = x.coordinates()
            else:
                try:
                    if x == 0:
                        return self.zero()
                except (TypeError, ValueError):
                    pass
                try:
                    x, y = x
                except (TypeError, ValueError):
                    raise TypeError("a twisted Edwards point must have two coordinates")
        return self.element_class(self, (x, y), check=check)

    def _an_element_(self):
        return self(0, -1)

    def some_elements(self):
        points = [self.zero(), self(0, -1)]
        if self._a == 1:
            points.extend([self(1, 0), self(-1, 0)])
        return list(dict.fromkeys(points))

    @cached_method
    def zero(self):
        r"""Return the identity point ``(0, 1)``."""
        return self.element_class(self, (0, 1), check=False)

    def to_elliptic_curve(self):
        r"""Return the corresponding Weierstrass elliptic curve.

        If ``delta = a-d``, the returned curve is

        .. MATH::

            Y^2 = X^3 + ABX^2 + B^2X,

        with ``A=2(a+d)/delta`` and ``B=4/delta``.  Point conversion is
        provided by :meth:`TwistedEdwardsPoint.to_elliptic_curve` and
        :meth:`from_elliptic_curve`.
        """
        B = self._montgomery_B()
        A = self._montgomery_A()
        return EllipticCurve([0, A * B, 0, B**2, 0])

    def _montgomery_A(self):
        return 2 * (self._a + self._d) / (self._a - self._d)

    def _montgomery_B(self):
        return 4 / (self._a - self._d)

    def from_elliptic_curve(self, point):
        r"""Map a point on :meth:`to_elliptic_curve` back to this curve."""
        E = self.to_elliptic_curve()
        try:
            curve = point.curve()
        except AttributeError:
            raise TypeError("the argument must be a point on a Weierstrass curve")
        if curve != E:
            raise ValueError("the point must be on this curve's Weierstrass model")
        if point.is_zero():
            return self.zero()

        B = self._montgomery_B()
        u = point[0] / B
        v = point[1] / B**2
        if u == 0 and v == 0:
            return self(0, -1)
        if v == 0 or u == -1:
            raise ValueError("the point is outside the affine Edwards chart")
        y = (u - 1) / (u + 1)
        x = u / v
        return self(x, y)
