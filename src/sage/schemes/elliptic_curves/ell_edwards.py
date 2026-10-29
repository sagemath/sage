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

Scalar multiplication ``n * P`` is computed in extended homogeneous
coordinates ``(X : Y : Z : T)`` with ``x = X/Z``, ``y = Y/Z`` and
``T = XY/Z``, using the complete addition law of Hisil--Wong--Carter--Dawson.
Only a single field inversion is performed, when the affine coordinates of
the result are reconstructed.  Over prime fields ``GF(p)`` the arithmetic of
the scalar-multiplication inner loop runs directly on Python integers, and
points that are scalar-multiplied repeatedly build a cached fixed-base comb
table, which makes repeated multiplications such as ``n * B`` for the
Ed25519 base point substantially faster::

    sage: from sage.crypto.ed25519 import Ed25519, Ed25519BasePoint
    sage: B = Ed25519BasePoint()
    sage: k = 2^100 + 12345
    sage: (k * B).to_elliptic_curve() == k * B.to_elliptic_curve()
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
from sage.rings.infinity import infinity
from sage.rings.integer_ring import ZZ
from sage.schemes.elliptic_curves.constructor import EllipticCurve
from sage.structure.element import AdditiveGroupElement, Element
from sage.structure.parent import Parent
from sage.structure.richcmp import op_EQ, op_NE
from sage.structure.unique_representation import UniqueRepresentation

# Scalar-multiplication parameters.  One-shot multiplications use
# double-and-add with a small table of ``[j]P`` for ``j < 2^_EDWARDS_MUL_WINDOW``.
# Points that are scalar-multiplied ``_EDWARDS_COMB_THRESHOLD`` times over a
# prime field invest in a fixed-base comb table with window width
# ``_EDWARDS_COMB_WINDOW`` (radix 256, so a 255-bit scalar costs about 32
# mixed additions and no doublings).  ``_EDWARDS_COMB_MAX_BITS`` bounds the
# scalars served from the comb table; larger scalars fall back to the
# windowed method.
_EDWARDS_MUL_WINDOW = 4
_EDWARDS_COMB_WINDOW = 8
_EDWARDS_COMB_RADIX = 1 << _EDWARDS_COMB_WINDOW
_EDWARDS_COMB_MASK = _EDWARDS_COMB_RADIX - 1
_EDWARDS_COMB_MAX_BITS = 512
_EDWARDS_COMB_THRESHOLD = 2


def _edwards_ext_add(a, d, neg_a, one_a, P, Q):
    """Return the sum of two points in extended coordinates over any field.

    The unified addition law of Hisil--Wong--Carter--Dawson
    (add-2008-hwcd) is complete for the curves accepted by
    :class:`TwistedEdwardsCurve`, so no special cases are needed.
    ``P`` and ``Q`` are ``(X, Y, Z, T)`` tuples with ``T = XY/Z``.
    """
    X1, Y1, Z1, T1 = P
    X2, Y2, Z2, T2 = Q
    A = X1 * X2
    B = Y1 * Y2
    C = d * T1 * T2
    D = Z1 * Z2
    E = (X1 + Y1) * (X2 + Y2) - A - B
    F = D - C
    G = D + C
    if neg_a:
        H = B + A
    elif one_a:
        H = B - A
    else:
        H = B - a * A
    return (E * F, G * H, F * G, E * H)


def _edwards_ext_dbl(a, neg_a, one_a, P):
    """Return the double of ``P`` in extended coordinates (dbl-2008-hwcd)."""
    X1, Y1, Z1 = P[0], P[1], P[2]
    A = X1 * X1
    B = Y1 * Y1
    C = (Z1 + Z1) * Z1
    if neg_a:
        D = -A
    elif one_a:
        D = A
    else:
        D = a * A
    E = (X1 + Y1) * (X1 + Y1) - A - B
    G = D + B
    F = G - C
    H = D - B
    return (E * F, G * H, F * G, E * H)


def _edwards_ext_add_int(p, a, d, neg_a, one_a, P, Q):
    """Extended-coordinate addition on Python integers modulo the prime ``p``."""
    X1, Y1, Z1, T1 = P
    X2, Y2, Z2, T2 = Q
    A = X1 * X2 % p
    B = Y1 * Y2 % p
    C = T1 * d % p * T2 % p
    D = Z1 * Z2 % p
    E = ((X1 + Y1) * (X2 + Y2) - A - B) % p
    F = (D - C) % p
    G = (D + C) % p
    if neg_a:
        H = (B + A) % p
    elif one_a:
        H = (B - A) % p
    else:
        H = (B - a * A) % p
    return (E * F % p, G * H % p, F * G % p, E * H % p)


def _edwards_ext_dbl_int(p, a, neg_a, one_a, P):
    """Extended-coordinate doubling on Python integers modulo the prime ``p``."""
    X1, Y1, Z1 = P[0], P[1], P[2]
    A = X1 * X1 % p
    B = Y1 * Y1 % p
    C = (Z1 + Z1) * Z1 % p
    if neg_a:
        D = p - A
    elif one_a:
        D = A
    else:
        D = a * A % p
    E = ((X1 + Y1) * (X1 + Y1) - A - B) % p
    G = (D + B) % p
    F = (G - C) % p
    H = (D - B) % p
    return (E * F % p, G * H % p, F * G % p, E * H % p)


def _edwards_mul_windowed_int(point, k, p):
    """Windowed scalar multiplication of ``point`` on Python integers.

    ``k`` must be a positive integer.  Builds the table of small multiples
    ``[j]P``, ``1 <= j < 16``, and performs one inversion at the end.
    """
    C = point.parent()
    a = int(C.a().lift())
    d = int(C.d().lift())
    neg_a = (a == p - 1)
    one_a = (a == 1)
    x0 = int(point._x.lift())
    y0 = int(point._y.lift())
    base = (x0, y0, 1, x0 * y0 % p)

    k = int(k)
    radix = 1 << _EDWARDS_MUL_WINDOW
    digits = []
    while k:
        digits.append(k & (radix - 1))
        k >>= _EDWARDS_MUL_WINDOW
    digits.reverse()

    table = [None] * radix
    table[1] = base
    for j in range(2, radix):
        table[j] = _edwards_ext_add_int(p, a, d, neg_a, one_a, table[j - 1], base)

    acc = table[digits[0]]
    for digit in digits[1:]:
        for _ in range(_EDWARDS_MUL_WINDOW):
            acc = _edwards_ext_dbl_int(p, a, neg_a, one_a, acc)
        if digit:
            acc = _edwards_ext_add_int(p, a, d, neg_a, one_a, acc, table[digit])

    X, Y, Z, _ = acc
    zinv = pow(Z, -1, p)
    return C.element_class(C, (X * zinv % p, Y * zinv % p), check=False)


def _edwards_build_comb(point, p, bits):
    """Build the fixed-base comb table for ``point`` over ``GF(p)``.

    Returns ``(bits, rows, p, a, d, neg_a, one_a)`` where
    ``rows[i][j]`` stores the affine triple ``(x, y, x + y)`` together with
    ``d * x * y`` for the point ``[j * 16^i]P``.  All inversions are batched
    with Montgomery's trick, so building the table costs a single modular
    inversion per row.
    """
    C = point.parent()
    a = int(C.a().lift())
    d = int(C.d().lift())
    neg_a = (a == p - 1)
    one_a = (a == 1)
    x0 = int(point._x.lift())
    y0 = int(point._y.lift())

    positions = max(1, (bits + _EDWARDS_COMB_WINDOW - 1) // _EDWARDS_COMB_WINDOW)
    radix = 1 << _EDWARDS_COMB_WINDOW
    rows = []
    cur = (x0, y0, 1, x0 * y0 % p)
    for _ in range(positions):
        # Extended-coordinate multiples [j]([16^i]P) for j = 1..15.
        ext = [(0, 1, 1, 0)]
        for j in range(1, radix):
            ext.append(_edwards_ext_add_int(p, a, d, neg_a, one_a, ext[j - 1], cur))
        # Batch-invert the Z coordinates (never zero on a complete curve).
        zs = [entry[2] for entry in ext[1:]]
        prefix = []
        acc = 1
        for z in zs:
            acc = acc * z % p
            prefix.append(acc)
        invs = [0] * len(zs)
        inv = pow(acc, -1, p)
        for i in range(len(zs) - 1, -1, -1):
            invs[i] = inv * (prefix[i - 1] if i else 1) % p
            inv = inv * zs[i] % p
        row = [None]
        for (X, Y, Z, T), zinv in zip(ext[1:], invs):
            x = X * zinv % p
            y = Y * zinv % p
            row.append((x, y, (x + y) % p, d * (x * y % p) % p))
        rows.append(row)
        for _ in range(_EDWARDS_COMB_WINDOW):
            cur = _edwards_ext_dbl_int(p, a, neg_a, one_a, cur)
    return (bits, rows, p, a, d, neg_a, one_a)


def _edwards_comb_mul(point, k, comb):
    """Evaluate ``[k]point`` from a fixed-base comb table built by
    :func:`_edwards_build_comb`.

    Every nonzero base-16 digit of ``k`` contributes one mixed addition with
    an affine table entry, so no doublings are needed at evaluation time.
    """
    bits, rows, p, a, d, neg_a, one_a = comb
    k = int(k)
    X, Y, Z, T = 0, 1, 1, 0
    pos = 0
    while k:
        digit = k & _EDWARDS_COMB_MASK
        if digit:
            x2, y2, s2, c2 = rows[pos][digit]
            A = X * x2 % p
            B = Y * y2 % p
            CC = T * c2 % p
            E = ((X + Y) * s2 - A - B) % p
            F = (Z - CC) % p
            G = (Z + CC) % p
            if neg_a:
                H = (B + A) % p
            elif one_a:
                H = (B - A) % p
            else:
                H = (B - a * A) % p
            X = E * F % p
            Y = G * H % p
            Z = F * G % p
            T = E * H % p
        k >>= _EDWARDS_COMB_WINDOW
        pos += 1
    zinv = pow(Z, -1, p)
    C = point.parent()
    return C.element_class(C, (X * zinv % p, Y * zinv % p), check=False)


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
        denominators are nonzero for every pair of affine points.  The
        products ``x_1y_2`` and ``y_1x_2`` are shared with the denominator
        ``t = d x_1x_2y_1y_2`` to keep the number of field operations
        minimal.
        """
        C = self.parent()
        x1, y1 = self._x, self._y
        x2, y2 = other._x, other._y
        u = x1 * y2
        v = y1 * x2
        t = C.d() * u * v
        dx = 1 + t
        dy = 1 - t
        if dx == 0 or dy == 0:
            raise ZeroDivisionError("the Edwards addition denominator is zero")
        x3 = (u + v) / dx
        y3 = (y1 * y2 - C.a() * x1 * x2) / dy
        return type(self)(C, (x3, y3), check=False)

    def _neg_(self):
        C = self.parent()
        return type(self)(C, (-self._x, self._y), check=False)

    def _acted_upon_(self, other, self_on_left):
        r"""Return ``[other]self`` computed in extended Edwards coordinates.

        This hook provides the scalar multiplication ``n * P`` and ``P * n``
        for integer scalars, replacing the generic repeated-addition
        default of the coercion model.  The computation uses the complete
        extended-coordinate addition law, performs a single field inversion
        at the very end, and over prime fields caches a fixed-base comb
        table on points that are multiplied repeatedly.

        Scalars must coerce canonically into ``ZZ``: Python and Sage
        integers are accepted, while rational or finite-field elements are
        rejected with a :exc:`TypeError` rather than silently truncated to
        an integer (a rational action cached by the coercion model must not
        turn later fractional scalars into a silent ``None`` result)::

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.ell_edwards import TwistedEdwardsCurve
            sage: C = TwistedEdwardsCurve(GF(19), 1, 2)
            sage: P = C(1, 0)
            sage: 7 * P
            (18 : 0)
            sage: 7 * P == sum([P]*7, C.zero())
            True
            sage: all(k * P == sum([P]*int(k), C.zero()) for k in range(1, 20))
            True
            sage: 0 * P
            (0 : 1)
            sage: 1 * P is P
            True
            sage: (-3) * P == -(3 * P)
            True
            sage: P * 5 == 5 * P
            True

        The result agrees with scalar multiplication on the birationally
        equivalent Weierstrass model::

            sage: E = C.to_elliptic_curve()
            sage: all(k * P.to_elliptic_curve() == (k * P).to_elliptic_curve()
            ....:     for k in [0, 1, 2, 3, 7, 100, 2^60 + 5])
            True

        A curve with a generic parameter ``a`` takes the same path::

            sage: C17 = TwistedEdwardsCurve(GF(17), 4, 3)
            sage: Q = C17(3, 6)
            sage: all(k * Q == sum([Q]*int(k), C17.zero()) for k in range(1, 16))
            True
            sage: 12345 * Q == Q + 12344 * Q
            True

        Over fields that are not prime fields the same algorithm runs on
        the field's own elements::

            sage: F25.<z> = GF(25)
            sage: C25 = TwistedEdwardsCurve(F25, 1, z + 2)
            sage: P25 = C25(z + 3, 2)
            sage: 12 * P25 == C25.zero()
            True
            sage: (12345 * P25).to_elliptic_curve() == 12345 * P25.to_elliptic_curve()
            True

        A relative extension of prime degree is *not* a prime field; its
        points take the field-element path and still agree with the
        Weierstrass model::

            sage: K = F25.extension(PolynomialRing(F25, 't').gen() + 1, 'b')
            sage: CK = TwistedEdwardsCurve(K, 1, z + 2)
            sage: PK = CK(z + 3, 2)
            sage: (3 * PK).to_elliptic_curve() == 3 * PK.to_elliptic_curve()
            True

        Non-integer scalars are rejected::

            sage: QQ(1/2) * P
            Traceback (most recent call last):
            ...
            TypeError: scalar multiplication of a twisted Edwards point requires an integer scalar, got an element of Rational Field
            sage: GF(3)(2) * P
            Traceback (most recent call last):
            ...
            TypeError: scalar multiplication of a twisted Edwards point requires an integer scalar, got an element of Finite Field of size 3
        """
        try:
            k = ZZ.coerce(other)
        except TypeError:
            if not isinstance(other, Element):
                return None
            raise TypeError(
                "scalar multiplication of a twisted Edwards point requires "
                "an integer scalar, got an element of %s" % other.parent())
        if k == 1:
            return self
        if k == 0:
            return self.parent().zero()
        if k < 0:
            Q = self._edwards_mul(-k)._neg_()
        else:
            Q = self._edwards_mul(k)
        n = getattr(self, "_order", None)
        if n is not None and n != infinity:
            Q._order = n // n.gcd(k)
        return Q

    def additive_order(self):
        r"""Return the order of this point in the group law.

        The birational map :meth:`to_elliptic_curve` is a group
        isomorphism, so the order is computed on the Weierstrass model
        (where Sage's optimized algorithms are available) and cached on
        this point.  Points of order dividing 8 -- the cofactor torsion
        shared by all complete twisted Edwards curves -- are detected
        directly without touching the Weierstrass model, and scalar
        multiplications propagate a cached order to their result using
        Lagrange's theorem, mirroring the behaviour of Weierstrass
        points.  Infinite-order points over non-finite fields keep
        infinite order under nonzero multiplication.

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.ell_edwards import TwistedEdwardsCurve
            sage: C = TwistedEdwardsCurve(GF(19), 1, 2)
            sage: P = C(1, 0)
            sage: P.additive_order()
            4
            sage: P.order()
            4
            sage: P.order() is P.order()
            True
            sage: C.zero().order()
            1
            sage: C(0, -1).order()
            2
            sage: (3 * P).order()
            4
            sage: E = C.to_elliptic_curve()
            sage: all(Q.to_elliptic_curve().order() == Q.order()
            ....:     for Q in C.some_elements())
            True

        Over non-finite fields, points of infinite order stay of infinite
        order::

            sage: Cq = TwistedEdwardsCurve(QQ, 1, QQ(1/3))
            sage: Pq = Cq(2, 3)
            sage: Pq.order()
            +Infinity
            sage: (2 * Pq).order()
            +Infinity

        The Ed25519 base point has prime order ``L``::

            sage: from sage.crypto.ed25519 import ED25519_SUBGROUP_ORDER, Ed25519BasePoint
            sage: B = Ed25519BasePoint()
            sage: B.order() == ED25519_SUBGROUP_ORDER        # long time
            True
            sage: C25 = Ed25519()
            sage: C25(0, -1).order()
            2
        """
        o = self.__dict__.get("_order")
        if o is not None:
            return o

        if not self:
            o = ZZ.one()
        else:
            R = self + self
            if R.is_zero():
                o = ZZ(2)
            else:
                R = R + R
                if R.is_zero():
                    o = ZZ(4)
                else:
                    R = R + R
                    if R.is_zero():
                        o = ZZ(8)
                    else:
                        o = self.to_elliptic_curve().additive_order()
        self._order = o
        return o

    def _edwards_mul(self, k):
        r"""Return ``[k]self`` for an integer ``k > 1``.

        Dispatches between the generic field-element implementation and the
        Python-integer fast path for prime base fields, and manages the
        per-point fixed-base comb table.

        The comb table is built the second time the same point is
        multiplied, as a one-time cost: for a 255-bit scalar it takes on
        the order of 50 ms and a few megabytes of memory.  It pays off
        after roughly a hundred multiplications of the same point; one-shot
        multiplications of many different points stay on the windowed path,
        which costs about a millisecond per multiplication.
        """
        C = self.parent()
        p = C._edwards_prime_modulus()
        if p is None:
            return self._edwards_mul_field(k)
        bits = int(k).bit_length()
        if bits > _EDWARDS_COMB_MAX_BITS:
            return _edwards_mul_windowed_int(self, k, p)
        count = self.__dict__.get("_edwards_mult_count", 0) + 1
        self.__dict__["_edwards_mult_count"] = count
        if count < _EDWARDS_COMB_THRESHOLD:
            return _edwards_mul_windowed_int(self, k, p)
        comb = self.__dict__.get("_edwards_comb")
        if comb is None or comb[0] < bits:
            comb = _edwards_build_comb(self, p, bits)
            self.__dict__["_edwards_comb"] = comb
        return _edwards_comb_mul(self, k, comb)

    def _edwards_mul_field(self, k):
        r"""Return ``[k]self`` using the base field's own arithmetic.

        Windowed left-to-right scalar multiplication in extended
        coordinates with one inversion at the end.  Used over fields that
        are not prime fields.
        """
        C = self.parent()
        R = C.base_ring()
        a = C.a()
        d = C.d()
        one = R.one()
        neg_a = (a == -one)
        one_a = (a == one)
        x0, y0 = self._x, self._y
        base = (x0, y0, one, x0 * y0)

        digits = []
        while k:
            digits.append(k & 15)
            k >>= _EDWARDS_MUL_WINDOW
        digits.reverse()

        table = [None] * (1 << _EDWARDS_MUL_WINDOW)
        table[1] = base
        for j in range(2, 1 << _EDWARDS_MUL_WINDOW):
            table[j] = _edwards_ext_add(a, d, neg_a, one_a, table[j - 1], base)

        acc = table[digits[0]]
        for digit in digits[1:]:
            for _ in range(_EDWARDS_MUL_WINDOW):
                acc = _edwards_ext_dbl(a, neg_a, one_a, acc)
            if digit:
                acc = _edwards_ext_add(a, d, neg_a, one_a, acc, table[digit])

        X, Y, Z, _ = acc
        zinv = ~Z
        return type(self)(C, (X * zinv, Y * zinv), check=False)

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

    @cached_method
    def _edwards_prime_modulus(self):
        r"""Return the integer ``p`` if the base ring is ``GF(p)``, else ``None``.

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.ell_edwards import TwistedEdwardsCurve
            sage: TwistedEdwardsCurve(GF(19), 1, 2)._edwards_prime_modulus()
            19
            sage: TwistedEdwardsCurve(QQ, 1, 2)._edwards_prime_modulus() is None
            True

            sage: F25 = GF(25, 'z')
            sage: K = F25.extension(PolynomialRing(F25, 't').gen() + 1, 'b')
            sage: TwistedEdwardsCurve(K, 1, F25.gen() + 2)._edwards_prime_modulus() is None
            True
        """
        F = self.base_ring()
        try:
            # The integer fast path serves exactly the prime fields: the
            # cardinality must equal the characteristic.  A relative
            # extension of prime degree has degree() == 1 but is not
            # GF(p), so degree() alone is not enough here.
            p = F.characteristic()
            if F.cardinality() != p:
                return None
            return int(p)
        except (AttributeError, NotImplementedError, TypeError, ValueError):
            return None

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
