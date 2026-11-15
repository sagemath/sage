r"""
Arithmetic on the Jacobian

This module implements the group operation in the Picard group of a
hyperelliptic curve, represented as divisors in Mumford
representation, using Cantor's algorithm.

AUTHORS:

- Sabrina Kunzweiler, Gareth Ma, Giacomo Pope (2024): adapt to smooth model
"""

# ****************************************************************************
#       Copyright (C) 2025 Sabrina Kunzweiler <sabrina.kunzweiler@math.u-bordeaux.fr>
#                     2025 Gareth Ma <grhkm21@gmail.com>
#                     2025 Giacomo Pope <giacomopope@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from __future__ import annotations

from typing import TYPE_CHECKING, Self

if TYPE_CHECKING:
    from collections.abc import Iterator

    from sage.rings.integer import Integer

from sage.groups.generic import order_from_multiple
from sage.misc.cachefunc import cached_method
from sage.rings.finite_rings.finite_field_base import FiniteField as FiniteField_generic
from sage.categories.fields import Fields
from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.polynomial_element import Polynomial
from sage.schemes.generic.morphism import SchemeMorphism
from sage.structure.element import AdditiveGroupElement
from sage.structure.richcmp import richcmp


class MumfordDivisorClassField(AdditiveGroupElement, SchemeMorphism):
    r"""
    An element of a Jacobian defined over a field, i.e. in
    `J(K) = \mathrm{Pic}^0_K(C)`.
    """

    def __init__(self, parent, u, v, check=True) -> None:
        """
        Create an element of the Jacobian of a hyperelliptic curve.

        EXAMPLES::

            sage: R.<x> = GF(19)[]
            sage: H = HyperellipticCurve(x^5 + x, x^2 + 1)
            sage: J = H.jacobian()
            sage: D = J(x^2 + 14*x + 16, 3*x + 4); D # indirect doctest
            (x^2 + 14*x + 16, 3*x + 4)
        """
        SchemeMorphism.__init__(self, parent)

        if not isinstance(u, Polynomial) or not isinstance(v, Polynomial):
            raise TypeError(f"arguments u={u} and v={v} must be polynomials")

        # TODO:
        # 1. allow elements of the base field as input
        #   (in particular something like (u,v) = (x-alpha, 0))

        # Ensure the divisor is valid
        if check:
            f, h = parent.curve().hyperelliptic_polynomials()
            assert (v**2 + v * h - f) % u == 0, (
                f"u={u}, v={v} do not define a divisor on the Jacobian"
            )

            # TODO: should we automatically do reduction here if the degree of u is
            # too large?

        self._parent = parent
        self._u: Polynomial = u
        self._v: Polynomial = v

    def scheme(self):
        r"""
        Return the scheme this morphism maps to; or, where this divisor lives.

        .. WARNING::

           Although a pointset is defined over a specific field, the
           scheme returned may be over a different (usually smaller)
           field.  The example below demonstrates this: the pointset
           is determined over a number field of absolute degree 2 but
           the scheme returned is defined over the rationals.

        EXAMPLES::

            sage: x = QQ['x'].gen()
            sage: f = x^5 + x
            sage: H = HyperellipticCurve(f)
            sage: F.<a> = NumberField(x^2 - 2, 'a')                                     # needs sage.rings.number_field
            sage: J = H.jacobian()(F); J                                                # needs sage.rings.number_field
            Abelian group of points over Number Field in a with defining polynomial x^2 - 2 on Jacobian of Hyperelliptic Curve over Rational Field defined by y^2 = x^5 + x


        ::

            sage: P = J(H.lift_x(F(1)))                                                 # needs sage.rings.number_field
            sage: P.scheme()                                                            # needs sage.rings.number_field
            Jacobian of Hyperelliptic Curve over Rational Field defined by y^2 = x^5 + x
        """
        return self.codomain()

    def _repr_(self) -> str:
        r"""
        Return the Mumford presentation of the divisor class.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(QQ)
            sage: H = HyperellipticCurve(x^5 + x^2 + x + 1)
            sage: P = H([0,-1])
            sage: Q = H([1,-2])
            sage: J = H.jacobian()
            sage: D = J(P,Q); D # indirect doctest
            (x^2 - x, 3*x - 1)
        """
        return f"({self._u}, {self._v})"

    def _latex_(self) -> str:
        r"""
        Return the Mumford presentation of the divisor class in LaTeX.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(QQ)
            sage: H = HyperellipticCurve(x^5 + x^2 + x + 1)
            sage: P = H([0,-1])
            sage: Q = H([1,-2])
            sage: J = H.jacobian()
            sage: D = J(P, Q)
            sage: latex(D)
            (x^{2} - x, 3 x - 1)
        """
        return f"({self._u._latex_()}, {self._v._latex_()})"

    def uv(self) -> tuple[Polynomial, Polynomial]:
        r"""
        Return the `u` and `v` component of this Mumford divisor.

        EXAMPLES::

            sage: x = polygen(GF(1993))
            sage: H = HyperellipticCurve(x^7 + 3 * x + 1)
            sage: J = H.jacobian(); J
            Jacobian of Hyperelliptic Curve over Finite Field of size 1993 defined by y^2 = x^7 + 3*x + 1
            sage: u, v = x^3 + 1570*x^2 + 1930*x + 81, 368*x^2 + 1478*x + 256
            sage: P = J(u, v); P
            (x^3 + 1570*x^2 + 1930*x + 81, 368*x^2 + 1478*x + 256)
            sage: P.uv() == (u, v)
            True
        """
        return (self._u, self._v)

    def _richcmp_(self, other, op) -> bool:
        r"""
        Method for rich comparison.

        TESTS::

            sage: x = polygen(GF(23))
            sage: H = HyperellipticCurve(x^7 + x + 1)
            sage: J = H.jacobian(); J
            Jacobian of Hyperelliptic Curve over Finite Field of size 23 defined by y^2 = x^7 + x + 1
            sage: P = J.random_element()
            sage: P == P
            True
            sage: P != P
            False

            sage: Z = J(0); Z
            (1, 0)
            sage: Z == 0
            True
            sage: Z != 0
            False
        """
        return richcmp(tuple(self), tuple(other), op)

    def __iter__(self) -> Iterator[Polynomial | int | Integer]:
        r"""
        TESTS:

        Indirect tests::

            sage: R.<x> = GF(13)[]
            sage: H = HyperellipticCurve(x^5 - x + 1)
            sage: J = H.jacobian()
            sage: P = J(x + 5, 12)
            sage: list(P)
            [x + 5, 12]
            sage: tuple(P)
            (x + 5, 12)
        """
        yield from [self._u, self._v]

    def __getitem__(self, n) -> Polynomial:
        r"""
        Return the n-th item in the Mumford presentation of the divisor.

        TESTS::

            sage: R.<x> = GF(23)[]
            sage: H = HyperellipticCurve(x^5 + x + 1, x^2 + x +  1)
            sage: J = H.jacobian()
            sage: D = J(x^2 + 21*x + 10, 14*x + 22)
            sage: D[0] # indirect doctest
            x^2 + 21*x + 10
            sage: D[1] # indirect doctest
            14*x + 22
        """
        return (self._u, self._v)[n]

    def __hash__(self) -> int:
        r"""
        Compute the hash value of this element.

        TESTS::

            sage: R.<x> = GF(23)[]
            sage: H = HyperellipticCurve(x^5 + x + 1, x^2 + x +  1)
            sage: J = H.jacobian()
            sage: D = J(x^2 + 21*x + 10, 14*x + 22)
            sage: hash(D) == hash(2*D)
            False
        """
        return hash(tuple(self))

    def __reduce__(self) -> tuple:
        r"""
        Return data for pickling.

        TESTS::

            sage: x = polygen(GF(5))
            sage: H = HyperellipticCurve(x^5 + 3*x + 1)
            sage: J = H.jacobian()
            sage: P = J.random_element()
            sage: loads(dumps(P)) == P
            True
        """
        return (self._parent, tuple(self))

    def __bool__(self) -> bool:
        r"""
        Return ``True`` if this is not the zero element of the Jacobian.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: H = HyperellipticCurve(x^2 + x, x^3 + 1)
            sage: J = H.jacobian()
            sage: bool(J(x+1,0))
            True
            sage: bool(J.zero())
            False

        ::

            sage: x = polygen(GF(419))
            sage: C = HyperellipticCurve(x^5 + x)
            sage: bool(C.jacobian().zero())
            False
        """
        return not (self._u.is_one() and self._v.is_zero())

    @cached_method
    def order(self) -> Integer:
        r"""
        Returns the order of self.
        This is only implemented over finite fields.

        EXAMPLES::

            sage: K = FiniteField(7)
            sage: R.<x> = K[]
            sage: H = HyperellipticCurve(x^6 + 3*x + 2)
            sage: JK = Jacobian(H)(K)
            sage: D = JK(x^2 + 5*x + 6, 6*x + 3)
            sage: D.order()
            38
        """
        if not isinstance(self.base_ring(), FiniteField_generic):
            raise NotImplementedError(
                "this is only implemented for Jacobians over a finite field"
            )
        n = self.parent().order()
        return order_from_multiple(self, n)

    def degree(self) -> Integer:
        r"""
        Returns the degree of the affine part of the divisor.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: H = HyperellipticCurve(6*x^5 + 9*x^4 - x^3 - 3*x^2, 1)
            sage: J = H.jacobian()
            sage: J.zero().degree()
            0
            sage: J(x, 0).degree()
            1
            sage: J(x^2 + 1/2*x, -x - 1).degree()
            2
        """
        return self._u.degree()

    def _add_(self, other) -> Self:
        r"""
        Add `other` to the divisor `self`.

        TESTS::

            sage: R.<x> = QQ[]
            sage: H = HyperellipticCurve(x^5 - x^4 + x^2 - x, 1)
            sage: J = H.jacobian()
            sage: D1 = J(x^2 + x, x)
            sage: D2 = J(x^2, x - 1)
            sage: D1 + D2
            (x^2 + x, -1)

            sage: Jp = J.change_ring(GF(13))
            sage: D1, D2, D3 = [Jp.random_element() for i in range(3)]
            sage: D1 + D2 + D3 == D3 + D1 + D2
            True
        """
        # Collect data from HyperellipticCurve
        H = self.parent().curve()
        g = H.genus()

        # Extract out mumford coordinates
        u1, v1 = self.uv()
        u2, v2 = other.uv()

        # Step one: cantor composition of the two divisors
        u3, v3 = self._parent.cantor_composition(u1, v1, u2, v2)

        # Step two: cantor reduction of the above to ensure
        # the degree of u is smaller than g + 1
        while u3.degree() > g:
            u3, v3 = self._parent.cantor_reduction(u3, v3)
        v3 %= u3

        return self._parent(u3, v3, check=False)

    def _neg_(self) -> Self:
        r"""
        Multiply the divisor by `[-1]`.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: H = HyperellipticCurve(x^5 + 2*x^4 + 3*x^3 + 2*x^2 + x, 1)
            sage: P = H(1,0,0)
            sage: Q = H(0,0)
            sage: J = H.jacobian()
            sage: D = J(P,Q); D
            (x, -1)
            sage: -D #indirect doctest
            (x, 0)
            sage: J(Q,P) == - D #indirect doctest
            True
        """
        _, h = self.parent().curve().hyperelliptic_polynomials()
        u0, v0 = self.uv()

        if h.is_zero():
            return self._parent(u0, -v0, check=False)

        v1 = (-h - v0) % u0
        return self._parent(u0, v1, check=False)

    # ------------------------------------------------------------------
    # Miller functions and pairings
    # ------------------------------------------------------------------

    def _miller_add(self, other) -> tuple:
        r"""
        Return the sum ``self + other`` together with a function `h`
        on the curve such that

        .. MATH::

            \mathrm{div}(h) = D_{\mathrm{self}} + D_{\mathrm{other}} - D_{\mathrm{self} + \mathrm{other}},

        where `D_A` denotes the degree-zero divisor encoded by the
        Mumford representation of `A` (including its part at infinity).

        The function `h` is returned as a list of triples ``(a, b, c)`` of
        polynomials, each representing the function `(a(x) + b(x) y) / c(x)`,
        which give `h` when multiplied together.

        EXAMPLES::

            sage: R.<x> = GF(13)[]
            sage: H = HyperellipticCurve(x^5 + 2*x + 1)
            sage: J = H.jacobian()
            sage: D1 = J(x^2 + 1, 10*x + 6)
            sage: D2 = J(x + 5, 5)
            sage: D3, h = D1._miller_add(D2)
            sage: D3 == D1 + D2
            True
            sage: h
            [(x + 5, 0, 1)]
            sage: D3, h = D1._miller_add(D1)
            sage: D3 == 2*D1
            True
            sage: h
            [(1, 0, 1), (4*x^2 + 3*x + 11, 1, x + 10)]
        """
        parent = self._parent
        g = parent.curve().genus()
        R = self._u.parent()

        u1, v1 = self.uv()
        u2, v2 = other.uv()

        # Composition (affine parts): D1 + D2 = D3 + div(s)
        u3, v3, s = parent._cantor_composition_generic(u1, v1, u2, v2)
        funs = [(s, R.zero(), R.one())]

        # Reduction (affine parts): D3 - D4 = div((y - v3) / u4)
        while u3.degree() > g:
            u4, v4 = parent.cantor_reduction(u3, v3)
            funs.append((-v3, R.one(), u4))
            u3, v3 = u4, v4
        v3 %= u3

        return parent(u3, v3, check=False), funs

    def _miller_values(self, n, Ds) -> list:
        r"""
        Return the values `f_{n,P}(D)` for each `D` in ``Ds``, where
        `P` is ``self`` and `f_{n,P}` is the Miller function with divisor
        `n D_P - D_{nP}`.

        Here `D_A` denotes the divisor encoded by the Mumford representation
        of `A`. The divisors `D` are effective affine divisors given by
        their Mumford coordinates ``(u, v)``. Since the Miller function is
        only well-defined up to a constant, the output is only meaningful for
        combinations of degree zero.

        This method raises a :exc:`ZeroDivisionError` if one of the divisors
        meets the support of one of the intermediate functions.

        EXAMPLES::

            sage: R.<x> = GF(101)[]
            sage: H = HyperellipticCurve(x^5 + 2*x + 1)
            sage: J = H.jacobian()
            sage: P = J(x^2 + x + 8, 80*x + 51)
            sage: P.order()
            2399
            sage: E1 = J(x^2 + 78*x + 14, 51*x + 38)
            sage: E2 = J(x^2 + 74*x + 90, 54*x + 53)
            sage: P._miller_values(2399, [E1.uv(), E2.uv()])
            [62, 33]

        Since `P` itself is in the support of `f_{n,P}`, we cannot
        evaluate the Miller function at `P`::

            sage: P._miller_values(2399, [P.uv()])
            Traceback (most recent call last):
            ...
            ZeroDivisionError: divisor meets the support of a Miller function
        """
        n = ZZ(n)
        if n <= 0:
            raise ValueError("n must be a positive integer")

        K = self._parent.value_ring()
        Ds = [(u.monic(), v) for u, v in Ds]
        num = [K.one()] * len(Ds)
        den = [K.one()] * len(Ds)

        def update(funs):
            for i, (u, v) in enumerate(Ds):
                if u.degree() == 0:
                    continue
                for a, b, c in funs:
                    num[i] *= u.resultant((a + b * v) % u)
                    den[i] *= u.resultant(c % u)
            if not all(num + den):
                raise ZeroDivisionError("divisor meets the support of a Miller function")

        T = self
        for bit in reversed(n.bits()[:-1]):
            num = [z**2 for z in num]
            den = [z**2 for z in den]
            T, funs = T._miller_add(T)
            update(funs)
            if bit:
                T, funs = T._miller_add(self)
                update(funs)

        return [a / b for a, b in zip(num, den)]

    def _random_auxiliary_divisor(self, other=None, tries=10):
        r"""
        Return a random element `R` of the Jacobian such that the Mumford
        coordinates `u` of `R` and of ``self + R`` both have degree `g`.

        For such `R`, the divisor class of ``self`` is represented by the
        affine degree-zero divisor `E_{\mathrm{self} + R} - E_R`, where
        `E_A` denotes the effective affine divisor defined by the Mumford
        coordinates of `A`. If ``other`` is an element of the parent, we
        also require that the supports of `E_{\mathrm{self} + R}` and `E_R`
        are disjoint from those of `E_{\mathrm{other}}`.

        OUTPUT: `R`

        A :exc:`ValueError` is raised if no such `R` is found after
        ``tries`` random samples, which may happen over very small fields.

        EXAMPLES::

            sage: R.<x> = GF(13)[]
            sage: H = HyperellipticCurve(x^5 + 2*x + 1)
            sage: J = H.jacobian()
            sage: P = J(x + 5, 5)
            sage: PR, R = P._random_auxiliary_divisor()
            sage: PR - R == P
            True
            sage: PR[0].degree() == R[0].degree() == 2
            True

        Over `\FF_2`, the Jacobian below has no element with Mumford
        coordinate `u` of degree `2`, so there is no suitable auxiliary
        divisor::

            sage: R.<x> = GF(2)[]
            sage: H = HyperellipticCurve(x^5 + 1, x^2 + 1)
            sage: J = H.jacobian()
            sage: J.points()
            [(1, 0), (x + 1, 0)]
            sage: J.zero()._random_auxiliary_divisor()
            Traceback (most recent call last):
            ...
            ValueError: unable to find a suitable auxiliary divisor

        TESTS:

        This curve has no rational points; check that we don't hang::

            sage: K.<a> = GF(4)
            sage: R.<x> = K[]
            sage: H = HyperellipticCurve(x^6 + a*x^4 + a*x^3 + a*x + a, a*x^3 + 1)
            sage: H.is_inert()
            True
            sage: H.rational_points()
            []
            sage: J = H.jacobian()(K)
            sage: J.zero()._random_auxiliary_divisor()
            Traceback (most recent call last):
            ...
            ValueError: unable to find a suitable auxiliary divisor
        """
        parent = self._parent
        H = parent.extended_curve()
        g = H.genus()
        for _ in range(tries):
            try:
                R = parent.random_element(fast=True)
            except ValueError:
                # curve is pointless
                break
            if R._u.degree() != g:
                continue
            PR = self + R
            if PR._u.degree() != g:
                continue
            if other is not None and not other._u.gcd(PR._u * R._u).is_one():
                continue
            return PR, R
        raise ValueError("unable to find a suitable auxiliary divisor")

    def change_ring(self, R):
        r"""
        Return ``self`` in the group of points of the Jacobian
        over a different ring.

        The argument ``R`` can either be a ring or a coercion
        map into another ring.

        EXAMPLES::

            sage: R.<x> = GF(13)[]
            sage: H = HyperellipticCurve(x^6 + x + 1)
            sage: J = H.jacobian()
            sage: P = J(x + 3, 8, 0)
            sage: L, phi = GF(13).extension(3, map=True)
            sage: PL = P.change_ring(L); PL
            (x + 3, 8 : 0)
            sage: PL.parent()
            Abelian group of points over Finite Field in z3 of size 13^3
              on Jacobian of Hyperelliptic Curve over Finite Field of size 13
              defined by y^2 = x^6 + x + 1
            sage: PL = P.change_ring(phi); PL
            (x + 3, 8 : 0)
            sage: PL.parent()
            Abelian group of points over Finite Field in z3 of size 13^3
              on Jacobian of Hyperelliptic Curve over Finite Field of size 13
              defined by y^2 = x^6 + x + 1
        """
        if R in Fields():
            K = self._parent.value_ring()
            phi = R.coerce_map_from(K)
        else:
            phi = R
            R = phi.codomain()
        J = self._parent.codomain()(R)
        coords = [c.map_coefficients(phi, R) if isinstance(c, Polynomial) else c for c in self]
        return J(*coords)

    def _pairing_extensions(self, n):
        r"""
        Iterate over embeddings of the field of definition `K` of ``self``
        into extensions `L / K` which can be used to compute pairings on
        the curve when `K` is too small to find suitable auxiliary divisors.

        The degrees `[L : K]` run over the odd primes not dividing `n`:
        odd-degree extensions preserve inertness, and degrees coprime to
        `n` permit recovering the Tate pairing over `K`.

        EXAMPLES::

            sage: R.<x> = GF(13)[]
            sage: H = HyperellipticCurve(x^5 + 2*x + 1)
            sage: P = H.jacobian()(x + 5, 5)
            sage: it = P._pairing_extensions(15)
            sage: [next(it).codomain() for _ in range(2)]
            [Finite Field in z7 of size 13^7, Finite Field in z11 of size 13^11]
        """
        K = self._parent.value_ring()
        d = ZZ(2)
        while True:
            d = d.next_prime()
            if n % d:
                yield K.extension(d, map=True)[1]

    def _weil_pairing_value(self, Q, n, tries=5):
        r"""
        Return the Weil pairing `e_n(P, Q)` of `P` = ``self`` and ``Q``,
        or ``None`` if no suitable auxiliary divisors were found within
        ``tries`` attempts.

        See :meth:`weil_pairing`.

        EXAMPLES::

            sage: R.<x> = GF(2)[]
            sage: H = HyperellipticCurve(x^5 + 1, x^2 + 1)
            sage: P = H.jacobian()(x + 1, 0)
            sage: P._weil_pairing_value(P, 2) is None
            True
        """
        P = self
        for _ in range(tries):
            try:
                PS, S = P._random_auxiliary_divisor()
                QR, R = Q._random_auxiliary_divisor(PS)
                if not S._u.gcd(QR._u * R._u).is_one():
                    continue
                # f_P = f_{n,P+S} / f_{n,S} has divisor n (E_{P+S} - E_S)
                a1, a2 = PS._miller_values(n, [QR.uv(), R.uv()])
                b1, b2 = S._miller_values(n, [QR.uv(), R.uv()])
                # f_Q = f_{n,Q+R} / f_{n,R} has divisor n (E_{Q+R} - E_R)
                c1, c2 = QR._miller_values(n, [PS.uv(), S.uv()])
                d1, d2 = R._miller_values(n, [PS.uv(), S.uv()])
            except (ValueError, ZeroDivisionError):
                continue
            return (a1 * b2 * c2 * d1) / (a2 * b1 * c1 * d2)
        return None

    def _tate_miller_value(self, Q, n, tries=5):
        r"""
        Return `f_P(D_Q)`, where `P` = ``self``, `f_P` is a function with
        divisor `n D_P` and `D_Q` is a random degree-zero divisor over the
        field of definition `K` representing ``Q``, or ``None`` if no
        suitable `D_Q` was found in ``tries`` attempts.

        The output is well-defined modulo `n`-th powers in `K^*`. See
        :meth:`tate_pairing`.

        EXAMPLES::

            sage: R.<x> = GF(157)[]
            sage: H = HyperellipticCurve(x^5 + 66*x^4 + 47*x^3 + 15*x^2 + 65*x + 7)
            sage: J = H.jacobian()
            sage: P = J(x^2 + 45*x + 15, 75*x + 136)
            sage: Q = J(x^2 + 84*x + 87, 81*x + 25)
            sage: P._tate_miller_value(Q, 13)^12  # random
            101
            sage: P._tate_miller_value(Q, 13)^12 == P.tate_pairing(Q, 13)
            True
        """
        for _ in range(tries):
            try:
                QR, R = Q._random_auxiliary_divisor(self)
                a1, a2 = self._miller_values(n, [QR.uv(), R.uv()])
            except (ValueError, ZeroDivisionError):
                continue
            return a1 / a2
        return None

    def _check_pairing_input(self, other, n):
        r"""
        Check the input of :meth:`weil_pairing` and :meth:`tate_pairing`.

        TESTS::

            sage: R.<x> = GF(13)[]
            sage: H = HyperellipticCurve(x^5 + 2*x + 1)
            sage: J = H.jacobian()
            sage: P = J(x + 5, 5)
            sage: P._check_pairing_input(P, 0)
            Traceback (most recent call last):
            ...
            ValueError: n must be a positive integer
            sage: P._check_pairing_input(P, 2)
            Traceback (most recent call last):
            ...
            ValueError: the element P = (x + 5, 5) must be 2-torsion
            sage: P._check_pairing_input(J.change_ring(GF(13^2))(x + 5, 5), 2)
            Traceback (most recent call last):
            ...
            ValueError: elements must lie in the same group of points of the Jacobian
        """
        n = ZZ(n)
        if n <= 0:
            raise ValueError("n must be a positive integer")
        if other.parent() is not self.parent():
            raise ValueError("elements must lie in the same group of points of the Jacobian")
        if not isinstance(self._parent.value_ring(), FiniteField_generic):
            raise NotImplementedError("pairings are only implemented for Jacobians over finite fields")
        if n * self:
            raise ValueError(f"the element P = {self} must be {n}-torsion")
        return n

    def weil_pairing(self, Q, n):
        r"""
        Compute the Weil pairing `e_n(P, Q)` of `P` = ``self`` and ``Q``.

        For `n`-torsion elements `P`, `Q` of the Jacobian, choose degree-zero
        divisors `D_P`, `D_Q` with disjoint supports, and functions `f_P`, `f_Q`
        with divisors `n D_P` and `n D_Q`. The Weil pairing then equals

        .. MATH::

            e_n(P, Q) = \frac{f_P(D_Q)}{f_Q(D_P)}.

        For elliptic curves this agrees with the convention of
        :meth:`sage.schemes.elliptic_curves.ell_point.EllipticCurvePoint_field.weil_pairing`.

        The functions are computed using Miller's algorithm, where each
        step of Cantor's algorithm contributes a function on the curve.
        The divisors `D_P` and `D_Q` are chosen at random from their class.
        If no suitable divisors are found after a few attempts, which may
        happen over very small fields, the computation is carried out over
        an extension field instead.

        INPUT:

        - ``Q`` -- an element of the same group of points as ``self``

        - ``n`` -- positive integer such that `nP = nQ = 0`

        OUTPUT: an `n`-th root of unity in the base field

        EXAMPLES:

        The Jacobian of `y^2 = x^5 + x` over `\FF_7` is supersingular;
        over `\FF_{7^2}` the Frobenius acts as `-7`, so the full
        `8`-torsion is rational there::

            sage: R.<x> = GF(7)[]
            sage: H = HyperellipticCurve(x^5 + x)
            sage: K.<z2> = GF(7^2)
            sage: J = H.jacobian()(K)
            sage: J.order().factor()
            2^12
            sage: P = J(x^2 + (3*z2 + 1)*x + 2*z2 + 2, (3*z2 + 6)*x + z2)
            sage: Q = J(x^2 + (3*z2 + 5)*x + 3*z2 + 3, (5*z2 + 1)*x + 2*z2)
            sage: e = P.weil_pairing(Q, 8); e
            5*z2 + 3
            sage: e.multiplicative_order()
            8

        The Weil pairing is bilinear, alternating and Galois equivariant::

            sage: P.weil_pairing(P, 8)
            1
            sage: Q.weil_pairing(P, 8) == e^-1
            True
            sage: (3*P).weil_pairing(5*Q, 8) == e^15
            True
            sage: (P + Q).weil_pairing(Q, 8) == e
            True
            sage: frob = lambda D: J(*[c.map_coefficients(lambda a: a^7) for c in D])
            sage: frob(P).weil_pairing(frob(Q), 8) == e^7
            True

        A split example of genus 3::

            sage: R.<x> = GF(23)[]
            sage: H = HyperellipticCurve(x^8 + 1)
            sage: H.is_split()
            True
            sage: K.<z2> = GF(23^2)
            sage: J = H.jacobian()(K)
            sage: J.order().factor()
            2^18 * 3^6
            sage: P = J(x^3 + (5*z2 + 4)*x^2 + (5*z2 + 8)*x + 19*z2 + 13,
            ....:       (19*z2 + 12)*x^2 + (6*z2 + 15)*x + 17*z2 + 6)
            sage: Q = J(x^3 + (16*z2 + 9)*x^2 + 15*z2*x + 4*z2 + 7,
            ....:       5*x^2 + (8*z2 + 12)*x + z2 + 17)
            sage: P.weil_pairing(Q, 3)
            19*z2 + 15

        A larger example: the genus-3 curve `y^2 = x^7 + x` is superspecial
        modulo `p`, and the `(p+1)`-torsion is rational over `\FF_{p^2}`.
        We pair elements of large prime order `\ell` dividing `p+1`, as well
        as elements of order `p+1`::

            sage: p = 10000000000000001119
            sage: R.<x> = GF(p)[]
            sage: H = HyperellipticCurve(x^7 + x)
            sage: K.<z2> = GF(p^2)
            sage: J = H.jacobian()(K)
            sage: (p + 1) * J.random_element()
            (1, 0)
            sage: factor(p + 1)
            2^5 * 5 * 197 * 317258883248731
            sage: l = 317258883248731
            sage: P = ((p + 1) // l) * J.random_element()
            sage: Q = ((p + 1) // l) * J.random_element()
            sage: e = P.weil_pairing(Q, l)
            sage: e.multiplicative_order() == l  # random (False with probability about 1/l)
            True
            sage: (2*P).weil_pairing(3*Q, l) == e^6
            True
            sage: Q.weil_pairing(P, l) == e^-1
            True
            sage: frob = lambda D: J(*[c.map_coefficients(lambda a: a^p) for c in D])
            sage: frob(P).weil_pairing(frob(Q), l) == e^p
            True
            sage: P.tate_pairing(Q, l) / Q.tate_pairing(P, l) == e^((p^2-1) // l)
            True
            sage: P = J.random_element()
            sage: Q = J.random_element()
            sage: e = P.weil_pairing(Q, p + 1)
            sage: e^(p + 1)
            1
            sage: (P + Q).weil_pairing(Q, p + 1) == e
            True

        For elliptic curves, the result matches the Weil pairing
        implementation in
        :meth:`sage.schemes.elliptic_curves.ell_point.EllipticCurvePoint_field.weil_pairing`::

            sage: K.<z2> = GF(11^2)
            sage: E = EllipticCurve(K, [0, 1])
            sage: R.<x> = GF(11)[]
            sage: J = HyperellipticCurve(x^3 + 1).jacobian()(K)
            sage: P, Q = E.torsion_basis(12)
            sage: to_J = lambda T: J(x - T.x(), T.y())
            sage: to_J(P).weil_pairing(to_J(Q), 12) == P.weil_pairing(Q, 12)
            True

        Pairings can also be computed over very small fields, where the
        computation may be carried out over an extension field. On the
        curve `y^2 = x^5 - x` over `\FF_5`, the `2`-torsion is
        rational and the classes `T_a = (a, 0) - \infty` satisfy
        `e_2(T_a, T_b) = -1` for `a \neq b`::

            sage: S.<t> = GF(5)[]
            sage: J5 = HyperellipticCurve(t^5 - t).jacobian()
            sage: T = [J5(t - a, 0) for a in GF(5)]
            sage: matrix([[A.weil_pairing(B, 2) for B in T] for A in T])
            [1 4 4 4 4]
            [4 1 4 4 4]
            [4 4 1 4 4]
            [4 4 4 1 4]
            [4 4 4 4 1]

        TESTS:

        Both elements must be `n`-torsion::

            sage: to_J(P).weil_pairing(to_J(Q), 6)
            Traceback (most recent call last):
            ...
            ValueError: the element P = ... must be 6-torsion
            sage: to_J(6*P).weil_pairing(to_J(Q), 6)
            Traceback (most recent call last):
            ...
            ValueError: the element Q = ... must be 6-torsion

        Very small fields, where there may be no suitable auxiliary divisors
        over the base field::

            sage: R.<x> = GF(2)[]
            sage: J = HyperellipticCurve(x^5 + 1, x^2 + 1).jacobian()
            sage: J.points()
            [(1, 0), (x + 1, 0)]
            sage: J(x + 1, 0).weil_pairing(J(x + 1, 0), 2)
            1
            sage: K.<a> = GF(4)
            sage: R.<x> = K[]
            sage: H = HyperellipticCurve(x^6 + a*x^4 + a*x^3 + a*x + a, a*x^3 + 1)
            sage: H.is_inert(), H.rational_points()
            (True, [])
            sage: J = H.jacobian()(K)
            sage: P = J(x^2 + a*x + a, 0); 5*P
            (1, 0 : 1)
            sage: P.weil_pairing(2*P, 5)
            1
        """
        P = self
        n = P._check_pairing_input(Q, n)
        if n * Q:
            raise ValueError(f"the element Q = {Q} must be {n}-torsion")

        e = P._weil_pairing_value(Q, n)
        if e is not None:
            return e

        # The base field is too small: compute the pairing over an extension.
        # The value is the same, since the Weil pairing is intrinsic.
        for phi in P._pairing_extensions(n):
            e = P.change_ring(phi)._weil_pairing_value(Q.change_ring(phi), n)
            if e is not None:
                return phi.section()(e)

    def tate_pairing(self, Q, n, k=None, q=None):
        r"""
        Return the (reduced) Tate pairing of `P` = ``self`` and ``Q``.

        Let `P` be an `n`-torsion element of the Jacobian over
        `K = \FF_{q^k}`, where `n` divides `q^k-1`, and let `Q` be any
        element of the Jacobian over `K`. Choose a function `f_P` with
        divisor `n D_P`, and a divisor `D_Q` with support disjoint
        from `D_P`. The reduced Tate pairing then equals

        .. MATH::

            t_n(P, Q) = f_P(D_Q)^{(q^k-1)/n}.

        The Tate pairing is bilinear and nondegenerate as a pairing
        `J[n](K) \times J(K)/nJ(K) \to \mu_n`.

        INPUT:

        - ``Q`` -- an element of the same group of points of the Jacobian as
          ``self``

        - ``n`` -- positive integer such that `nP = 0`

        - ``k`` -- (optional) positive integer; the embedding degree

        - ``q`` -- (optional) positive integer; size of the base field
          (the "big" field is `\FF_{q^k}`). By default, `q^k` is the
          size of the field of definition of the elements.

        OUTPUT: an `n`-th root of unity in the base field

        ALGORITHM:

        The function `f_P` is computed using Miller's algorithm and evaluated
        at a random divisor `D_Q`. If no suitable divisor is found after a
        few attempts, which may happen over very small fields, we compute
        `u = f_P(D_Q')` for a divisor `D_Q'` over an extension `L / K` of
        degree `d` coprime to `n`. Then `N_{L/K}(u)` agrees with
        `f_P(D_Q)^d` up to `n`-th powers in `K^*`, so we recover
        `t_n(P, Q) = (N_{L/K}(u)^{(q^k-1)/n})^{1/d}`.

        EXAMPLES:

        Since `n` must divide `q^k-1`, the Tate pairing can be
        nondegenerate on elements defined over a field which does
        not contain the full `n`-torsion::

            sage: R.<x> = GF(157)[]
            sage: H = HyperellipticCurve(x^5 + 66*x^4 + 47*x^3 + 15*x^2 + 65*x + 7)
            sage: J = H.jacobian()
            sage: J.order().factor()
            2 * 7 * 13 * 137
            sage: P = J(x^2 + 45*x + 15, 75*x + 136)
            sage: Q = J(x^2 + 84*x + 87, 81*x + 25)
            sage: t = P.tate_pairing(Q, 13); t
            101
            sage: t.multiplicative_order()
            13
            sage: P.tate_pairing(2*Q, 13) == t^2
            True
            sage: P.tate_pairing(13*Q, 13)
            1

        The Weil pairing can be computed from the Tate pairing::

            sage: R.<x> = GF(7)[]
            sage: H = HyperellipticCurve(x^5 + x)
            sage: K.<z2> = GF(7^2)
            sage: J = H.jacobian()(K)
            sage: P = J(x^2 + (3*z2 + 1)*x + 2*z2 + 2, (3*z2 + 6)*x + z2)
            sage: Q = J(x^2 + (3*z2 + 5)*x + 3*z2 + 3, (5*z2 + 1)*x + 2*z2)
            sage: P.tate_pairing(Q, 8)
            2*z2 + 1
            sage: P.tate_pairing(Q, 8, k=2) == P.tate_pairing(Q, 8, q=7) == P.tate_pairing(Q, 8)
            True
            sage: P.tate_pairing(Q, 8) / Q.tate_pairing(P, 8) == P.weil_pairing(Q, 8)^6
            True

        For elliptic curves, we recover
        :meth:`sage.schemes.elliptic_curves.ell_point.EllipticCurvePoint_field.tate_pairing`::

            sage: K.<z2> = GF(11^2)
            sage: E = EllipticCurve(K, [0, 1])
            sage: R.<x> = GF(11)[]
            sage: J = HyperellipticCurve(x^3 + 1).jacobian()(K)
            sage: to_J = lambda T: J(x - T.x(), T.y())
            sage: P = E.torsion_basis(12)[0]
            sage: Q = E.random_point()
            sage: to_J(P).tate_pairing(to_J(Q), 12) == P.tate_pairing(Q, 12, 2)
            True

        Over very small fields, the computation may be carried out over an
        extension field::

            sage: S.<t> = GF(3)[]
            sage: J3 = HyperellipticCurve(t^5 + t^4 + 2*t^2 + t + 2).jacobian()
            sage: J3.order()
            8
            sage: P3 = J3(t + 1, 0)
            sage: Q3 = J3(t + 2, 2)
            sage: Q3.order()
            8
            sage: P3.tate_pairing(Q3, 2)
            2
            sage: P3.tate_pairing(2*Q3, 2)
            1

        TESTS::

            sage: to_J(P).tate_pairing(to_J(Q), 5)
            Traceback (most recent call last):
            ...
            ValueError: the element P = ... must be 5-torsion
            sage: to_J(P).tate_pairing(to_J(Q), 12, q=5)
            Traceback (most recent call last):
            ...
            ValueError: the size of the field of definition must be a power of q = 5
            sage: to_J(P).tate_pairing(to_J(Q), 12, k=1, q=5)
            Traceback (most recent call last):
            ...
            ValueError: n = 12 does not divide q^k-1 = 4
        """
        P = self
        n = P._check_pairing_input(Q, n)

        K = P.parent().value_ring()
        if q is None:
            if k is None:
                q, k = K.order(), 1
            elif K.degree().is_one():
                q = K.order()
            elif K.degree() == k:
                q = K.base_ring().order()
            else:
                raise ValueError("unexpected field degree: set keyword argument q equal to the "
                                 f"size of the base field (big field is GF(q^{k}))")
        elif k is None:
            try:
                k = ZZ(K.order().log(q))
            except TypeError:
                raise ValueError(f"the size of the field of definition must be a power of q = {q}")
        e, r = (q**k-1).quo_rem(n)
        if r:
            raise ValueError(f"n = {n} does not divide q^k-1 = {q**k-1}")

        u = P._tate_miller_value(Q, n)
        if u is not None:
            return u**e

        # The base field is too small: compute over an extension of degree coprime to n.
        for phi in P._pairing_extensions(n):
            u = P.change_ring(phi)._tate_miller_value(Q.change_ring(phi), n)
            if u is not None:
                L = phi.codomain()
                d = L.degree() // K.degree()
                u = phi.section()(u**((L.order()-1) // (K.order()-1)))
                return (u**e)**d.inverse_mod(n)


# =======================================================================
# Child classes for representation for ramified, inert and split models
# =======================================================================


class MumfordDivisorClassFieldRamified(MumfordDivisorClassField):
    def __init__(self, parent, u, v, check=True) -> None:
        r"""
        Create an element of the Jacobian of a ramified
        hyperelliptic curve.

        TESTS::

            sage: R.<x> = GF(5)[]
            sage: H = HyperellipticCurve(x^7 + 3*x + 2, x^2 + 1)
            sage: H.is_ramified()
            True
            sage: J = Jacobian(H)
            sage: D = J(x^3 + x^2 + 3*x + 1, 4*x^2 + 2*x + 3)
            sage: type(D)
            <class 'sage.schemes.hyperelliptic_curves.jacobian_morphism.MumfordDivisorClassFieldRamified'>
        """
        if not parent.curve().is_ramified():
            raise TypeError("hyperelliptic curve must be ramified")

        super().__init__(parent, u, v, check=check)


class MumfordDivisorClassFieldInert(MumfordDivisorClassField):
    def __init__(self, parent, u, v, check=True) -> None:
        r"""
        Create an element of the Jacobian of an inert
        hyperelliptic curve.

        TESTS::

            sage: R.<x> = GF(5)[]
            sage: H = HyperellipticCurve(2*x^6 + 1)
            sage: H.is_inert()
            True
            sage: J = Jacobian(H)
            sage: D = J(x^2 + x + 4, 4*x)
            sage: type(D)
            <class 'sage.schemes.hyperelliptic_curves.jacobian_morphism.MumfordDivisorClassFieldInert'>
        """
        if not parent.curve().is_inert():
            raise TypeError("hyperelliptic curve must be inert")

        if u.degree() % 2:
            raise ValueError(f"mumford coordinate {u} must have even degree")

        g = parent.curve().genus()
        self._n = (g - u.degree()) // 2
        super().__init__(parent, u, v, check=check)

    def _repr_(self) -> str:
        """
        Return a representation of the element.

        TESTS::

            sage: R.<x> = GF(5)[]
            sage: H = HyperellipticCurve(2*x^6 + 1)
            sage: H.is_inert()
            True
            sage: J = Jacobian(H)
            sage: D = J(x^2 + x + 4, 4*x); D # indirect doctest
            (x^2 + x + 4, 4*x : 0)
        """
        return f"({self._u}, {self._v} : {self._n})"

    def _latex_(self) -> str:
        """
        Return a LaTeX representation of the element.

        TESTS::

            sage: R.<x> = GF(5)[]
            sage: H = HyperellipticCurve(2*x^6 + 1)
            sage: H.is_inert()
            True
            sage: J = Jacobian(H)
            sage: D = J(x^2 + x + 4, 4*x)
            sage: latex(D)
            (x^{2} + x + 4, 4 x : 0)
        """
        return f"({self._u._latex_()}, {self._v._latex_()} : {self._n})"

    def _add_(self, *args, **kwds):
        r"""
        Thin wrapper around :meth:`MumfordDivisorClassField._add_` which
        fails in case of odd genus since that case is not implemented.

        EXAMPLES::

            sage: R.<x> = GF(5)[]
            sage: H = HyperellipticCurve(2*x^6 + 1)
            sage: H.genus()
            2
            sage: H.is_inert()
            True
            sage: J = Jacobian(H)
            sage: D = J(x^2 + x + 4, 4*x)
            sage: D + D
            (x^2 + x + 4, x : 0)

        ::

            sage: H = HyperellipticCurve(2*x^8 + 1)
            sage: H.genus()
            3
            sage: H.is_inert()
            True
            sage: J = Jacobian(H)
            sage: D = J(x^2 + 3, 3*x)
            sage: D + D
            Traceback (most recent call last):
            ...
            NotImplementedError: unable to perform arithmetic for inert models of odd genus; consider extending the base field to adjoin the points at infinity
        """
        if self.parent().curve().genus() % 2:
            raise NotImplementedError('unable to perform arithmetic for inert models of odd genus; '
                                      'consider extending the base field to adjoin the points at infinity')
        return super()._add_(*args, **kwds)

    def _neg_(self, *args, **kwds):
        r"""
        Thin wrapper around :meth:`MumfordDivisorClassField._neg_` which
        fails in case of odd genus since that case is not implemented.

        EXAMPLES::

            sage: R.<x> = GF(5)[]
            sage: H = HyperellipticCurve(2*x^6 + 1)
            sage: H.genus()
            2
            sage: H.is_inert()
            True
            sage: J = Jacobian(H)
            sage: D = J(x^2 + x + 4, 4*x)
            sage: -D
            (x^2 + x + 4, x : 0)

        ::

            sage: H = HyperellipticCurve(2*x^8 + 1)
            sage: H.genus()
            3
            sage: H.is_inert()
            True
            sage: J = Jacobian(H)
            sage: D = J(x^2 + 3, 3*x)
            sage: -D
            Traceback (most recent call last):
            ...
            NotImplementedError: unable to perform arithmetic for inert models of odd genus; consider extending the base field to adjoin the points at infinity
        """
        if self.parent().curve().genus() % 2:
            raise NotImplementedError('unable to perform arithmetic for inert models of odd genus; '
                                      'consider extending the base field to adjoin the points at infinity')
        return super()._neg_(*args, **kwds)


class MumfordDivisorClassFieldSplit(MumfordDivisorClassField):
    def __init__(self, parent, u, v, n=0, check=True) -> None:
        r"""
        Create an element of the Jacobian of a split
        hyperelliptic curve.

        TESTS::

            sage: R.<x> = GF(5)[]
            sage: H = HyperellipticCurve(x^6 + 1)
            sage: H.is_split()
            True
            sage: J = Jacobian(H)
            sage: D = J(x^2 + 3, 3)
            sage: type(D)
            <class 'sage.schemes.hyperelliptic_curves.jacobian_morphism.MumfordDivisorClassFieldSplit'>
        """
        if not parent.curve().is_split():
            raise TypeError("hyperelliptic curve must be split")

        # Ensure the weight is set correctly
        g = parent.curve().genus()
        assert 0 <= n <= (g - u.degree())
        self._n = n
        self._m = g - u.degree() - n

        super().__init__(parent, u, v, check=check)

    def _repr_(self) -> str:
        r"""
        Return a representation of the element.

        TESTS::

            sage: R.<x> = GF(5)[]
            sage: H = HyperellipticCurve(x^6 + 1)
            sage: J = Jacobian(H)
            sage: D = J(x^2 + 3, 3); D
            (x^2 + 3, 3 : 0)
        """
        return f"({self._u}, {self._v} : {self._n})"

    def _latex_(self) -> str:
        r"""
        Return a LaTeX representation of the element.

        TESTS::

            sage: R.<x> = GF(5)[]
            sage: H = HyperellipticCurve(x^6 + 1)
            sage: J = Jacobian(H)
            sage: D = J(x^2 + 3, 3)
            sage: latex(D)
            (x^{2} + 3, 3 : 0)
        """
        return f"({self._u._latex_()}, {self._v._latex_()} : {self._n})"

    def __bool__(self) -> bool:
        r"""
        Return ``True`` if this is not the zero element of the Jacobian.

        TESTS:

        We consider a genus-3 hyperelliptic curve with two points
        at infinity ``P_0``, ``P_1``. Elements of the Jacobian are
        represented as ``[D - 2P_0 - P_1]``.
        In particular the zero element has Mumford presentation
        `(1, 0 : 2)`::

            sage: R.<x> = QQ[]
            sage: H = HyperellipticCurve(x^8  + 1)
            sage: J = Jacobian(H)
            sage: D = J(1,0,2); D
            (1, 0 : 2)
            sage: D.is_zero()
            True

        There are more elements that are only supported at infinity,
        but are non-zero::

            sage: [P0,P1] = H.points_at_infinity()
            sage: D = J(P0,P1); D
            (1, 0 : 3)
            sage: D.is_zero()
            False
        """
        g = self._parent.curve().genus()
        return not (self._u.is_one() and self._v.is_zero() and self._n == (g + 1) // 2)

    def __iter__(self) -> Iterator[Polynomial | int | Integer]:
        r"""
        TESTS:

        Indirect tests::

            sage: R.<x> = GF(7)[]
            sage: H = HyperellipticCurve(x^6 - x + 1)
            sage: J = H.jacobian()
            sage: P = J(x + 2, 5)
            sage: list(P)
            [x + 2, 5, 1]
            sage: tuple(P)
            (x + 2, 5, 1)
            sage: Q = J(x + 2, 5, 0)
            sage: list(Q)
            [x + 2, 5, 0]
            sage: tuple(Q)
            (x + 2, 5, 0)
            sage: P == Q
            False
        """
        yield from [self._u, self._v, self._n]

    def __hash__(self) -> int:
        r"""
        Compute the hash value of this element.

        TESTS::

            sage: R.<x> = GF(23)[]
            sage: H = HyperellipticCurve(x^6 + 1)
            sage: J = H.jacobian()
            sage: D1 = J(x + 6, 17, 0)
            sage: D2 = J(x + 6, 17, 1)
            sage: hash(D1) == hash(D2)
            False
        """
        return hash(tuple(self))

    def _add_(self, other) -> Self:
        r"""
        Return the sum of the two elements.

        Follows algorithm 3.7 of [Mireles2008]_.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: H = HyperellipticCurve(x^5 - x^4 + x^2 - x, x^3 + 1)
            sage: J = Jacobian(H)
            sage: D1 = J(x^2 + x, 0)
            sage: D2 = J(x^2, -x)
            sage: D1 + D2
            (x^2 - 3/2*x + 1/2, -5/2*x + 1/2 : 0)
            sage: D3 = J(x - 1, -2)
            sage: D1 + D3
            (x^2 - 1/2*x, 5/4*x - 1 : 0)
        """
        # Collect data from HyperellipticCurve
        H = self.parent().curve()
        g = H.genus()

        # Extract out mumford coordinates
        u1, v1 = self.uv()
        u2, v2 = other.uv()

        # Extract out integers for weights
        n1, n2 = self._n, other._n

        # Step one: cantor composition of the two divisors
        u3, v3, n3 = self._parent.cantor_composition(u1, v1, n1, u2, v2, n2)

        # Step two: cantor reduction of the above to ensure
        # the degree of u is smaller than g + 1
        while u3.degree() > (g + 1):
            u3, v3, n3 = self._parent.cantor_reduction(u3, v3, n3)

        # Step three: compose and then reduce at infinity to ensure
        # unique representation of D
        while n3 < 0 or n3 > g - u3.degree():
            plus = n3 >= 0
            u3, v3, n3, _ = self._parent.cantor_compose_at_infinity(u3, v3, n3, plus=plus)

        return self._parent(u3, v3, n3, check=False)

    def _miller_add(self, other) -> tuple:
        r"""
        Return the sum ``self + other`` together with a function `h`
        on the curve such that

        .. MATH::

            \mathrm{div}(h) = D_{\mathrm{self}} + D_{\mathrm{other}} - D_{\mathrm{self} + \mathrm{other}}`.

        See :meth:`MumfordDivisorClassField._miller_add` for the output format.

        Here `D_A = E_A + n \infty_+ + m \infty_- - D_\infty` is the
        balanced divisor of `A`, so that the functions arising from
        the steps of :meth:`_add_` also take into account the points
        at infinity.

        EXAMPLES::

            sage: R.<x> = GF(13)[]
            sage: H = HyperellipticCurve(x^6 + x + 1)
            sage: J = H.jacobian()
            sage: D1 = J(x + 3, 8, 0)
            sage: D2 = J(x^2 + 7, 7*x + 1, 0)
            sage: D3, h = D1._miller_add(D2)
            sage: D3 == D1 + D2
            True

        Here, the composition does not cancel any points, and the
        divisor is then reduced by a composition at infinity::

            sage: h
            [(1, 0, 1), (x^3 + 11*x^2 + 11, 1, x^2 + 9*x + 6)]
        """
        parent = self._parent
        H = parent.extended_curve()
        g = H.genus()
        R = self._u.parent()

        u1, v1 = self.uv()
        u2, v2 = other.uv()
        n1, n2 = self._n, other._n

        # Composition: D1 + D2 = D3 + div(s)
        u3, v3, s = parent._cantor_composition_generic(u1, v1, u2, v2)
        n3 = n1 + n2 + s.degree() - ((g + 1) // 2)
        funs = [(s, R.zero(), R.one())]

        # Reduction: D3 - D4 = div((y - v3) / u4)
        while u3.degree() > (g + 1):
            u4, v4, n4 = parent.cantor_reduction(u3, v3, n3)
            funs.append((-v3, R.one(), u4))
            u3, v3, n3 = u4, v4, n4

        # Composition at infinity: D3 - D4 = div((y - w) / u4)
        while n3 < 0 or n3 > g - u3.degree():
            plus = n3 >= 0
            u4, v4, n4, w = parent.cantor_compose_at_infinity(u3, v3, n3, plus=plus)
            funs.append((-w, R.one(), u4))
            u3, v3, n3 = u4, v4, n4

        return parent(u3, v3, n3, check=False), funs

    def _neg_(self):
        r"""
        Return the divisor multiplied by -1.

        Follows algorithm 3.8 of [Mireles2008]_.

        TESTS::

            sage: R.<x> = GF(101)[]
            sage: H = HyperellipticCurve(x^5 + 1, x^3 + x)
            sage: assert H.is_split()
            sage: J = Jacobian(H)
            sage: [P0,P1] = H.points_at_infinity()
            sage: P = H.random_point()
            sage: Q = H.hyperelliptic_involution(P)
            sage: - J(P,P0) == J(Q,P1)
            True
        """
        # Collect data from HyperellipticCurve
        H = self.parent().curve()
        _, h = H.hyperelliptic_polynomials()
        g = H.genus()

        u0, v0 = self.uv()
        n0, m0 = self._n, self._m

        # Case for even genus
        if not (g % 2):
            v1 = (-h - v0) % u0
            u1 = u0
            n1 = m0
        # Odd genus, positive n0
        elif n0 > 0:
            v1 = (-h - v0) % u0
            u1 = u0
            # There is a typo for n1 in [GHM2008]_ and [Mireles2008]_.
            # This is corrected on page 617 of [Gal2018]_.
            n1 = m0 + 1
        else:
            # Composition at infinity always with plus=True.
            # want to "subtract" \infty_+ - \infty_-
            u1, v1, n1, _ = self._parent.cantor_compose_at_infinity(u0, -h - v0, n0, plus=True)
            n1 = n1 - n0 + m0 + 1  # Shouldn't this always be 0? (See Alg 3.8 of [Mireles2008]_)
            assert n1 == 0

        return self._parent(u1, v1, n1, check=False)
