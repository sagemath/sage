r"""
Rational point sets on a Jacobian of a hyperelliptic curve (split case)

Uses the balanced divisors technique described in [Mireles2008]_, [GHM2008]_, and [Gal2018]_.

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

from sage.rings.integer import Integer
from sage.schemes.hyperelliptic_curves.jacobian_homset_generic import (
    HyperellipticJacobianHomset,
)
from sage.schemes.hyperelliptic_curves.jacobian_morphism import (
    MumfordDivisorClassFieldSplit,
)
from sage.structure.element import parent


class HyperellipticJacobianHomsetSplit(HyperellipticJacobianHomset):
    Element = MumfordDivisorClassFieldSplit

    # Split models accept an explicit weight at infinity as a third argument.
    _max_constructor_args = 3

    def __init__(self, Y, X, **kwds):
        r"""
        Create the Jacobian Hom-set of a hyperelliptic curve with
        two rational points at infinity.

        TESTS::

            sage: R.<x> = GF(7)[]
            sage: H = HyperellipticCurve(x^6 + 2*x^2 + 1)
            sage: assert H.is_split()
            sage: JK = Jacobian(H)(GF(7))
            sage: type(JK)
            <class 'sage.schemes.hyperelliptic_curves.jacobian_g2_homset_split.HyperellipticJacobianHomsetSplit_g2_with_category'>
            sage: TestSuite(JK).run(skip='_test_elements')
        """
        super().__init__(Y, X, **kwds)
        self._morphism_element = MumfordDivisorClassFieldSplit

    def zero(self, check=True):
        r"""
        Return the zero element of the Jacobian

        EXAMPLES ::

            sage: R.<x> = QQ[]
            sage: H = HyperellipticCurve(x^8 + 1)
            sage: J = Jacobian(H)
            sage: J.zero()
            (1, 0 : 2)
        """
        g = self.extended_curve().genus()
        R = self.extended_curve().polynomial_ring()
        n = (g + 1) // 2
        return self._morphism_element(self, R.one(), R.zero(), n)

    def point_to_mumford_coordinates(self, P):
        r"""
        On input a point ``P``, return the Mumford coordinates
        of (the affine part of) the divisor `[P]` and an integer `n`,
        where

        - `n = 1` if ``P`` is the point `\infty_+`
        - `n = 0` otherwise .

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: H = HyperellipticCurve(x^6 - 8*x^4 + 6*x^3 + 8*x^2 - 4*x + 1)
            sage: P = H([-1, 46, 3]); P
            (-1/3 : 46/27 : 1)
            sage: O = H([1,1,0])
            sage: JQ = H.jacobian()(QQ)
            sage: JQ.point_to_mumford_coordinates(P)
            (x + 1/3, 46/27, 0)
            sage: JQ.point_to_mumford_coordinates(O)
            (1, 0, 1)
        """

        R, x = self.extended_curve().polynomial_ring().objgen()
        [X, Y, Z] = P._coords
        if Z == 0:
            alpha = Y / X
            if alpha == self.extended_curve().roots_at_infinity()[0]:
                n = 1
            else:
                n = 0
            return R.one(), R.zero(), n
        u = x - X
        v = R(Y)
        return u, v, 0

    def _cantor_compose_points(self, P1, P2):
        r"""
        Return the Mumford coordinates ``(u, v, n)`` of the class `[P1 - P2]`.

        Here ``n`` is the multiplicity of the point `\infty_+` in the balanced
        representation. This is the split-model hook for
        :meth:`~sage.schemes.hyperelliptic_curves.jacobian_homset_generic.HyperellipticJacobianHomset._element_constructor_`.

        We use `P_2 + \iota(P_2) \sim \infty_+ + \infty_-` to write
        `[P_1 - P_2] = [P_1] + [\iota(P_2)] - \infty_+ - \infty_-`. Composing the
        affine parts of `[P_1]` and `[\iota(P_2)]` removes ``s_deg`` pairs of
        mutually inverse points, each contributing `\infty_+ + \infty_-`, so the
        multiplicity of `\infty_+` is increased by ``s_deg``.

        EXAMPLES::

            sage: R.<x> = GF(13)[]
            sage: H = HyperellipticCurve(x^8 + x + 1)
            sage: H.is_split()
            True
            sage: J = Jacobian(H)
            sage: JH = J.point_homset()
            sage: P = H.lift_x(1)
            sage: Q = H.lift_x(2)
            sage: D1 = JH(P); D1
            (x + 12, 4 : 1)
            sage: D2 = JH(Q); D2
            (x + 11, 5 : 1)
            sage: D = JH(P, Q); D
            (x^2 + 10*x + 2, 4*x : 1)
            sage: D == D1 - D2
            True
            sage: JH(x^2 + 10*x + 2, 4*x, 1) == D
            True
            sage: JH(x^2 + 10*x + 2, 4*x, 0) == D
            False

        The points at infinity may also be embedded into the Jacobian::

            sage: [P0, P1] = H.points_at_infinity()
            sage: JH(P0)
            (1, 0 : 2)
            sage: JH(P1)
            (1, 0 : 1)

        For any point ``P``, the class ``[P - P]`` must be the identity
        of the Jacobian (see :issue:`42349`)::

            sage: R.<x> = GF(101)[]
            sage: H = HyperellipticCurve(x^8 + x + 1)
            sage: J = Jacobian(H)
            sage: P = H.lift_x(3)
            sage: J(P, P)
            (1, 0 : 2)
            sage: J(P, P) == J.zero()
            True
            sage: J(P, P).is_zero()
            True

        This also holds on genus 1 split models, where the affine part of
        ``[P - Q]`` may exceed the genus before reduction::

            sage: R.<x> = GF(101)[]
            sage: H = HyperellipticCurve(x^4 + x + 1)
            sage: J = Jacobian(H)
            sage: J(H.lift_x(2), H.lift_x(3)) == J(H.lift_x(2)) - J(H.lift_x(3))
            True
            sage: J(H.lift_x(2), H.lift_x(2)).is_zero()
            True

        Construction also works over an extension `L` of the base field, using
        the curve base-extended to `L`::

            sage: R.<x> = GF(13)[]
            sage: H = HyperellipticCurve(x^8 + x + 1)
            sage: L.<z2> = GF(13^2)
            sage: JL = Jacobian(H)(L)
            sage: HL = JL.extended_curve()
            sage: P = HL.lift_x(z2 + 2); Q = HL.lift_x(2*z2 + 1)
            sage: JL(P)
            (x + 12*z2 + 11, 5*z2 + 11 : 1)
            sage: JL(P, Q)
            (x^2 + (10*z2 + 10)*x + 7*z2 + 11, (8*z2 + 4)*x + 3*z2 + 6 : 1)
            sage: JL(P, Q) == JL(P) - JL(Q)
            True
            sage: JL(P, P).is_zero()
            True
        """
        H = self.extended_curve()
        g = H.genus()
        u1, v1, n1 = self.point_to_mumford_coordinates(P1)
        P2_inv = H.hyperelliptic_involution(P2)
        u2, v2, n2 = self.point_to_mumford_coordinates(P2_inv)
        u, v, s_deg = self._cantor_composition_generic(u1, v1, u2, v2)
        n = (g + 1) // 2 - 1 + n1 + n2 + s_deg
        # Reduce to the canonical balanced representative. This is a no-op when
        # the divisor is already reduced (genus >= 2), and handles the genus 1
        # case where ``deg(u)`` can exceed ``g``.
        while u.degree() > g + 1:
            u, v, n = self.cantor_reduction(u, v, n)
        while n < 0 or n > g - u.degree():
            u, v, n = self.cantor_compose_at_infinity(u, v, n, plus=(n >= 0))
        return u, v, n

    def _mumford_from_coordinates(self, args):
        r"""
        Return the Mumford coordinates ``(u, v, n)`` parsed from ``args``.

        This is the split-model hook for
        :meth:`~sage.schemes.hyperelliptic_curves.jacobian_homset_generic.HyperellipticJacobianHomset._element_constructor_`.
        On input two polynomials ``(u, v)``, the weight defaults to
        `\lceil (g - \deg u) / 2 \rceil`; an explicit weight may be supplied as
        a third argument.

        EXAMPLES::

            sage: R.<x> = GF(13)[]
            sage: H = HyperellipticCurve(x^8 + x + 1)
            sage: J = Jacobian(H)
            sage: JH = J.point_homset()
            sage: D = JH(x^2 + 10*x + 2, 4*x); D
            (x^2 + 10*x + 2, 4*x : 1)
            sage: JH(x^2 + 10*x + 2, 4*x, 1) == D
            True
            sage: JH(x^2 + 10*x + 2, 4*x, 0) == D
            False

        The various spellings of the zero element agree::

            sage: J() == J(0) == J(1, 0) == J.zero() == JH(0) == 0
            True
        """
        R = self.extended_curve().polynomial_ring()
        P1, P2 = args[0], args[1]
        if not (R.coerce_map_from(parent(P1)) and R.coerce_map_from(parent(P2))):
            raise ValueError(
                "the input must consist of one or two points, or Mumford coordinates"
            )
        u = R(P1)
        v = R(P2)
        if len(args) == 3 and isinstance(args[2], (int, Integer)):
            n = args[2]
        else:
            g = self.extended_curve().genus()
            n = (g - u.degree() + 1) // 2
        return u, v, n

    def cantor_composition(self, u1, v1, n1, u2, v2, n2):
        r"""
        Return the Cantor composition of the divisors represented by
        ``(u1, v1, n1)`` and ``(u2, v2, n2)``.
        Here ``n1`` and ``n2`` denote the multiplicity of the point
        `\infty_+`.

        Follows algorithm 3.4 of [Mireles2008]_.

        TODO: when h = 0 we can speed this up.

        EXAMPLES::

            sage: R.<x> = GF(7)[]
            sage: H = HyperellipticCurve(x^8 + 3*x + 2)
            sage: JF = Jacobian(H).point_homset()
            sage: D1 = [x^2 + 4*x + 3, 2*x + 2, 1]
            sage: assert JF(D1)
            sage: D2 = [x^3 + 6*x^2 + 6*x, 6*x^2 + 6*x + 3, 0]
            sage: assert JF(D2)
            sage: D3 = JF.cantor_composition(*D1, *D2); D3
            (x^5 + 3*x^4 + 5*x^3 + 4*x, 3*x^3 + 3*x^2 + 3*x + 3, -1)
        """
        # Collect data from HyperellipticCurve
        H = self.extended_curve()
        g = H.genus()

        # Cantor composition
        u3, v3, s_deg = self._cantor_composition_generic(u1, v1, u2, v2)

        # Compute new weight
        n3 = n1 + n2 + s_deg - ((g + 1) // 2)

        return u3, v3, n3

    def cantor_reduction(self, u0, v0, n0):
        r"""
        Compute the Cantor reduction of ``(u0,v0,n0)``,
        where ``(u0,v0)`` represent an affine semi-reduced divisor and
        ``n0`` is the multiplicity of the point `\infty_+`.

        Follows algorithm 3.5 of [Mireles2008]_.

        EXAMPLES::

            sage: R.<x> = GF(7)[]
            sage: H = HyperellipticCurve(x^8 + 3*x + 2)
            sage: JF = Jacobian(H).point_homset()
            sage: D1 = [x^2 + 4*x + 3, 2*x + 2, 1]
            sage: D2 = [x^3 + 6*x^2 + 6*x, 6*x^2 + 6*x + 3, 0]
            sage: D3 = JF.cantor_composition(*D1, *D2); D3
            (x^5 + 3*x^4 + 5*x^3 + 4*x, 3*x^3 + 3*x^2 + 3*x + 3, -1)
            sage: JF.cantor_reduction(*D3)
            (x^3 + 4*x^2 + 2*x + 5, 2*x^2 + 3*x + 5, 0)
        """
        # Collect data from HyperellipticCurve
        H = self.extended_curve()
        g = H.genus()

        # Perform regular cantor reduction
        u1, v1 = self._cantor_reduction_generic(u0, v0)

        # Compute the counter weights
        d0 = u0.degree()
        d1 = u1.degree()
        a_plus, a_minus = H.roots_at_infinity()

        if v0.degree() <= g + 1:
            leading_coefficient = v0[g + 1]  # check coefficient of x^(g+1)
            if leading_coefficient == a_plus:
                n1 = n0 + d0 - g - 1
            elif leading_coefficient == a_minus:
                n1 = n0 + g + 1 - d1
            else:
                n1 = n0 + (d0 - d1) // 2
        else:
            n1 = n0 + (d0 - d1) // 2
        return u1, v1, n1

    def cantor_compose_at_infinity(self, u0, v0, n0, plus=True):
        r"""
        Compute the composition of `(u_0,v_0,n_0)` with a divisor supported
        at `\infty_+` (default) or `\infty_-` , and apply a reduction step.

        Follows algorithm 3.6 of [Mireles2008]_.

        EXAMPLES::

            sage: R.<x> = GF(7)[]
            sage: H = HyperellipticCurve(x^8 + 3*x + 2)
            sage: JF = Jacobian(H).point_homset()
            sage: D1 = [x^2 + 4*x + 3, 2*x + 2, 1]

        Composing at `\infty_+` decreases the value of `n_0` ,
        while composing at `\infty_-` increases that value::

            sage: JF.cantor_compose_at_infinity(x^2 + 4*x + 3, 2*x + 2, 1)
            (x^2 + 3*x + 6, 5*x + 5, -1)
            sage: JF.cantor_compose_at_infinity(x^2 + 4*x + 3, 2*x + 2, 1, plus=False)
            (x^3 + 6*x^2 + x + 4, 5*x + 5, 2)
        """
        # Collect data from HyperellipticCurve
        H = self.extended_curve()
        f, h = H.hyperelliptic_polynomials()
        g = H.genus()

        # Pick either G_plus or G_minus for reduction
        G_plus, G_minus = H.split_G_plus_minus()
        if plus:
            G = G_plus
        else:
            G = G_minus

        v1_prime = G + ((v0 - G) % u0)
        u1 = (v1_prime**2 + h * v1_prime - f) // u0
        u1 = u1.monic()
        v1 = (-h - v1_prime) % u1

        # Compute the counter weights
        if plus:
            n1 = n0 + u0.degree() - g - 1
        else:
            n1 = n0 + g + 1 - u1.degree()

        return u1, v1, n1
