"""
Rational point sets on a Jacobian of a hyperelliptic curve (split case)
"""
from sage.rings.integer import Integer
from sage.schemes.hyperelliptic_curves_smooth_model.jacobian_homset_generic import (
    HyperellipticJacobianHomset,
)
from sage.schemes.hyperelliptic_curves_smooth_model.jacobian_morphism import (
    MumfordDivisorClassFieldSplit,
)

from sage.schemes.weighted_projective.weighted_projective_point import (
    SchemeMorphism_point_weighted_projective_ring,
)
from sage.structure.element import parent


class HyperellipticJacobianHomsetSplit(HyperellipticJacobianHomset):
    def __init__(self, Y, X, **kwds):
        """
        Create the Jacobian Hom-set of a hyperelliptic curve with
        two rational points at infinity.

        TESTS::

            sage: R.<x> = GF(7)[]
            sage: H = HyperellipticCurveSmoothModel(x^6 + 2*x^2 + 1)
            sage: assert H.is_split()
            sage: JK = Jacobian(H)(GF(7))
            sage: type(JK)
            <class 'sage.schemes.hyperelliptic_curves_smooth_model.jacobian_g2_homset_split.HyperellipticJacobianHomsetSplit_g2_with_category'>
        """
        super().__init__(Y, X, **kwds)
        self._morphism_element = MumfordDivisorClassFieldSplit

    def zero(self, check=True):
        """
        Return the zero element of the Jacobian

        EXAMPLES ::

            sage: R.<x> = QQ[]
            sage: H = HyperellipticCurveSmoothModel(x^8 + 1)
            sage: J = Jacobian(H)
            sage: J.zero()
            (1, 0 : 2)
        """
        g = self.curve().genus()
        R = self.curve().polynomial_ring()
        n = (g + 1) // 2
        return self._morphism_element(self, R.one(), R.zero(), n)

    def point_to_mumford_coordinates(self, P):
        """
        On input a point P, return the Mumford coordinates
        of (the affine part of) the divisor [P] and an integer n,
        where
        * n = 1 if P is the point oo+
        * n = 0 otherwise .

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: H = HyperellipticCurveSmoothModel(x^6 - 8*x^4 + 6*x^3 + 8*x^2 - 4*x + 1)
            sage: P = H([-1, 46, 3]); P
            (-1/3 : 46/27 : 1)
            sage: O = H([1,1,0])
            sage: JQ = H.jacobian()(QQ)
            sage: JQ.point_to_mumford_coordinates(P)
            (x + 1/3, 46/27, 0)
            sage: JQ.point_to_mumford_coordinates(O)
            (1, 0, 1)
        """

        R, x = self.curve().polynomial_ring().objgen()
        [X, Y, Z] = P._coords
        if Z == 0:
            alpha = Y / X
            if alpha == self.curve().roots_at_infinity()[0]:
                n = 1
            else:
                n = 0
            return R.one(), R.zero(), n
        u = x - X
        v = R(Y)
        return u, v, 0

    def __call__(self, *args, check=True):
        r"""
        Return a rational point in the abstract Homset `J(K)`, given:

        1. No arguments or the integer `0`; return `0 \in J`;

        2. A point `P` on `J = Jac(C)`, return `P`;
        3. Polynomials u,v such that `v^2 + h*v - f = 0 mod u`,
           returning J(u,v,n) with n = ((g - deg(u))/2).ceil().

        4. Polynomials u,v and an integer n such that `v^2 + h*v - f = 0 mod u`
           and 0 <= n <= g/2,
           returning J(u,v,n).

        Return a rational point in the abstract Homset `J(K)`, given:

        1. No arguments or the integer `0`; return `0 \in J`;

        2. A point `P` on `J = Jac(C)`, return `P`;

        3. A point `P` on the curve `H` such that `J = Jac(H)`;
           return `[P - P_0]`, where `P_0` is the distinguished point of `H`.
           By default, `P_0 = \infty`;

        4. Two points `P, Q` on the curve `H` such that `J = Jac(H)`;
           return `[P - Q]`;

        5. Polynomials `(u, v)` such that `v^2 + hv - f \equiv 0 \pmod u`;
           reutrn `[(u(x), y - v(x)) : \lceil (g - \deg(u)) / 2 \rceil]`;

        6. Polynomials `(u, v)` and an integer `n` such that
           `v^2 + hv - f \equiv 0 \pmod u` and `0 \leq n \leq g / 2`;
           return `[u, v : n]`.

        EXAMPLES::

            sage: R.<x> = GF(13)[]
            sage: H = HyperellipticCurveSmoothModel(x^8 + x + 1)
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
            sage: D = JH(P,Q); D
            (x^2 + 10*x + 2, 4*x : 1)
            sage: D == D1 - D2
            True
            sage: JH(x^2+10*x+2, 4*x, 1) == D
            True
            sage: JH(x^2+10*x+2, 4*x, 0) == D
            False

        The points at infinity may also be embedded into the Jacobian::

            sage: [P0, P1] = H.points_at_infinity()
            sage: JH(P0)
            (1, 0 : 2)
            sage: JH(P1)
            (1, 0 : 1)

        TESTS::

            sage: R.<x> = GF(13)[]
            sage: H = HyperellipticCurveSmoothModel(x^8 + x + 1)
            sage: J = Jacobian(H)
            sage: JH = J.point_homset()
            sage: J() == J(0) == J(1, 0) == J.zero() == JH(0) == 0
            True

        TODO: Merge this code with that of `HyperellipticJacobianHomset`
        """
        g = self.curve().genus()
        R = self.curve().polynomial_ring()

        if len(args) > 3:
            raise ValueError("at most three arguments are allowed as input")

        if len(args) == 0 or (len(args) == 1 and args[0] == ()):
            return self._morphism_element(
                self, R.one(), R.zero(), n=(g + 1) // 2, check=check
            )

        if len(args) == 1 and isinstance(args[0], (list, tuple)):
            args = args[0]

        if len(args) == 1:
            P1 = args[0]
            if P1 == 0:
                u = R.one()
                v = R.zero()
                n = (g + 1) // 2
            elif isinstance(P1, self._morphism_element):
                return P1
            elif isinstance(P1, SchemeMorphism_point_weighted_projective_ring):
                # TODO: Test this path when args is a tuple
                args = args + (
                    self.curve().distinguished_point(),
                )  # this case will now be handled below.
            else:
                raise ValueError(
                    "the input must consist of one or two points, or Mumford coordinates"
                )

        if len(args) == 2 or len(args) == 3:
            P1 = args[0]
            P2 = args[1]
            if isinstance(P1, SchemeMorphism_point_weighted_projective_ring) and isinstance(
                P2, SchemeMorphism_point_weighted_projective_ring
            ):
                if len(args) == 3:
                    raise ValueError("the input must consist of at most two points")
                u1, v1, n1 = self.point_to_mumford_coordinates(P1)
                P2_inv = self.curve().hyperelliptic_involution(P2)
                u2, v2, n2 = self.point_to_mumford_coordinates(P2_inv)
                u, v, _ = self._cantor_composition_generic(u1, v1, u2, v2)
                n = (g + 1) // 2 - 1 + n1 + n2  # this solution is a bit hacky
            # This checks whether P1 and P2 can be interpreted as polynomials
            elif R.coerce_map_from(parent(P1)) and R.coerce_map_from(parent(P2)):
                u = R(P1)
                v = R(P2)
                if len(args) == 3 and isinstance(args[2], (int, Integer)):
                    n = args[2]
                else:
                    n = (
                        g - u.degree() + 1
                    ) // 2  # TODO: do we really want to allow this input?
            else:
                raise ValueError(
                    "the input must consist of one or two points, or Mumford coordinates"
                )

        return self._morphism_element(self, u, v, n=n, check=check)

    def cantor_composition(self, u1, v1, n1, u2, v2, n2):
        r"""
        Return the Cantor composition of the divisors represented by
        ``(u1, v1, n1)`` and ``(u2, v2, n2)``.
        Here ``n1`` and ``n2`` denote the multiplicity of the point
        `\infty_+`.

        Follows algorithm 3.4 of

        Efficient Arithmetic on Hyperelliptic Curves With Real Representation
        David J. Mireles Morales (2008)
        https://www.math.auckland.ac.nz/~sgal018/Dave-Mireles-Full.pdf

        TODO: when h = 0 we can speed this up.

        EXAMPLES::

            sage: R.<x> = GF(7)[]
            sage: H = HyperellipticCurveSmoothModel(x^8 + 3*x + 2)
            sage: JF = Jacobian(H).point_homset()
            sage: D1 = [x^2 + 4*x + 3, 2*x + 2, 1]
            sage: assert JF(D1)
            sage: D2 = [x^3 + 6*x^2 + 6*x, 6*x^2 + 6*x + 3, 0]
            sage: assert JF(D2)
            sage: D3 = JF.cantor_composition(*D1, *D2); D3
            (x^5 + 3*x^4 + 5*x^3 + 4*x, 3*x^3 + 3*x^2 + 3*x + 3, -1)
        """
        # Collect data from HyperellipticCurve
        H = self.curve()
        g = H.genus()

        # Cantor composition
        u3, v3, s_deg = self._cantor_composition_generic(u1, v1, u2, v2)

        # Compute new weight
        n3 = n1 + n2 + s_deg - ((g + 1) // 2)

        return u3, v3, n3

    def cantor_reduction(self, u0, v0, n0):
        """
        Compute the Cantor reduction of ``(u0,v0,n0)``,
        where ``(u0,v0)`` represent an affine semi-reduced divisor and
        n0 is the multiplicity of the point infty+.

        Follows algorithm 3.5 of

        Efficient Arithmetic on Hyperelliptic Curves With Real Representation
        David J. Mireles Morales (2008)
        https://www.math.auckland.ac.nz/~sgal018/Dave-Mireles-Full.pdf

        EXAMPLES::

            sage: R.<x> = GF(7)[]
            sage: H = HyperellipticCurveSmoothModel(x^8 + 3*x + 2)
            sage: JF = Jacobian(H).point_homset()
            sage: D1 = [x^2 + 4*x + 3, 2*x + 2, 1]
            sage: D2 = [x^3 + 6*x^2 + 6*x, 6*x^2 + 6*x + 3, 0]
            sage: D3 = JF.cantor_composition(*D1, *D2); D3
            (x^5 + 3*x^4 + 5*x^3 + 4*x, 3*x^3 + 3*x^2 + 3*x + 3, -1)
            sage: JF.cantor_reduction(*D3)
            (6*x^3 + 3*x^2 + 5*x + 2, 2*x^2 + 3*x + 5, 0)
        """
        # Collect data from HyperellipticCurve
        H = self.curve()
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

        Follows algorithm 3.6 of

        Efficient Arithmetic on Hyperelliptic Curves With Real Representation
        David J. Mireles Morales (2008)
        https://www.math.auckland.ac.nz/~sgal018/Dave-Mireles-Full.pdf

        EXAMPLES::

            sage: R.<x> = GF(7)[]
            sage: H = HyperellipticCurveSmoothModel(x^8 + 3*x + 2)
            sage: JF = Jacobian(H).point_homset()
            sage: D1 = [x^2 + 4*x + 3, 2*x + 2, 1]

        Composing at `\infty_+` decreases the value of `n_0` ,
        while composing at `\infty_-` increases that value.

            sage: JF.cantor_compose_at_infinity(x^2 + 4*x + 3, 2*x + 2, 1)
            (x^2 + 3*x + 6, 5*x + 5, -1)
            sage: JF.cantor_compose_at_infinity(x^2 + 4*x + 3, 2*x + 2, 1, plus=False)
            (x^3 + 6*x^2 + x + 4, 5*x + 5, 2)
        """
        # Collect data from HyperellipticCurve
        H = self.curve()
        f, h = H.hyperelliptic_polynomials()
        g = H.genus()

        # Pick either G_plus or G_minus for reduction
        G_plus, G_minus = H.infinite_polynomials()
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
