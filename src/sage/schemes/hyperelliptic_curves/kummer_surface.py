"""
Kummer surfaces over a general ring

AUTHORS:

- David Kohel (2006): initial version
- Sabrina Kunzweiler, Gareth Ma, Giacomo Pope (2026): add map from the Jacobian
"""

# *****************************************************************************
#  Copyright (C) 2006 David Kohel <kohel@maths.usyd.edu>
#                2026 Sabrina Kunzweiler <sabrina.kunzweiler@math.u-bordeaux.fr>
#                2026 Gareth Ma <grhkm21@gmail.com>
#                2026 Giacomo Pope <giacomopope@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
# *****************************************************************************

from sage.schemes.hyperelliptic_curves.jacobian_morphism import MumfordDivisorClassField
from sage.schemes.projective.projective_space import ProjectiveSpace
from sage.schemes.projective.projective_subscheme import (
    AlgebraicScheme_subscheme_projective,
)


class KummerSurface(AlgebraicScheme_subscheme_projective):
    r"""
    Kummer surface of the Jacobian of a genus-2 curve.

    EXAMPLES::

        sage: R.<x> = GF(13)[]
        sage: C = HyperellipticCurve(x^6+2*x+1)
        sage: J = C.jacobian()
        sage: K = J.kummer_surface()

    Points can be constructed by providing explicit coordinates,
    or as the images of elements of the Jacobian::

        sage: P = K([3,6,0,1]); P
        (3 : 6 : 0 : 1)
        sage: Q = J([x^2+10*x+2,11])
        sage: K(Q)
        (9 : 1 : 5 : 1)
    """

    def __init__(self, J):
        r"""
        Constructor for a Kummer surface of the Jacobian of a genus-2 curve.

        The equation for the Kummer surface is based on the code provided in
        https://people.maths.ox.ac.uk/flynn/genus2/kummer/ with the modifications
        outlined in [Mue2010]_.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: f = x^5 + x + 1
            sage: C = HyperellipticCurve(f)
            sage: J = Jacobian(C)
            sage: K = KummerSurface(J); K
            Kummer Surface of Jacobian of Hyperelliptic Curve over Rational Field defined by y^2 = x^5 + x + 1.
            The defining equation is X0^4 - 4*X0*X1^3 + 4*X0^2*X1*X2 - 4*X0*X1^2*X2 + 2*X0^2*X2^2 + X2^4 - 4*X0^3*X3 - 2*X0^2*X1*X3 - 2*X1*X2^2*X3 + X1^2*X3^2 - 4*X0*X2*X3^2
        """
        R = J.base_ring()
        self._base_ring = R
        self._curve = J.curve()
        self._jacobian = J

        PP = ProjectiveSpace(3, R, ["X0", "X1", "X2", "X3"])
        X0, X1, X2, X3 = PP.gens()
        f, h = self._curve.hyperelliptic_polynomials()
        [h0, h1, h2, h3] = [h[i] for i in range(4)]
        [f0, f1, f2, f3, f4, f5, f6] = [f[i] for i in range(7)]

        K2 = X1**2 - 4 * X0 * X2

        K1 = (
            (4 * f0 + h0**2) * X0**3
            + (2 * f1 + h0 * h1) * X0**2 * X1
            + (h0 * h2) * X0 * X1**2
            + (h0 * h3) * X1**3
            + (4 * f2 - 2 * h0 * h2 + h1**2) * X0**2 * X2
            + (2 * f3 - 3 * h0 * h3 + h1 * h2) * X0 * X1 * X2
            + (h1 * h3) * X1**2 * X2
            + (4 * f4 - 2 * h1 * h3 + h2**2) * X0 * X2**2
            + (2 * f5 + h2 * h3) * X1 * X2**2
            + (4 * f6 + h3**2) * X2**3
        )

        K0 = (
            (-4 * f0 * f2 - f0 * h1**2 + f1**2 + f1 * h0 * h1 - f2 * h0**2) * X0**4
            + (-4 * f0 * f3 - 2 * f0 * h1 * h2 + f1 * h0 * h2 - f3 * h0**2) * X0**3 * X1
            + (-4 * f0 * f4 - 2 * f0 * h1 * h3 - f0 * h2**2 + f1 * h0 * h3 - f4 * h0**2)
            * X0**2
            * X1**2
            + (-4 * f0 * f5 - 2 * f0 * h2 * h3 - f5 * h0**2) * X0 * X1**3
            + (-4 * f0 * f6 - f0 * h3**2 - f6 * h0**2) * X1**4
            + (
                2 * f0 * h1 * h3
                - 2 * f1 * f3
                - f1 * h0 * h3
                - f1 * h1 * h2
                + 2 * f2 * h0 * h2
                - f3 * h0 * h1
            )
            * X0**3
            * X2
            + (
                4 * f0 * f5
                + 2 * f0 * h2 * h3
                - 4 * f1 * f4
                - f1 * h1 * h3
                - f1 * h2**2
                + 2 * f2 * h0 * h3
                + f3 * h0 * h2
                - 2 * f4 * h0 * h1
                + f5 * h0**2
            )
            * X0**2
            * X1
            * X2
            + (
                8 * f0 * f6
                + 2 * f0 * h3**2
                - 4 * f1 * f5
                - 2 * f1 * h2 * h3
                + f3 * h0 * h3
                - 2 * f5 * h0 * h1
                + 2 * f6 * h0**2
            )
            * X0
            * X1**2
            * X2
            + (-4 * f1 * f6 - f1 * h3**2 - 2 * f6 * h0 * h1) * X1**3 * X2
            + (
                -4 * f0 * f6
                - f0 * h3**2
                + 2 * f1 * f5
                + f1 * h2 * h3
                - 4 * f2 * f4
                - f2 * h2**2
                + f3**2
                + f3 * h0 * h3
                + f3 * h1 * h2
                - f4 * h1**2
                + f5 * h0 * h1
                - f6 * h0**2
            )
            * X0**2
            * X2**2
            + (
                4 * f1 * f6
                + f1 * h3**2
                - 4 * f2 * f5
                - 2 * f2 * h2 * h3
                + f3 * h1 * h3
                + 2 * f4 * h0 * h3
                - f5 * h0 * h2
                - f5 * h1**2
                + 2 * f6 * h0 * h1
            )
            * X0
            * X1
            * X2**2
            + (-4 * f2 * f6 - f2 * h3**2 + f5 * h0 * h3 - 2 * f6 * h0 * h2 - f6 * h1**2)
            * X1**2
            * X2**2
            + (
                -2 * f3 * f5
                - f3 * h2 * h3
                + 2 * f4 * h1 * h3
                - f5 * h0 * h3
                - f5 * h1 * h2
                + 2 * f6 * h0 * h2
            )
            * X0
            * X2**3
            + (-4 * f3 * f6 - f3 * h3**2 + f5 * h1 * h3 - 2 * f6 * h1 * h2) * X1 * X2**3
            + (-4 * f4 * f6 - f4 * h3**2 + f5**2 + f5 * h2 * h3 - f6 * h2**2) * X2**4
        )

        F = K2 * X3**2 - K1 * X3 + K0
        self._defining_equation = F
        AlgebraicScheme_subscheme_projective.__init__(self, PP, F)

        J._kummer_surface = self

    def __repr__(self):
        r"""
        String representation of the Kummer surface.

        EXAMPLES::

            sage: R.<x> = GF(19)[]
            sage: C = HyperellipticCurve(x^5-1)
            sage: J = Jacobian(C)
            sage: K = KummerSurface(J); K
            Kummer Surface of Jacobian of Hyperelliptic Curve over Finite Field of size 19 defined by y^2 = x^5 + 18.
            The defining equation is 4*X0*X1^3 - 4*X0^2*X1*X2 + X2^4 + 4*X0^3*X3 - 2*X1*X2^2*X3 + X1^2*X3^2 - 4*X0*X2*X3^2
        """
        return f"Kummer Surface of {self._jacobian}.\nThe defining equation is {self._defining_equation}"

    def _mumford_to_kummer(self, P):
        r"""
        Given a point P on the Jacobian J of a genus-2 curve,
        return the Kummer coordinates of P.

        TESTS::

            sage: R.<x> = GF(19)[]
            sage: H = HyperellipticCurve(x^6 + 3*x + 5, x^2 + x)
            sage: J = Jacobian(H)
            sage: K = J.kummer_surface()
            sage: P = J([x^2 + 11*x + 13, 18*x + 10])
            sage: K._mumford_to_kummer(P)
            (7, 18, 15, 15)
            sage: Q = J([x, 9])
            sage: K._mumford_to_kummer(Q)
            (0, 1, 0, 18)
            sage: K._mumford_to_kummer(J.zero())
            (0, 0, 0, 1)
            sage: S = J.random_element()
            sage: K(S) == K(-S)
            True
        """

        J = self.jacobian()
        C = J.curve()
        R = J.base_ring()
        f, h = C.hyperelliptic_polynomials()
        [h0, h1, h2, h3] = [h[i] for i in range(4)]
        [f0, f1, f2, f3, f4, f5, f6] = [f[i] for i in range(7)]
        u, v = P.uv()
        [u0, u1, u2] = [u[i] for i in range(3)]
        [v0, v1] = [v[i] for i in range(2)]

        # divisors of the form [P + Q - D_infty]
        if u2 == R.one():
            # auxiliary values
            F0 = (
                2 * u0**3 * f6
                - u0**2 * u1 * f5
                + 2 * u0**2 * f4
                - u0 * u1 * f3
                + 2 * u0 * f2
                - u1 * f1
                + 2 * f0
            )
            y1y2 = u0 * v1**2 - u1 * v0 * v1 + v0**2
            h1h2 = (
                h3**2 * u0**3
                + (-h2 * h3) * u0**2 * u1
                + (h2**2 - 2 * h1 * h3) * u0**2
                + h1 * h3 * u0 * u1**2
                + (-h1 * h2 + 3 * h0 * h3) * u0 * u1
                + (h1**2 - 2 * h0 * h2) * u0
                + (-h0 * h3) * u1**3
                + h0 * h2 * u1**2
                + (-h0 * h1) * u1
                + h0**2
            )

            if u1**2 - 4 * u0 == 0:
                # divisor of the form [2*P - D_infty] with P = [-u1/2,v0]
                # TODO: Is this formula correct?
                denom = 4 * (
                    (4 * f6 + h3**2) * u1**6
                    + (-8 * f5 - 4 * h2 * h3) * u1**5
                    + (16 * f4 + 8 * h1 * h3 + 4 * h2**2) * u1**4
                    + (-32 * f3 - 16 * h0 * h3 - 16 * h1 * h2) * u1**3
                    + (64 * f2 + 32 * h0 * h2 + 16 * h1**2) * u1**2
                    + (-128 * f1 - 64 * h0 * h1) * u1
                    + (256 * f0 + 64 * h0**2)
                )
                x0 = denom
                x1 = -u1 * denom
                x2 = u0 * denom
                x3 = (-1) * (
                    (4 * f4 * f6 + f4 * h3**2 - f5**2 - f5 * h2 * h3 + f6 * h2**2)
                    * u1**8
                    + (
                        -16 * f3 * f6
                        - 4 * f3 * h3**2
                        + 4 * f5 * h1 * h3
                        - 8 * f6 * h1 * h2
                    )
                    * u1**7
                    + (
                        64 * f2 * f6
                        + 16 * f2 * h3**2
                        + 8 * f3 * f5
                        + 4 * f3 * h2 * h3
                        - 8 * f4 * h1 * h3
                        - 12 * f5 * h0 * h3
                        + 4 * f5 * h1 * h2
                        + 24 * f6 * h0 * h2
                        + 16 * f6 * h1**2
                    )
                    * u1**6
                    + (
                        -192 * f1 * f6
                        - 48 * f1 * h3**2
                        - 64 * f2 * f5
                        - 32 * f2 * h2 * h3
                        + 16 * f3 * h1 * h3
                        + 32 * f4 * h0 * h3
                        - 16 * f5 * h0 * h2
                        - 16 * f5 * h1**2
                        - 96 * f6 * h0 * h1
                    )
                    * u1**5
                    + (
                        576 * f0 * f6
                        + 144 * f0 * h3**2
                        + 224 * f1 * f5
                        + 112 * f1 * h2 * h3
                        + 64 * f2 * f4
                        + 16 * f2 * h2**2
                        - 16 * f3**2
                        - 80 * f3 * h0 * h3
                        - 16 * f3 * h1 * h2
                        + 16 * f4 * h1**2
                        + 112 * f5 * h0 * h1
                        + 144 * f6 * h0**2
                    )
                    * u1**4
                    + (
                        -768 * f0 * f5
                        - 384 * f0 * h2 * h3
                        - 256 * f1 * f4
                        - 64 * f1 * h1 * h3
                        - 64 * f1 * h2**2
                        + 128 * f2 * h0 * h3
                        + 64 * f3 * h0 * h2
                        - 128 * f4 * h0 * h1
                        - 192 * f5 * h0**2
                    )
                    * u1**3
                    + (
                        1024 * f0 * f4
                        + 384 * f0 * h1 * h3
                        + 256 * f0 * h2**2
                        + 128 * f1 * f3
                        - 192 * f1 * h0 * h3
                        + 64 * f1 * h1 * h2
                        - 128 * f2 * h0 * h2
                        + 64 * f3 * h0 * h1
                        + 256 * f4 * h0**2
                    )
                    * u1**2
                    + (
                        -1024 * f0 * f3
                        - 512 * f0 * h1 * h2
                        + 256 * f1 * h0 * h2
                        - 256 * f3 * h0**2
                    )
                    * u1
                    + (
                        1024 * f0 * f2
                        + 256 * f0 * h1**2
                        - 256 * f1**2
                        - 256 * f1 * h0 * h1
                        + 256 * f2 * h0**2
                    )
                )
            elif y1y2 == R.zero():
                # divisor of the form [(a,0) + (b,0) - D_infty]
                denom = u1**2 - 4 * u0
                x0 = denom
                x1 = -u1 * denom
                x2 = u0 * denom
                x3 = F0 + h1h2 - y1y2
            else:
                # the generic case
                denom = 4 * y1y2 * (u1**2 - 4 * u0)
                x0 = denom
                x1 = -u1 * denom
                x2 = u0 * denom
                term = (
                    (-(f5**2) + 4 * f4 * f6) * u0**4
                    + (-4 * f3 * f6) * u0**3 * u1
                    + 2 * f3 * f5 * u0**3
                    + 4 * f2 * f6 * u0**2 * u1**2
                    + (-4 * f2 * f5 + 4 * f1 * f6) * u0**2 * u1
                    + (-(f3**2) + 4 * f2 * f4 - 2 * f1 * f5 + 4 * f0 * f6) * u0**2
                    + (-4 * f1 * f6) * u0 * u1**3
                    + (4 * f1 * f5 - 8 * f0 * f6) * u0 * u1**2
                    + (-4 * f1 * f4 + 4 * f0 * f5) * u0 * u1
                    + 2 * f1 * f3 * u0
                    + 4 * f0 * f6 * u1**4
                    + (-4 * f0 * f5) * u1**3
                    + 4 * f0 * f4 * u1**2
                    + (-4 * f0 * f3) * u1
                    - f1**2
                    + 4 * f0 * f2
                )
                x3 = 4 * (F0 + h1h2 - y1y2) * y1y2 - (term * (u1**2 - 4 * u0) + F0**2)
        elif u1 == R.one():
            # In this case, the divisor is of the form [P + O0 - D_{\infty}],
            # where P = (-u0,v0) and O0 is a point at infinity.
            if C.is_ramified():
                O0 = C.points_at_infinity()[0]
            if C.is_split():
                [P_plus, P_minus] = C.points_at_infinity()
                if P._n == 1:
                    O0 = P_plus
                else:
                    O0 = P_minus
            # coordinates of the points in the support of the divisor
            X1, Y1, Z1 = -u0, v0, R.one()
            X2, Y2, Z2 = O0

            x0 = R.zero()
            x1 = X2**3 * Z1**3
            x2 = X1 * X2**3 * Z1**2
            x3 = -(
                2 * Y1 * Y2
                + h3 * Y1 * X2**3
                + Y2 * (h3 * X1**3 + h2 * X1**2 * Z1 + h1 * X1 * Z1**2 + h0 * Z1**3)
                - 2 * f6 * X1**3 * X2**3
                - f5 * X1**2 * X2**3 * Z1
            )
        elif C.is_split() and P._n != 1:
            # divisor [P_infty+ - P_infty-]
            # TODO: Is this the correct image?
            x0 = R.zero()
            x1 = R.zero()
            x2 = -(h3**2) - 4 * f6
            x3 = f6 * h2**2 - f5 * h2 * h3 + f4 * h3**2 - f5**2 + 4 * f4 * f6
        else:  # P = 0
            x0, x1, x2, x3 = R.zero(), R.zero(), R.zero(), R.one()
        return (x0, x1, x2, x3)

    def __call__(self, P):
        r"""
        Create a point on the Kummer surface.

        INPUT: ``P`` -- either a point `P` on the Jacobian or
        coordinates `(X0, X1, X2, X3)` defining a point on the Kummer surface

        OUTPUT: A point on the Kummer surface.

        EXAMPLES::

            sage: R.<x> = FiniteField(13)[]
            sage: H = HyperellipticCurve(x^6 + 11*x^4 + 6*x^3 + 10*x^2 + 11*x + 1)
            sage: J = Jacobian(H)
            sage: K = KummerSurface(J)
            sage: D1 = K([1,11,9,1]); D1
            (1 : 11 : 9 : 1)
            sage: P1 = J(x^2 + 2*x + 9, 4*x + 4)
            sage: K(P1)
            (1 : 11 : 9 : 1)
        """
        if P == 0:
            coords = (0, 0, 0, 1)
        # construct point from Mumford coordinates
        elif isinstance(P, MumfordDivisorClassField):
            coords = self._mumford_to_kummer(P)
        else:
            if not len(P) == 4:
                raise ValueError("The input must consist of 4 coordinates")
            coords = P

        return self.point(coords)

    def defining_equation(self):
        r"""
        Return the defining equation of the Kummer surface in P^3.

        EXAMPLES::

            sage: R.<x> = GF(19)[]
            sage: C = HyperellipticCurve(x^5-1)
            sage: J = Jacobian(C)
            sage: K = KummerSurface(J); K.defining_equation()
            4*X0*X1^3 - 4*X0^2*X1*X2 + X2^4 + 4*X0^3*X3 - 2*X1*X2^2*X3 + X1^2*X3^2 - 4*X0*X2*X3^2
            sage: K.ambient_space()
            Projective Space of dimension 3 over Finite Field of size 19
        """
        return self._defining_equation

    def jacobian(self):
        r"""
        Return the Jacobain of the genus-2 curve lying above the Kummer surface.

        EXAMPLES::

            sage: R.<x> = GF(19)[]
            sage: C = HyperellipticCurve(x^5-1)
            sage: J = Jacobian(C)
            sage: K = KummerSurface(J); K.jacobian()
            Jacobian of Hyperelliptic Curve over Finite Field of size 19 defined by y^2 = x^5 + 18
        """
        return self._jacobian
