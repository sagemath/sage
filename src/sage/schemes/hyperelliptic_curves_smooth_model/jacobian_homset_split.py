from sage.rings.integer import Integer
from sage.rings.polynomial.polynomial_element import Polynomial
from sage.schemes.hyperelliptic_curves_smooth_model.jacobian_homset_generic import (
    HyperellipticJacobianHomset,
)
from sage.schemes.hyperelliptic_curves_smooth_model.jacobian_morphism import (
    MumfordDivisorClassFieldSplit,
)

"""
TODO

should we make a hyperelliptic point class?
at the moment, this is the type we get from calling a point from the projective model
"""
from sage.schemes.toric.morphism import SchemeMorphism_point_toric_field


class HyperellipticJacobianHomsetSplit(HyperellipticJacobianHomset):
    def __init__(self, Y, X, **kwds):
        super().__init__(Y, X, **kwds)
        self._morphism_element = MumfordDivisorClassFieldSplit

    def zero(self):
        """
        Return the zero element of the Jacobian
        """
        g = self.curve().genus()
        R = self.curve().polynomial_ring()
        n = (g / 2).ceil()
        return self._morphism_element(self, R.one(), R.zero(), n)

    def point_to_mumford_coordinates(self, P):
        """
        On input a point P, return the Mumford coordinates
        of (the affine part of) the divisor [P] and an integer n,
        where 
        * n = 1 if P is the point oo+ 
        * n = 0 otherwise .
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
        """
        Return a rational point in the abstract Homset J(K), given:

        0. The integer 0, returning 0 in J;

        1. A point P in J = Jac(C), returning P;

        2. A point P on the curve H such that J = Jac(H),
           returning [P - P0], where P0 is the distinguished point of H
           (by default P0 = oo+)

        2. Two points P, Q on the curve H such that J = Jac(H),
           returning [P-Q];

        3. Polynomials u,v such that `v^2 + h*v - f = 0 mod u`,
           returning J(u,v,n) with n = ((g - deg(u))/2).ceil().

        4. Polynomials u,v and an integer n such that `v^2 + h*v - f = 0 mod u`
           and 0 <= n <= g/2,
           returning J(u,v,n).


        EXAMPLES::

            sage: from hyperelliptic_constructor import HyperellipticCurveSmoothModel
            sage: R.<x> = PolynomialRing(GF(13))
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

        The points at infinity may also be embedded into the Jacobian.
            sage: [P0,P1] = H.points_at_infinity()
            sage: JH(P0)
            (1, 0 : 2)
            sage: JH(P1)
            (1, 0 : 1)

        TODO: code reuse
        """

        g = self.curve().genus()

        if len(args) == 1 and isinstance(args[0], (list,tuple)):
            args = args[0] #unpack
        # if len(args) == 1:
        #     P1 = args[0]
        #     if P1 == 0:
        #         # TODO: what is R here?!
        #         u = R.one()
        #         v = R.zero()
        #         n = (g / 2).ceil()
        #     elif isinstance(P1, self._morphism_element):
        #         return P1
        #     elif isinstance(P1, SchemeMorphism_point_toric_field):
        #         args = args + (self.curve().distinguished_point(),) #this case will now be handeled below.
        if len(args) == 2 or len(args) == 3:
            P1 = args[0]
            P2 = args[1]
            if isinstance(P1, SchemeMorphism_point_toric_field) and isinstance(P2, SchemeMorphism_point_toric_field):
                u1, v1, n1 = self.point_to_mumford_coordinates(P1)
                P2_inv = self.curve().hyperelliptic_involution(P2)
                u2, v2, n2 = self.point_to_mumford_coordinates(P2_inv)
                u, v, _ = self._cantor_composition_generic(u1, v1, u2, v2)
                n = (g / 2).ceil() - 1 + n1 + n2 #this solution is a bit hacky
                if len(args) == 3:
                    raise ValueError("Only one or two points are allowed as input.")
            elif isinstance(P1, Polynomial) and isinstance(P2, Polynomial):
                u = P1
                v = P2
                if len(args) == 3 and isinstance(args[2], Integer):
                    n = args[2]
                else:
                    n = ((g - u.degree()) / 2).ceil() #TODO: do we really want to allow this input?
            else:
                raise ValueError("The input must consist of points or polynomials.")
        if len(args) > 3:
            raise ValueError("At most three arguments are allowed as input.")

        return self._morphism_element(self, u, v, n=n, check=check)
        

    def cantor_composition(self, u1, v1, n1, u2, v2, n2):
        """
        Follows algorithm 3.4 of

        Efficient Arithmetic on Hyperelliptic Curves With Real Representation
        David J. Mireles Morales (2008)
        https://www.math.auckland.ac.nz/~sgal018/Dave-Mireles-Full.pdf

        TODO: when h = 0 we can speed this up.
        """
        # Collect data from HyperellipticCurve
        H = self.curve()
        g = H.genus()

        # Cantor composition
        u3, v3, s_deg = self._cantor_composition_generic(u1, v1, u2, v2)

        # Compute new weight
        n3 = n1 + n2 + s_deg - (g / 2).ceil()

        return u3, v3, n3

    def cantor_reduction(self, u0, v0, n0):
        """
        Follows algorithm 3.5 of

        Efficient Arithmetic on Hyperelliptic Curves With Real Representation
        David J. Mireles Morales (2008)
        https://www.math.auckland.ac.nz/~sgal018/Dave-Mireles-Full.pdf
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
        """
        Follows algorithm 3.6 of

        Efficient Arithmetic on Hyperelliptic Curves With Real Representation
        David J. Mireles Morales (2008)
        https://www.math.auckland.ac.nz/~sgal018/Dave-Mireles-Full.pdf
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
