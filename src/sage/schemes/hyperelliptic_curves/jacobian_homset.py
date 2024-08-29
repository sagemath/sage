"""
Rational point sets on a Jacobian

EXAMPLES::

    sage: x = QQ['x'].0
    sage: f = x^5 + x + 1
    sage: C = HyperellipticCurve(f); C
    Hyperelliptic Curve over Rational Field defined by y^2 = x^5 + x + 1
    sage: C(QQ)
    Set of rational points of Hyperelliptic Curve over Rational Field
     defined by y^2 = x^5 + x + 1
    sage: P = C([0,1,1])
    sage: J = C.jacobian(); J
    Jacobian of Hyperelliptic Curve over Rational Field defined by y^2 = x^5 + x + 1
    sage: Q = J(QQ)(P); Q
    (x, y - 1)
    sage: Q + Q
    (x^2, y - 1/2*x - 1)
    sage: Q*3
    (x^2 - 1/64*x + 1/8, y + 255/512*x + 65/64)

::

    sage: F.<a> = GF(3)
    sage: R.<x> = F[]
    sage: f = x^5 - 1
    sage: C = HyperellipticCurve(f)
    sage: J = C.jacobian()
    sage: X = J(F)
    sage: a = x^2 - x + 1; b = -x + 1; c = x - 1; d = 0
    sage: D1 = X([a,b]); D1
    (x^2 + 2*x + 1, y + x + 2)
    sage: D2 = X([c,d]); D2
    (x + 2, y)
    sage: D1 + D2
    (x^2 + 2*x + 2, y + 2*x + 1)
"""
# ****************************************************************************
#  Copyright (C) 2006 David Kohel <kohel@maths.usyd.edu>
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.misc.lazy_import import lazy_import
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.polynomial_element import Polynomial
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.schemes.generic.homset import SchemeHomset_points
from sage.schemes.hyperelliptic_curves.jacobian_morphism import JacobianMorphism_divisor_class_field

lazy_import('sage.schemes.generic.morphism', 'SchemeMorphism')


class JacobianHomset_divisor_classes(SchemeHomset_points):
    def __init__(self, Y, X, **kwds):
        R = X.base_ring()
        S = Y.coordinate_ring()
        SchemeHomset_points.__init__(self, Y, X, **kwds)
        P2 = X.curve()._printing_ring
        if S != R:
            y = str(P2.gen())
            x = str(P2.base_ring().gen())
            P1 = PolynomialRing(S, name=x)
            P2 = PolynomialRing(P1, name=y)
        self._printing_ring = P2

    def __call__(self, P):
        r"""
        Return a rational point P in the abstract Homset J(K), given:

        0. A point P in J = Jac(C), returning P;

        1. A point P on the curve C such that J = Jac(C), where C is
           an odd degree model, returning [P - oo];

        2. A pair of points (P, Q) on the curve C such that J = Jac(C),
           returning [P-Q];

        3. A list of polynomials (a,b) such that `b^2 + h*b - f = 0 mod a`,
           returning [(a(x),y-b(x))].

        EXAMPLES::

            sage: P.<x> = PolynomialRing(QQ)
            sage: f = x^5 - x + 1; h = x
            sage: C = HyperellipticCurve(f,h,'u,v')
            sage: P = C(0,1,1)
            sage: J = C.jacobian()
            sage: Q = J(QQ)(P)
            sage: for i in range(6): i*Q
            (1)
            (u, v - 1)
            (u^2, v + u - 1)
            (u^2, v + 1)
            (u, v + 1)
            (1)

        ::

            sage: F.<a> = GF(3)
            sage: R.<x> = F[]
            sage: f = x^5 - 1
            sage: C = HyperellipticCurve(f)
            sage: J = C.jacobian()
            sage: X = J(F)
            sage: a = x^2 - x + 1; b = -x + 1; c = x - 1; d = 0
            sage: D1 = X([a,b]); D1
            (x^2 + 2*x + 1, y + x + 2)
            sage: D2 = X([c,d]); D2
            (x + 2, y)
            sage: D1 + D2
            (x^2 + 2*x + 2, y + 2*x + 1)
        """
        if isinstance(P, (Integer, int)) and P == 0:
            R = PolynomialRing(self.value_ring(), 'x')
            return JacobianMorphism_divisor_class_field(self,
                                                        (R.one(), R.zero()))
        elif isinstance(P, (list, tuple)):
            if len(P) == 1 and P[0] == 0:
                R = PolynomialRing(self.value_ring(), 'x')
                return JacobianMorphism_divisor_class_field(self,
                                                            (R.one(), R.zero()))
            elif len(P) == 2:
                P1 = P[0]
                P2 = P[1]
                if isinstance(P1, Integer) and isinstance(P2, Integer):
                    R = PolynomialRing(self.value_ring(), 'x')
                    P1 = R(P1)
                    P2 = R(P2)
                    return JacobianMorphism_divisor_class_field(self, (P1, P2))
                if isinstance(P1, Integer) and isinstance(P2, Polynomial):
                    R = PolynomialRing(self.value_ring(), 'x')
                    P1 = R(P1)
                    return JacobianMorphism_divisor_class_field(self, (P1, P2))
                if isinstance(P2, Integer) and isinstance(P1, Polynomial):
                    R = PolynomialRing(self.value_ring(), 'x')
                    P2 = R(P2)
                    return JacobianMorphism_divisor_class_field(self, (P1, P2))
                if isinstance(P1, Polynomial) and isinstance(P2, Polynomial):
                    return JacobianMorphism_divisor_class_field(self, tuple(P))
                if isinstance(P1, SchemeMorphism) and isinstance(P2, SchemeMorphism):
                    return self(P1) - self(P2)
            raise TypeError("argument P (= %s) must have length 2" % P)
        elif isinstance(P, JacobianMorphism_divisor_class_field) and self == P.parent():
            return P
        elif isinstance(P, SchemeMorphism):
            x0 = P[0]
            y0 = P[1]
            R, x = PolynomialRing(self.value_ring(), 'x').objgen()
            return self((x - x0, R(y0)))
        raise TypeError("argument P (= %s) does not determine a divisor class" % P)

    def _morphism(self, *args, **kwds):
        return JacobianMorphism_divisor_class_field(*args, **kwds)

    def curve(self):
        return self.codomain().curve()

    def value_ring(self):
        """
        Return S for a homset X(T) where T = Spec(S).
        """
        from sage.schemes.generic.scheme import AffineScheme
        T = self.domain()
        if isinstance(T, AffineScheme):
            return T.coordinate_ring()
        else:
            raise TypeError("domain of argument must be of the form Spec(S)")

    def base_extend(self, R):
        if R != ZZ:
            raise NotImplementedError("Jacobian point sets viewed as modules over rings other than ZZ not implemented")
        return self
