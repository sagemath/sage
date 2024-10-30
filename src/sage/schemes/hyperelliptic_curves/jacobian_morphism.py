r"""
Jacobian 'morphism' as a class in the Picard group

This module implements the group operation in the Picard group of a
hyperelliptic curve, represented as divisors in Mumford
representation, using Cantor's algorithm.

A divisor on the hyperelliptic curve `y^2 + y h(x) = f(x)`
is stored in Mumford representation, that is, as two polynomials
`u(x)` and `v(x)` such that:

- `u(x)` is monic,

- `u(x)` divides `f(x) - h(x) v(x) - v(x)^2`,

- `deg(v(x)) < deg(u(x)) \le g`.

REFERENCES:

A readable introduction to divisors, the Picard group, Mumford
representation, and Cantor's algorithm:

- J. Scholten, F. Vercauteren. An Introduction to Elliptic and
  Hyperelliptic Curve Cryptography and the NTRU Cryptosystem. To
  appear in B. Preneel (Ed.) State of the Art in Applied Cryptography
  - COSIC '03, Lecture Notes in Computer Science, Springer 2004.

A standard reference in the field of cryptography:

- R. Avanzi, H. Cohen, C. Doche, G. Frey, T. Lange, K. Nguyen, and F.
  Vercauteren, Handbook of Elliptic and Hyperelliptic Curve
  Cryptography. CRC Press, 2005.

EXAMPLES: The following curve is the reduction of a curve whose
Jacobian has complex multiplication.

::

    sage: x = GF(37)['x'].gen()
    sage: H = HyperellipticCurve(x^5 + 12*x^4 + 13*x^3 + 15*x^2 + 33*x); H
    Hyperelliptic Curve over Finite Field of size 37 defined
     by y^2 = x^5 + 12*x^4 + 13*x^3 + 15*x^2 + 33*x

At this time, Jacobians of hyperelliptic curves are handled
differently than elliptic curves::

    sage: J = H.jacobian(); J
    Jacobian of Hyperelliptic Curve over Finite Field of size 37 defined
     by y^2 = x^5 + 12*x^4 + 13*x^3 + 15*x^2 + 33*x
    sage: J = J(J.base_ring()); J
    Set of rational points of Jacobian of Hyperelliptic Curve over Finite Field
     of size 37 defined by y^2 = x^5 + 12*x^4 + 13*x^3 + 15*x^2 + 33*x

Points on the Jacobian are represented by Mumford's polynomials.
First we find a couple of points on the curve::

    sage: P1 = H.lift_x(2); P1
    (2 : 11 : 1)
    sage: Q1 = H.lift_x(10); Q1
    (10 : 18 : 1)

Observe that 2 and 10 are the roots of the polynomials in x,
respectively::

    sage: P = J(P1); P
    (x + 35, y + 26)
    sage: Q = J(Q1); Q
    (x + 27, y + 19)

::

    sage: P + Q
    (x^2 + 25*x + 20, y + 13*x)
    sage: (x^2 + 25*x + 20).roots(multiplicities=False)
    [10, 2]

Frobenius satisfies

.. MATH::

    x^4 + 12*x^3 + 78*x^2 + 444*x + 1369

on the Jacobian of this reduction and the order of the Jacobian is
`N = 1904`.

::

    sage: 1904*P
    (1)
    sage: 34*P == 0
    True
    sage: 35*P == P
    True
    sage: 33*P == -P
    True

::

    sage: Q*1904
    (1)
    sage: Q*238 == 0
    True
    sage: Q*239 == Q
    True
    sage: Q*237 == -Q
    True
"""

#*****************************************************************************
#  Copyright (C) 2005 David Kohel <kohel@maths.usyd.edu.au>
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.latex import latex

from sage.structure.element import AdditiveGroupElement
from sage.structure.richcmp import richcmp, op_NE
from sage.schemes.generic.morphism import SchemeMorphism


def cantor_reduction_simple(a, b, f, genus):
    r"""
    Return the unique reduced divisor linearly equivalent to
    `(a, b)` on the curve `y^2 = f(x).`

    See the docstring of
    :mod:`sage.schemes.hyperelliptic_curves.jacobian_morphism` for
    information about divisors, linear equivalence, and reduction.

    EXAMPLES::

        sage: x = QQ['x'].gen()
        sage: f = x^5 - x
        sage: H = HyperellipticCurve(f); H
        Hyperelliptic Curve over Rational Field defined by y^2 = x^5 - x
        sage: J = H.jacobian()(QQ); J
        Set of rational points of Jacobian of Hyperelliptic Curve over Rational Field
         defined by y^2 = x^5 - x

    The following point is 2-torsion::

        sage: P = J(H.lift_x(-1)); P
        (x + 1, y)
        sage: 2 * P # indirect doctest
        (1)
    """
    a2 = (f - b**2) // a
    a2 = a2.monic()
    b2 = -b % (a2)
    if a2.degree() == a.degree():
        # XXX
        assert a2.degree() == genus + 1
        print("Returning ambiguous form of degree genus+1.")
        return (a2, b2)
    elif a2.degree() > genus:
        return cantor_reduction_simple(a2, b2, f, genus)
    return (a2, b2)


def cantor_reduction(a, b, f, h, genus):
    r"""
    Return the unique reduced divisor linearly equivalent to
    `(a, b)` on the curve `y^2 + y h(x) = f(x)`.

    See the docstring of
    :mod:`sage.schemes.hyperelliptic_curves.jacobian_morphism` for
    information about divisors, linear equivalence, and reduction.

    EXAMPLES::

        sage: x = QQ['x'].gen()
        sage: f = x^5 - x
        sage: H = HyperellipticCurve(f, x); H
        Hyperelliptic Curve over Rational Field defined by y^2 + x*y = x^5 - x
        sage: J = H.jacobian()(QQ); J
        Set of rational points of Jacobian of Hyperelliptic Curve over
         Rational Field defined by y^2 + x*y = x^5 - x

    The following point is 2-torsion::

        sage: Q = J(H.lift_x(0)); Q
        (x, y)
        sage: 2*Q # indirect doctest
        (1)

    The next point is not 2-torsion::

        sage: P = J(H.lift_x(-1)); P
        (x + 1, y)
        sage: 2 * J(H.lift_x(-1)) # indirect doctest
        (x^2 + 2*x + 1, y + 4*x + 4)
        sage: 3 * J(H.lift_x(-1)) # indirect doctest
        (x^2 - 487*x - 324, y + 10755*x + 7146)
    """
    assert a.degree() < 2*genus+1
    assert b.degree() < a.degree()
    k = f - h*b - b**2
    if 2*a.degree() == k.degree():
        # must adjust b to include the point at infinity
        g1 = a.degree()
        x = a.parent().gen()
        r = (x**2 + h[g1]*x - f[2*g1]).roots()[0][0]
        b = b + r*(x**g1 - (x**g1) % (a))
        k = f - h*b - b**2
    assert k % (a) == 0
    a = (k // a).monic()
    b = -(b+h) % (a)
    if a.degree() > genus:
        return cantor_reduction(a, b, f, h, genus)
    return (a, b)


def cantor_composition_simple(D1,D2,f,genus):
    r"""
    Given `D_1` and `D_2` two reduced Mumford
    divisors on the Jacobian of the curve `y^2 = f(x)`,
    computes a representative `D_1 + D_2`.

    .. warning::

       The representative computed is NOT reduced! Use
       :func:`cantor_reduction_simple` to reduce it.

    EXAMPLES::

        sage: x = QQ['x'].gen()
        sage: f = x^5 + x
        sage: H = HyperellipticCurve(f); H
        Hyperelliptic Curve over Rational Field defined by y^2 = x^5 + x

    ::

        sage: F.<a> = NumberField(x^2 - 2, 'a')                                         # needs sage.rings.number_field
        sage: J = H.jacobian()(F); J                                                    # needs sage.rings.number_field
        Set of rational points of Jacobian of Hyperelliptic Curve over
         Number Field in a with defining polynomial x^2 - 2 defined by y^2 = x^5 + x

    ::

        sage: # needs sage.rings.number_field
        sage: P = J(H.lift_x(F(1))); P
        (x - 1, y + a)
        sage: Q = J(H.lift_x(F(0))); Q
        (x, y)
        sage: 2*P + 2*Q # indirect doctest
        (x^2 - 2*x + 1, y + 3/2*a*x - 1/2*a)
        sage: 2*(P + Q) # indirect doctest
        (x^2 - 2*x + 1, y + 3/2*a*x - 1/2*a)
        sage: 3*P # indirect doctest
        (x^2 - 25/32*x + 49/32, y + 45/256*a*x + 315/256*a)
    """
    a1, b1 = D1
    a2, b2 = D2
    if a1 == a2 and b1 == b2:
        # Duplication law:
        d, h1, h3 = a1.xgcd(2*b1)
        a = (a1 // d)**2
        b = (b1 + h3*((f - b1**2) // d)) % (a)
    else:
        d0, _, h2 = a1.xgcd(a2)
        if d0 == 1:
            a = a1*a2
            b = (b2 + h2*a2*(b1-b2)) % (a)
        else:
            d, l, h3 = d0.xgcd(b1 + b2)
            a = (a1*a2) // (d**2)
            b = ((b2 + l*h2*(b1-b2)*(a2 // d)) + h3*((f - b2**2) // d)) % (a)
    a = a.monic()
    return (a, b)


def cantor_composition(D1,D2,f,h,genus):
    r"""
    EXAMPLES::

        sage: # needs sage.rings.finite_rings
        sage: F.<a> = GF(7^2, 'a')
        sage: x = F['x'].gen()
        sage: f = x^7 + x^2 + a
        sage: H = HyperellipticCurve(f, 2*x); H
        Hyperelliptic Curve over Finite Field in a of size 7^2
         defined by y^2 + 2*x*y = x^7 + x^2 + a
        sage: J = H.jacobian()(F); J
        Set of rational points of Jacobian of Hyperelliptic Curve over
         Finite Field in a of size 7^2 defined by y^2 + 2*x*y = x^7 + x^2 + a

    ::

        sage: Q = J(H.lift_x(F(1))); Q                                                  # needs sage.rings.finite_rings
        (x + 6, y + 5*a)
        sage: 10*Q  # indirect doctest                                                  # needs sage.rings.finite_rings
        (x^3 + (3*a + 1)*x^2 + (2*a + 5)*x + a + 5, y + (3*a + 2)*x^2 + (6*a + 1)*x + a + 4)
        sage: 7*8297*Q                                                                  # needs sage.rings.finite_rings
        (1)

    ::

        sage: Q = J(H.lift_x(F(a+1))); Q                                                # needs sage.rings.finite_rings
        (x + 6*a + 6, y + 2)
        sage: 7*8297*Q  # indirect doctest                                              # needs sage.rings.finite_rings
        (1)

        A test over a prime field:

        sage: # needs sage.rings.finite_rings
        sage: F = GF(next_prime(10^30))
        sage: x = F['x'].gen()
        sage: f = x^7 + x^2 + 1
        sage: H = HyperellipticCurve(f, 2*x); H
        Hyperelliptic Curve over Finite Field of size 1000000000000000000000000000057
         defined by y^2 + 2*x*y = x^7 + x^2 + 1
        sage: J = H.jacobian()(F); J
        Set of rational points of Jacobian of Hyperelliptic Curve
         over Finite Field of size 1000000000000000000000000000057
         defined by y^2 + 2*x*y = x^7 + x^2 + 1
        sage: Q = J(H.lift_x(F(1))); Q
        (x + 1000000000000000000000000000056, y + 1000000000000000000000000000056)
        sage: 10*Q  # indirect doctest
        (x^3 + 150296037169838934997145567227*x^2
             + 377701248971234560956743242408*x + 509456150352486043408603286615,
         y + 514451014495791237681619598519*x^2
           + 875375621665039398768235387900*x + 861429240012590886251910326876)
        sage: 7*8297*Q
        (x^3 + 35410976139548567549919839063*x^2
             + 26230404235226464545886889960*x + 681571430588959705539385624700,
         y + 999722365017286747841221441793*x^2
           + 262703715994522725686603955650*x + 626219823403254233972118260890)
    """
    a1, b1 = D1
    a2, b2 = D2
    if a1 == a2 and b1 == b2:
        # Duplication law:
        d, h1, h3 = a1.xgcd(2*b1 + h)
        a = (a1 // d)**2
        b = (b1 + h3*((f-h*b1-b1**2) // d)) % (a)
    else:
        d0, _, h2 = a1.xgcd(a2)
        if d0 == 1:
            a = a1 * a2
            b = (b2 + h2*a2*(b1-b2)) % (a)
        else:
            e0 = b1+b2+h
            if e0 == 0:
                a = (a1*a2) // (d0**2)
                b = (b2 + h2*(b1-b2)*(a2 // d0)) % (a)
            else:
                d, l, h3 = d0.xgcd(e0)
                a = (a1*a2) // (d**2)
                b = (b2 + l*h2*(b1-b2)*(a2 // d) + h3*((f-h*b2-b2**2) // d)) % (a)
    a = a.monic()
    return (a, b)


class JacobianMorphism_divisor_class_field(AdditiveGroupElement, SchemeMorphism):
    r"""
    An element of a Jacobian defined over a field, i.e. in
    `J(K) = \mathrm{Pic}^0_K(C)`.
    """
    def __init__(self, parent, polys, check=True):
        r"""
        Create a new Jacobian element in Mumford representation.

        INPUT:

        - ``parent`` -- the parent Homset
        - ``polys`` -- Mumford's `u` and `v` polynomials
        - ``check`` -- boolean (default: ``True``); if ``True``, ensure that
          polynomials define a divisor on the appropriate curve and are reduced

        .. warning::

           Not for external use! Use ``J(K)([u, v])`` instead.

        EXAMPLES::

            sage: x = GF(37)['x'].gen()
            sage: H = HyperellipticCurve(x^5 + 12*x^4 + 13*x^3 + 15*x^2 + 33*x)
            sage: J = H.jacobian()(GF(37));  J
            Set of rational points of Jacobian of Hyperelliptic Curve over
             Finite Field of size 37 defined by
            y^2 = x^5 + 12*x^4 + 13*x^3 + 15*x^2 + 33*x

        ::

            sage: P1 = J(H.lift_x(2)); P1  # indirect doctest
            (x + 35, y + 26)
            sage: P1.parent()
            Set of rational points of Jacobian of Hyperelliptic Curve over
             Finite Field of size 37 defined by
            y^2 = x^5 + 12*x^4 + 13*x^3 + 15*x^2 + 33*x
            sage: type(P1)
            <class 'sage.schemes.hyperelliptic_curves.jacobian_morphism.JacobianMorphism_divisor_class_field'>
        """
        SchemeMorphism.__init__(self, parent)
        if check:
            C = parent.curve()
            f, h = C.hyperelliptic_polynomials()
            a, b = polys
            if not (b**2 + h*b - f) % a == 0:
                raise ValueError("Argument polys (= %s) must be divisor on curve %s." % (
                    polys, C))
            genus = C.genus()
            if a.degree() > genus:
                polys = cantor_reduction(a, b, f, h, genus)
        self.__polys = polys

    def _printing_polys(self):
        r"""
        Internal function formatting Mumford polynomials for printing.

        TESTS::

            sage: # needs sage.rings.finite_rings
            sage: F.<a> = GF(7^2, 'a')
            sage: x = F['x'].gen()
            sage: f = x^7 + x^2 + a
            sage: H = HyperellipticCurve(f, 2*x)
            sage: J = H.jacobian()(F)

        ::

            sage: Q = J(H.lift_x(F(1))); Q  # indirect doctest                          # needs sage.rings.finite_rings
            (x + 6, y + 5*a)
        """
        a, b = self.__polys
        P = self.parent()._printing_ring
        y = P.gen()
        x = P.base_ring().gen()
        return (a(x), y - b(x))

    def _repr_(self):
        r"""
        Return a string representation of this Mumford divisor.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: F.<a> = GF(7^2, 'a')
            sage: x = F['x'].gen()
            sage: f = x^7 + x^2 + a
            sage: H = HyperellipticCurve(f, 2*x)
            sage: J = H.jacobian()(F)

        ::

            sage: Q = J(0); Q  # indirect doctest                                       # needs sage.rings.finite_rings
            (1)
            sage: Q = J(H.lift_x(F(1))); Q  # indirect doctest                          # needs sage.rings.finite_rings
            (x + 6, y + 5*a)
            sage: Q + Q  # indirect doctest                                             # needs sage.rings.finite_rings
            (x^2 + 5*x + 1, y + (4*a + 2)*x + a + 5)
        """
        if self.is_zero():
            return "(1)"
        a, b = self._printing_polys()
        return "(%s, %s)" % (a, b)

    def _latex_(self):
        r"""
        Return a LaTeX string representing this Mumford divisor.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: F.<alpha> = GF(7^2)
            sage: x = F['x'].gen()
            sage: f = x^7 + x^2 + alpha
            sage: H = HyperellipticCurve(f, 2*x)
            sage: J = H.jacobian()(F)

        ::

            sage: Q = J(0); print(latex(Q))  # indirect doctest                         # needs sage.rings.finite_rings
            \left(1\right)
            sage: Q = J(H.lift_x(F(1))); print(latex(Q))  # indirect doctest            # needs sage.rings.finite_rings
            \left(x + 6, y + 5 \alpha\right)

        ::

            sage: print(latex(Q + Q))                                                   # needs sage.rings.finite_rings
            \left(x^{2} + 5 x + 1, y + \left(4 \alpha + 2\right) x + \alpha + 5\right)
        """
        if self.is_zero():
            return "\\left(1\\right)"
        a, b = self._printing_polys()
        return "\\left(%s, %s\\right)" % (latex(a), latex(b))

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

            sage: # needs sage.rings.number_field
            sage: x = QQ['x'].gen()
            sage: f = x^5 + x
            sage: H = HyperellipticCurve(f)
            sage: F.<a> = NumberField(x^2 - 2, 'a')
            sage: J = H.jacobian()(F); J
            Set of rational points of Jacobian of Hyperelliptic Curve
             over Number Field in a with defining polynomial x^2 - 2
             defined by y^2 = x^5 + x
            sage: P = J(H.lift_x(F(1)))
            sage: P.scheme()
            Jacobian of Hyperelliptic Curve over Rational Field defined by y^2 = x^5 + x
        """
        return self.codomain()

    def point_of_jacobian_of_curve(self):
        r"""
        Return the point in the Jacobian of the curve.

        The Jacobian is the one attached to the projective curve
        corresponding to this hyperelliptic curve.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(GF(11))
            sage: f = x^6 + x + 1
            sage: H = HyperellipticCurve(f)
            sage: J = H.jacobian()
            sage: D = J(H.lift_x(1))
            sage: D  # divisor in Mumford representation
            (x + 10, y + 6)
            sage: jacobian_order = sum(H.frobenius_polynomial())
            sage: jacobian_order
            234
            sage: p = D.point_of_jacobian_of_curve(); p
            [Place (1/x0, 1/x0^3*x1 + 1)
             + Place (x0 + 10, x1 + 6)]
            sage: p  # Jacobian point represented by an effective divisor
            [Place (1/x0, 1/x0^3*x1 + 1)
             + Place (x0 + 10, x1 + 6)]
            sage: p.order()
            39
            sage: 234*p == 0
            True
            sage: G = p.parent()
            sage: G
            Group of rational points of Jacobian over Finite Field of size 11 (Hess model)
            sage: J = G.parent()
            sage: J
            Jacobian of Projective Plane Curve over Finite Field of size 11
             defined by x0^6 + x0^5*x1 + x1^6 - x0^4*x2^2 (Hess model)
            sage: C = J.curve()
            sage: C
            Projective Plane Curve over Finite Field of size 11
             defined by x0^6 + x0^5*x1 + x1^6 - x0^4*x2^2
            sage: C.affine_patch(0) == H.affine_patch(2)
            True
        """
        from sage.schemes.curves.constructor import Curve
        C = self.parent().curve()
        P = C.ambient_space()  # projective plane
        x0, x1, x2 = P.gens()

        # X is the curve positioned in the ambient space
        # such that x1 = x and x2 = y
        X = Curve(C.defining_ideal().gens(), P)
        X = X.affine_patch(2).projective_closure()

        u0, v0 = list(self)
        u1 = u0.subs(x1).homogenize(x0)
        v1 = (x2 - v0.subs(x1)).homogenize(x0)
        u2 = u1/x0**u1.degree()
        v2 = v1/x0**v1.degree()
        u = X.function(u2)
        v = X.function(v2)

        F = X.function_field()
        O = F.maximal_order()
        D = O.ideal([u,v]).divisor()

        Pinf = F.places_infinite()[0]
        assert Pinf.degree() == 1, "no rational point at infinity"

        J = X.jacobian(model='hess', base_div=F.genus()*Pinf)
        G = J.group(self.base_ring())
        return G(D - D.degree()*Pinf)

    def __list__(self):
        r"""
        Return a list `(a(x), b(x))` of the polynomials giving the
        Mumford representation of ``self``.

        TESTS::

            sage: x = QQ['x'].gen()
            sage: f = x^5 + x
            sage: H = HyperellipticCurve(f)
            sage: F.<a> = NumberField(x^2 - 2, 'a')                                     # needs sage.rings.number_field
            sage: J = H.jacobian()(F); J                                                # needs sage.rings.number_field
            Set of rational points of Jacobian of Hyperelliptic Curve
             over Number Field in a with defining polynomial x^2 - 2
             defined by y^2 = x^5 + x

        ::

            sage: P = J(H.lift_x(F(1)))                                                 # needs sage.rings.number_field
            sage: list(P)  # indirect doctest                                           # needs sage.rings.number_field
            [x - 1, -a]
        """
        return list(self.__polys)

    def __tuple__(self):
        r"""
        Return a tuple `(a(x), b(x))` of the polynomials giving the
        Mumford representation of ``self``.

        TESTS::

            sage: x = QQ['x'].gen()
            sage: f = x^5 + x
            sage: H = HyperellipticCurve(f)
            sage: F.<a> = NumberField(x^2 - 2, 'a')                                     # needs sage.rings.number_field
            sage: J = H.jacobian()(F); J                                                # needs sage.rings.number_field
            Set of rational points of Jacobian of Hyperelliptic Curve
             over Number Field in a with defining polynomial x^2 - 2
             defined by y^2 = x^5 + x

        ::

            sage: P = J(H.lift_x(F(1)))                                                 # needs sage.rings.number_field
            sage: tuple(P)  # indirect doctest                                          # needs sage.rings.number_field
            (x - 1, -a)
        """
        return tuple(self.__polys)

    def __getitem__(self, n):
        r"""
        Return the `n`-th item of the pair `(a(x), b(x))`
        of polynomials giving the Mumford representation of ``self``.

        TESTS::

            sage: x = QQ['x'].gen()
            sage: f = x^5 + x
            sage: H = HyperellipticCurve(f)
            sage: F.<a> = NumberField(x^2 - 2, 'a')                                     # needs sage.rings.number_field
            sage: J = H.jacobian()(F); J                                                # needs sage.rings.number_field
            Set of rational points of Jacobian of Hyperelliptic Curve
             over Number Field in a with defining polynomial x^2 - 2
             defined by y^2 = x^5 + x

        ::

            sage: # needs sage.rings.number_field
            sage: P = J(H.lift_x(F(1)))
            sage: P[0] # indirect doctest
            x - 1
            sage: P[1] # indirect doctest
            -a
            sage: P[-1] # indirect doctest
            -a
            sage: P[:1] # indirect doctest
            [x - 1]
        """
        return list(self.__polys)[n]

    def _richcmp_(self, other, op):
        r"""
        Compare ``self`` and ``other``.

        TESTS::

            sage: x = QQ['x'].gen()
            sage: f = x^5 - x
            sage: H = HyperellipticCurve(f); H
            Hyperelliptic Curve over Rational Field defined by y^2 = x^5 - x
            sage: J = H.jacobian()(QQ); J
            Set of rational points of Jacobian of Hyperelliptic Curve over
            Rational Field defined by y^2 = x^5 - x

        The following point is 2-torsion::

            sage: P = J(H.lift_x(-1)); P
            (x + 1, y)
            sage: 0 == 2 * P # indirect doctest
            True
            sage: P == P
            True

        ::

            sage: Q = J(H.lift_x(-1))
            sage: Q == P
            True

        ::

            sage: 2 == Q
            False
            sage: P == False
            False

        Let's verify the same "points" on different schemes are not equal::

            sage: x = QQ['x'].gen()
            sage: f = x^5 + x
            sage: H2 = HyperellipticCurve(f)
            sage: J2 = H2.jacobian()(QQ)

        ::

            sage: P1 = J(H.lift_x(0)); P1
            (x, y)
            sage: P2 = J2(H2.lift_x(0)); P2
            (x, y)
            sage: P1 == P2
            False
        """
        if self.scheme() != other.scheme():
            return op == op_NE
        # since divisors are internally represented as Mumford divisors,
        # comparing polynomials is well-defined
        return richcmp(self.__polys, other.__polys, op)

    def __bool__(self):
        r"""
        Return ``True`` if this divisor is not the additive identity element.

        EXAMPLES::

            sage: x = GF(37)['x'].gen()
            sage: H = HyperellipticCurve(x^5 + 12*x^4 + 13*x^3 + 15*x^2 + 33*x)
            sage: J = H.jacobian()(GF(37))

        ::

            sage: P1 = J(H.lift_x(2)); P1
            (x + 35, y + 26)
            sage: P1 == 0  # indirect doctest
            False
            sage: P1 - P1 == 0  # indirect doctest
            True
        """
        return self.__polys[0] != 1

    def __neg__(self):
        r"""
        Return the additive inverse of this divisor.

        EXAMPLES::

            sage: x = GF(37)['x'].gen()
            sage: H = HyperellipticCurve(x^5 + 12*x^4 + 13*x^3 + 15*x^2 + 33*x)
            sage: J = H.jacobian()(GF(37))
            sage: P1 = J(H.lift_x(2)); P1
            (x + 35, y + 26)
            sage: - P1  # indirect doctest
            (x + 35, y + 11)
            sage: P1 + (-P1)  # indirect doctest
            (1)

        ::

            sage: H2 = HyperellipticCurve(x^5 + 12*x^4 + 13*x^3 + 15*x^2 + 33*x, x)
            sage: J2 = H2.jacobian()(GF(37))
            sage: P2 = J2(H2.lift_x(2)); P2
            (x + 35, y + 24)
            sage: - P2  # indirect doctest
            (x + 35, y + 15)
            sage: P2 + (-P2)  # indirect doctest
            (1)

        TESTS:

        The following was fixed in :issue:`14264`::

            sage: # needs sage.rings.number_field
            sage: P.<x> = QQ[]
            sage: f = x^5 - x + 1; h = x
            sage: C = HyperellipticCurve(f, h, 'u,v')
            sage: J = C.jacobian()
            sage: K.<t> = NumberField(x^2 - 2)
            sage: R.<x> = K[]
            sage: Q = J(K)([x^2 - t, R(1)])
            sage: Q
            (u^2 - t, v - 1)
            sage: -Q
            (u^2 - t, v + u + 1)
            sage: Q + (-Q)  # indirect doctest
            (1)
        """
        if self.is_zero():
            return self
        polys = self.__polys
        X = self.parent()
        f, h = X.curve().hyperelliptic_polynomials()
        if h.is_zero():
            D = (polys[0],-polys[1])
        else:
            # It is essential that the modulus polys[0] can be converted into
            # the parent of the dividend h. This is not always automatically
            # the case (h can be a rational polynomial and polys[0] can a
            # non-constant polynomial over a number field different from
            # QQ). Hence, we force coercion into a common parent before
            # computing the modulus. See trac #14249
            D = (polys[0],-polys[1]-(h+polys[0]) % (polys[0]))
        return JacobianMorphism_divisor_class_field(X, D, check=False)

    def _add_(self,other):
        r"""
        Return a Mumford representative of the divisor ``self + other``.

        EXAMPLES::

            sage: x = GF(37)['x'].gen()
            sage: H = HyperellipticCurve(x^5 + 12*x^4 + 13*x^3 + 15*x^2 + 33*x)
            sage: J = H.jacobian()(GF(37))

        ::

            sage: P1 = J(H.lift_x(2)); P1
            (x + 35, y + 26)
            sage: P1 + P1  # indirect doctest
            (x^2 + 33*x + 4, y + 13*x)
        """
        X = self.parent()
        C = X.curve()
        f, h = C.hyperelliptic_polynomials()
        genus = C.genus()
        if h == 0:
            D = cantor_composition_simple(self.__polys, other.__polys, f, genus)
            if D[0].degree() > genus:
                D = cantor_reduction_simple(D[0], D[1], f, genus)
        else:
            D = cantor_composition(self.__polys, other.__polys, f, h, genus)
            if D[0].degree() > genus:
                D = cantor_reduction(D[0], D[1], f, h, genus)
        return JacobianMorphism_divisor_class_field(X, D, check=False)

    def _sub_(self, other):
        r"""
        Return a Mumford representative of the divisor ``self - other``.

        EXAMPLES::

            sage: x = GF(37)['x'].gen()
            sage: H = HyperellipticCurve(x^5 + 12*x^4 + 13*x^3 + 15*x^2 + 33*x)
            sage: J = H.jacobian()(GF(37))

        ::

            sage: P1 = J(H.lift_x(2)); P1
            (x + 35, y + 26)
            sage: P1 - P1  # indirect doctest
            (1)

        ::

            sage: P2 = J(H.lift_x(4)); P2
            (x + 33, y + 34)

        Observe that the `x`-coordinates are the same but the
        `y`-coordinates differ::

            sage: P1 - P2  # indirect doctest
            (x^2 + 31*x + 8, y + 7*x + 12)
            sage: P1 + P2  # indirect doctest
            (x^2 + 31*x + 8, y + 4*x + 18)
            sage: (P1 - P2) - (P1 + P2) + 2*P2  # indirect doctest
            (1)
        """
        return self + (-other)
