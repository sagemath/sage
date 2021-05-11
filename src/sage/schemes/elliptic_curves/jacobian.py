# -*- coding: utf-8 -*-
"""
Construct elliptic curves as Jacobians

An elliptic curve is a genus one curve with a designated point. The
Jacobian of a genus-one curve can be defined as the set of line
bundles on the curve, and is isomorphic to the original genus-one
curve. It is also an elliptic curve with the trivial line bundle as
designated point. The utility of this construction is that we can
construct elliptic curves without having to specify which point we
take as the origin.

EXAMPLES::

    sage: R.<u,v,w> = QQ[]
    sage: Jacobian(u^3+v^3+w^3)
    Elliptic Curve defined by y^2 = x^3 - 27/4 over Rational Field
    sage: Jacobian(u^4+v^4+w^2)
    Elliptic Curve defined by y^2 = x^3 - 4*x over Rational Field

    sage: C = Curve(u^3+v^3+w^3)
    sage: Jacobian(C)
    Elliptic Curve defined by y^2 = x^3 - 27/4 over Rational Field

    sage: P2.<u,v,w> = ProjectiveSpace(2, QQ)
    sage: C = P2.subscheme(u^3+v^3+w^3)
    sage: Jacobian(C)
    Elliptic Curve defined by y^2 = x^3 - 27/4 over Rational Field

One can also define Jacobians of varieties that are not genus-one
curves. These are not implemented in this module, but we call the
relevant functionality::

    sage: R.<x> = PolynomialRing(QQ)
    sage: f = x**5 + 1184*x**3 + 1846*x**2 + 956*x + 560
    sage: C = HyperellipticCurve(f)
    sage: Jacobian(C)
    Jacobian of Hyperelliptic Curve over Rational Field defined
    by y^2 = x^5 + 1184*x^3 + 1846*x^2 + 956*x + 560

REFERENCES:

- :wikipedia:`Jacobian_variety`
"""

##############################################################################
#       Copyright (C) 2013 Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
##############################################################################

from sage.schemes.elliptic_curves.constructor import EllipticCurve


def Jacobian(X, **kwds):
    """
    Return the Jacobian.

    INPUT:

    - ``X`` -- polynomial, algebraic variety, or anything else that
      has a Jacobian elliptic curve.

    - ``kwds`` -- optional keyword arguments.

    The input ``X`` can be one of the following:

    * A polynomial, see :func:`Jacobian_of_equation` for details.

    * A curve, see :func:`Jacobian_of_curve` for details.

    EXAMPLES::

        sage: R.<u,v,w> = QQ[]
        sage: Jacobian(u^3+v^3+w^3)
        Elliptic Curve defined by y^2 = x^3 - 27/4 over Rational Field

        sage: C = Curve(u^3+v^3+w^3)
        sage: Jacobian(C)
        Elliptic Curve defined by y^2 = x^3 - 27/4 over Rational Field

        sage: P2.<u,v,w> = ProjectiveSpace(2, QQ)
        sage: C = P2.subscheme(u^3+v^3+w^3)
        sage: Jacobian(C)
        Elliptic Curve defined by y^2 = x^3 - 27/4 over Rational Field

        sage: Jacobian(C, morphism=True)
        Scheme morphism:
          From: Closed subscheme of Projective Space of dimension 2 over Rational Field defined by:
          u^3 + v^3 + w^3
          To:   Elliptic Curve defined by y^2 = x^3 - 27/4 over Rational Field
          Defn: Defined on coordinates by sending (u : v : w) to
                (-u^4*v^4*w - u^4*v*w^4 - u*v^4*w^4 :
                1/2*u^6*v^3 - 1/2*u^3*v^6 - 1/2*u^6*w^3 + 1/2*v^6*w^3 + 1/2*u^3*w^6 - 1/2*v^3*w^6 :
                u^3*v^3*w^3)
    """
    try:
        return X.jacobian(**kwds)
    except AttributeError:
        pass

    morphism = kwds.pop('morphism', False)
    from sage.rings.polynomial.multi_polynomial_element import is_MPolynomial
    if is_MPolynomial(X):
        if morphism:
            from sage.schemes.curves.constructor import Curve
            return Jacobian_of_equation(X, curve=Curve(X), **kwds)
        else:
            return Jacobian_of_equation(X, **kwds)

    from sage.schemes.generic.scheme import is_Scheme
    if is_Scheme(X) and X.dimension() == 1:
        return Jacobian_of_curve(X, morphism=morphism, **kwds)


def Jacobian_of_curve(curve, morphism=False):
    """
    Return the Jacobian of a genus-one curve

    INPUT:

    - ``curve`` -- a one-dimensional algebraic variety of genus one.

    OUTPUT:

    Its Jacobian elliptic curve.

    EXAMPLES::

        sage: R.<u,v,w> = QQ[]
        sage: C = Curve(u^3+v^3+w^3)
        sage: Jacobian(C)
        Elliptic Curve defined by y^2 = x^3 - 27/4 over Rational Field
    """
    eqn = None
    try:
        eqn = curve.defining_polynomial()
    except AttributeError:
        pass
    if len(curve.defining_polynomials()) == 1:
        eqn = curve.defining_polynomials()[0]
    if eqn is not None:
        if morphism:
            return Jacobian_of_equation(eqn, curve=curve)
        else:
            return Jacobian_of_equation(eqn)
    raise NotImplementedError('Jacobian for this curve is not implemented')


def Jacobian_of_equation(polynomial, variables=None, curve=None):
    r"""
    Construct the Jacobian of a genus-one curve given by a polynomial.

    INPUT:

    - ``F`` -- a polynomial defining a plane curve of genus one. May
      be homogeneous or inhomogeneous.

    - ``variables`` -- list of two or three variables or ``None``
      (default). The inhomogeneous or homogeneous coordinates. By
      default, all variables in the polynomial are used.

    - ``curve`` -- the genus-one curve defined by ``polynomial`` or #
      ``None`` (default). If specified, suitable morphism from the
      jacobian elliptic curve to the curve is returned.

    OUTPUT:

    An elliptic curve in short Weierstrass form isomorphic to the
    curve ``polynomial=0``. If the optional argument ``curve`` is
    specified, a rational multicover from the Jacobian elliptic curve
    to the genus-one curve is returned.

    EXAMPLES::

        sage: R.<a,b,c> = QQ[]
        sage: f = a^3+b^3+60*c^3
        sage: Jacobian(f)
        Elliptic Curve defined by y^2 = x^3 - 24300 over Rational Field
        sage: Jacobian(f.subs(c=1))
        Elliptic Curve defined by y^2 = x^3 - 24300 over Rational Field

    If we specify the domain curve the birational covering is returned::

        sage: h = Jacobian(f, curve=Curve(f));  h
        Scheme morphism:
          From: Projective Plane Curve over Rational Field defined by a^3 + b^3 + 60*c^3
          To:   Elliptic Curve defined by y^2 = x^3 - 24300 over Rational Field
          Defn: Defined on coordinates by sending (a : b : c) to
                (-216000*a^4*b^4*c - 12960000*a^4*b*c^4 - 12960000*a*b^4*c^4 :
                108000*a^6*b^3 - 108000*a^3*b^6 - 6480000*a^6*c^3 + 6480000*b^6*c^3 + 388800000*a^3*c^6 - 388800000*b^3*c^6 :
                216000*a^3*b^3*c^3)

        sage: h([1,-1,0])
        (0 : 1 : 0)

    Plugging in the polynomials defining `h` allows us to verify that
    it is indeed a rational morphism to the elliptic curve::

        sage: E = h.codomain()
        sage: E.defining_polynomial()(h.defining_polynomials()).factor()
        (2519424000000000) * c^3 * b^3 * a^3 * (a^3 + b^3 + 60*c^3) *
        (a^9*b^6 + a^6*b^9 - 120*a^9*b^3*c^3 + 900*a^6*b^6*c^3 - 120*a^3*b^9*c^3 +
        3600*a^9*c^6 + 54000*a^6*b^3*c^6 + 54000*a^3*b^6*c^6 + 3600*b^9*c^6 +
        216000*a^6*c^9 - 432000*a^3*b^3*c^9 + 216000*b^6*c^9)


    By specifying the variables, we can also construct an elliptic
    curve over a polynomial ring::

        sage: R.<u,v,t> = QQ[]
        sage: Jacobian(u^3+v^3+t, variables=[u,v])
        Elliptic Curve defined by y^2 = x^3 + (-27/4*t^2) over
        Multivariate Polynomial Ring in u, v, t over Rational Field

    TESTS::

        sage: from sage.schemes.elliptic_curves.jacobian import Jacobian_of_equation
        sage: Jacobian_of_equation(f, variables=[a,b,c])
        Elliptic Curve defined by y^2 = x^3 - 24300 over Rational Field
    """
    from sage.schemes.toric.weierstrass import WeierstrassForm
    f, g = WeierstrassForm(polynomial, variables=variables)
    try:
        K = polynomial.base_ring()
        f = K(f)
        g = K(g)
    except (TypeError, ValueError):
        pass
    E = EllipticCurve([f, g])
    if curve is None:
        return E
    X, Y, Z = WeierstrassForm(polynomial, variables=variables, transformation=True)
    from sage.schemes.elliptic_curves.weierstrass_transform import WeierstrassTransformation
    return WeierstrassTransformation(curve, E, [X*Z, Y, Z**3], 1)
