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
    sage: Jacobian(u^3 + v^3 + w^3)
    Elliptic Curve defined by y^2 = x^3 - 27/4 over Rational Field
    sage: Jacobian(u^4 + v^4 + w^2)
    Elliptic Curve defined by y^2 = x^3 - 4*x over Rational Field

    sage: C = Curve(u^3 + v^3 + w^3)
    sage: Jacobian(C)
    Elliptic Curve defined by y^2 = x^3 - 27/4 over Rational Field

    sage: P2.<u,v,w> = ProjectiveSpace(2, QQ)
    sage: C = P2.subscheme(u^3 + v^3 + w^3)
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
      has a Jacobian elliptic curve

    - ``kwds`` -- optional keyword arguments

    The input ``X`` can be one of the following:

    * A polynomial, see :func:`Jacobian_of_equation` for details.

    * A curve, see :func:`Jacobian_of_curve` for details.

    EXAMPLES::

        sage: R.<u,v,w> = QQ[]
        sage: Jacobian(u^3 + v^3 + w^3)
        Elliptic Curve defined by y^2 = x^3 - 27/4 over Rational Field

        sage: C = Curve(u^3 + v^3 + w^3)
        sage: Jacobian(C)
        Elliptic Curve defined by y^2 = x^3 - 27/4 over Rational Field

        sage: P2.<u,v,w> = ProjectiveSpace(2, QQ)
        sage: C = P2.subscheme(u^3 + v^3 + w^3)
        sage: Jacobian(C)
        Elliptic Curve defined by y^2 = x^3 - 27/4 over Rational Field

        sage: Jacobian(C, morphism=True)
        Scheme morphism:
          From: Closed subscheme of Projective Space of dimension 2 over Rational Field defined by:
          u^3 + v^3 + w^3
          To:   Elliptic Curve defined by y^2 = x^3 - 27/4 over Rational Field
          Defn: Defined on coordinates by sending (u : v : w) to
                (v : -3/2*u + 3/2*w : -1/3*u - 1/3*w)

    TESTS:

    Check that the following doesn't segmentation fault
    (the error message could be improved)::

        sage: Jacobian(GF(11)['x,y'](3))
        Traceback (most recent call last):
        ...
        AssertionError
    """
    try:
        return X.jacobian(**kwds)
    except (AttributeError, TypeError):
        pass

    morphism = kwds.pop('morphism', False)
    from sage.rings.polynomial.multi_polynomial import MPolynomial
    if isinstance(X, MPolynomial):
        if morphism:
            from sage.schemes.curves.constructor import Curve
            return Jacobian_of_equation(X, curve=Curve(X), **kwds)
        return Jacobian_of_equation(X, **kwds)

    from sage.schemes.generic.scheme import Scheme
    if isinstance(X, Scheme) and X.dimension() == 1:
        return Jacobian_of_curve(X, morphism=morphism, **kwds)


def Jacobian_of_curve(curve, morphism=False):
    """
    Return the Jacobian of a genus-one curve.

    INPUT:

    - ``curve`` -- a one-dimensional algebraic variety of genus one

    OUTPUT: its Jacobian elliptic curve

    EXAMPLES::

        sage: R.<u,v,w> = QQ[]
        sage: C = Curve(u^3 + v^3 + w^3)
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
        return Jacobian_of_equation(eqn)
    raise NotImplementedError('Jacobian for this curve is not implemented')


def _plane_cubic_jacobian_morphism(polynomial, curve, jacobian):
    """
    Return an isomorphism from a plane cubic to its Jacobian, if available.

    This is only applicable to homogeneous ternary cubics with a rational
    flex.  If no such isomorphism can be constructed, ``None`` is returned
    so that callers can fall back to the toric multicover.
    """
    try:
        if (polynomial.parent().ngens() != 3 or
                polynomial.nvariables() != 3 or
                not polynomial.is_homogeneous() or
                polynomial.total_degree() != 3):
            return None
    except AttributeError:
        return None

    from sage.schemes.elliptic_curves.constructor import EllipticCurve_from_cubic
    from sage.schemes.elliptic_curves.weierstrass_morphism import baseWI
    from sage.schemes.elliptic_curves.weierstrass_transform import (
        WeierstrassTransformationWithInverse,
    )

    try:
        cubic_map = EllipticCurve_from_cubic(polynomial, morphism=True)
        iso = cubic_map.codomain().isomorphism_to(jacobian)
        # ``WeierstrassIsomorphism.__call__`` is overridden to act on curves
        # and points; we deliberately invoke the underlying ``baseWI.__call__``
        # so that ``iso`` acts on a list of polynomials via the formal
        # Weierstrass change of variables (u, r, s, t).
        fwd_polys = baseWI.__call__(iso, list(cubic_map.defining_polynomials()))
        inv_map = cubic_map.inverse()
        inv_iso_polys = baseWI.__call__(
            ~iso, list(jacobian.ambient_space().coordinate_ring().gens()))
        inv_polys = [p(inv_iso_polys) for p in inv_map.defining_polynomials()]

        # Both quotients below are nonzero scalars by a degree/irreducibility
        # argument, but guard against any unforeseen failure and fall back.
        fwd_post = 1 / (jacobian.defining_polynomial()(fwd_polys) / polynomial)
        inv_post = 1 / (polynomial(inv_polys) / jacobian.defining_polynomial())
    except (ArithmeticError, AttributeError, NotImplementedError,
            TypeError, ValueError, ZeroDivisionError):
        return None

    return WeierstrassTransformationWithInverse(
        curve, jacobian, fwd_polys, fwd_post, inv_polys, inv_post)


def Jacobian_of_equation(polynomial, variables=None, curve=None):
    r"""
    Construct the Jacobian of a genus-one curve given by a polynomial.

    INPUT:

    - ``F`` -- a polynomial defining a plane curve of genus one. May
      be homogeneous or inhomogeneous

    - ``variables`` -- list of two or three variables or ``None``
      (default). The inhomogeneous or homogeneous coordinates. By
      default, all variables in the polynomial are used.

    - ``curve`` -- the genus-one curve defined by ``polynomial`` or
      ``None`` (default). If specified, a suitable morphism from the
      curve to the jacobian elliptic curve is returned.

    OUTPUT:

    An elliptic curve in short Weierstrass form isomorphic to the
    curve ``polynomial=0``. If the optional argument ``curve`` is
    specified, a rational morphism from the genus-one curve to the
    Jacobian elliptic curve is returned.

    EXAMPLES::

        sage: R.<a,b,c> = QQ[]
        sage: f = a^3 + b^3 + 60*c^3
        sage: Jacobian(f)
        Elliptic Curve defined by y^2 = x^3 - 24300 over Rational Field
        sage: Jacobian(f.subs(c=1))
        Elliptic Curve defined by y^2 = x^3 - 24300 over Rational Field

    If we specify the domain curve, an isomorphism is returned when
    one can be constructed over the base field::

        sage: h = Jacobian(f, curve=Curve(f));  h
        Scheme morphism:
          From: Projective Plane Curve over Rational Field defined by a^3 + b^3 + 60*c^3
          To:   Elliptic Curve defined by y^2 = x^3 - 24300 over Rational Field
          Defn: Defined on coordinates by sending (a : b : c) to
                (-c : -3/2*a + 3/2*b : 1/180*a + 1/180*b)

        sage: h([1,-1,0])
        (0 : 1 : 0)

    Plugging in the polynomials defining `h` allows us to verify that
    it is indeed a rational morphism to the elliptic curve::

        sage: E = h.codomain()
        sage: E.defining_polynomial()(h.defining_polynomials()).factor()
        (1/60) * (a^3 + b^3 + 60*c^3)

    By specifying the variables, we can also construct an elliptic
    curve over a polynomial ring::

        sage: R.<u,v,t> = QQ[]
        sage: Jacobian(u^3 + v^3 + t, variables=[u,v])
        Elliptic Curve defined by y^2 = x^3 + (-27/4*t^2) over
         Multivariate Polynomial Ring in u, v, t over Rational Field

    TESTS::

        sage: from sage.schemes.elliptic_curves.jacobian import Jacobian_of_equation
        sage: Jacobian_of_equation(f, variables=[a,b,c])
        Elliptic Curve defined by y^2 = x^3 - 24300 over Rational Field

    Check that the returned morphism is an isomorphism, not the degree
    `9` toric multicover, for a plane cubic with a rational flex::

        sage: P2.<u,v,w> = ProjectiveSpace(2, QQ)
        sage: C = P2.subscheme(u^3 + v^3 + w^3)
        sage: phi = Jacobian(C, morphism=True)
        sage: Ps = [C([-1,0,1]), C([0,-1,1]), C([-1,1,0])]
        sage: [phi(P) for P in Ps]
        [(0 : 1 : 0), (3 : -9/2 : 1), (3 : 9/2 : 1)]
        sage: [phi.inverse()(phi(P)) for P in Ps] == Ps
        True

    The same isomorphism is returned for an asymmetric cubic, where the
    post-rescaling is a non-trivial constant::

        sage: R.<a,b,c> = QQ[]
        sage: g = a^3 + b^3 + 7*c^3
        sage: hg = Jacobian(g, curve=Curve(g))
        sage: Eg = hg.codomain()
        sage: Eg.defining_polynomial()(hg.defining_polynomials()).factor()
        (1/7) * (a^3 + b^3 + 7*c^3)
        sage: hg.inverse()(hg(Curve(g)([1,-1,0])))
        (-1 : 1 : 0)
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
    plane_cubic_morphism = _plane_cubic_jacobian_morphism(polynomial, curve, E)
    if plane_cubic_morphism is not None:
        return plane_cubic_morphism
    X, Y, Z = WeierstrassForm(polynomial, variables=variables, transformation=True)
    from sage.schemes.elliptic_curves.weierstrass_transform import WeierstrassTransformation
    return WeierstrassTransformation(curve, E, [X*Z, Y, Z**3], 1)
