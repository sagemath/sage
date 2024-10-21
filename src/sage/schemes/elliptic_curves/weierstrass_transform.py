# sage_setup: distribution = sagemath-schemes
"""
Morphism to bring a genus-one curve into Weierstrass form

You should use
:func:`~sage.schemes.elliptic_curves.constructor.EllipticCurve_from_cubic`
or
:func:`~sage.schemes.elliptic_curves.constructor.EllipticCurve_from_curve`
to construct the transformation starting with a cubic or with a genus
one curve.

EXAMPLES::

    sage: R.<u,v,w> = QQ[]
    sage: f = EllipticCurve_from_cubic(u^3 + v^3 + w^3, [1,-1,0], morphism=True);  f
    Scheme morphism:
      From: Projective Plane Curve over Rational Field defined by u^3 + v^3 + w^3
      To:   Elliptic Curve defined by y^2 - 9*y = x^3 - 27 over Rational Field
      Defn: Defined on coordinates by sending (u : v : w) to
            (-w : 3*u : 1/3*u + 1/3*v)

    sage: finv = f.inverse();  finv
    Scheme morphism:
      From: Elliptic Curve defined by y^2 - 9*y = x^3 - 27 over Rational Field
      To:   Projective Plane Curve over Rational Field defined by u^3 + v^3 + w^3
      Defn: Defined on coordinates by sending (x : y : z) to
            (1/3*y : -1/3*y + 3*z : -x)

    sage: (u^3 + v^3 + w^3)(f.inverse().defining_polynomials()) * f.inverse().post_rescaling()
    -x^3 + y^2*z - 9*y*z^2 + 27*z^3

    sage: E = finv.domain()
    sage: E.defining_polynomial()(f.defining_polynomials()) * f.post_rescaling()
    u^3 + v^3 + w^3

    sage: f([1,-1,0])
    (0 : 1 : 0)
    sage: f([1,0,-1])
    (3 : 9 : 1)
    sage: f([0,1,-1])
    (3 : 0 : 1)
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

from sage.schemes.generic.morphism import SchemeMorphism_polynomial


class WeierstrassTransformation(SchemeMorphism_polynomial):

    def __init__(self, domain, codomain, defining_polynomials, post_multiplication):
        r"""
        A morphism of a genus-one curve to/from the Weierstrass form.

        INPUT:

        - ``domain``, ``codomain`` -- two schemes, one of which is an
          elliptic curve

        - ``defining_polynomials`` -- triplet of polynomials that
          define the transformation

        - ``post_multiplication`` -- a polynomial to homogeneously
          rescale after substituting the defining polynomials

        EXAMPLES::

            sage: P2.<u,v,w> = ProjectiveSpace(2,QQ)
            sage: C = P2.subscheme(u^3 + v^3 + w^3)
            sage: E = EllipticCurve([2, -1, -1/3, 1/3, -1/27])
            sage: from sage.schemes.elliptic_curves.weierstrass_transform import WeierstrassTransformation
            sage: f = WeierstrassTransformation(C, E, [w, -v-w, -3*u-3*v], 1);  f
            Scheme morphism:
              From: Closed subscheme of Projective Space of dimension 2 over Rational Field defined by:
              u^3 + v^3 + w^3
              To:   Elliptic Curve defined by y^2 + 2*x*y - 1/3*y = x^3 - x^2 + 1/3*x - 1/27
                    over Rational Field
              Defn: Defined on coordinates by sending (u : v : w) to
                    (w : -v - w : -3*u - 3*v)

            sage: f([-1, 1, 0])
            (0 : 1 : 0)
            sage: f([-1, 0, 1])
            (1/3 : -1/3 : 1)
            sage: f([ 0,-1, 1])
            (1/3 : 0 : 1)

            sage: A2.<a,b> = AffineSpace(2,QQ)
            sage: C = A2.subscheme(a^3 + b^3 + 1)
            sage: f = WeierstrassTransformation(C, E, [1, -b-1, -3*a-3*b], 1);  f
            Scheme morphism:
              From: Closed subscheme of Affine Space of dimension 2 over Rational Field defined by:
              a^3 + b^3 + 1
              To:   Elliptic Curve defined by y^2 + 2*x*y - 1/3*y
                    = x^3 - x^2 + 1/3*x - 1/27 over Rational Field
              Defn: Defined on coordinates by sending (a, b) to
                    (1 : -b - 1 : -3*a - 3*b)
            sage: f([-1,0])
            (1/3 : -1/3 : 1)
            sage: f([0,-1])
            (1/3 : 0 : 1)
        """
        Hom = domain.Hom(codomain)
        super().__init__(Hom, defining_polynomials)
        self._post = post_multiplication

    def post_rescaling(self):
        """
        Return the homogeneous rescaling to apply after the coordinate
        substitution.

        OUTPUT: a polynomial; see the example below

        EXAMPLES::

            sage: R.<a,b,c> = QQ[]
            sage: cubic =  a^3+7*b^3+64*c^3
            sage: P = [2,2,-1]
            sage: f = EllipticCurve_from_cubic(cubic, P, morphism=True).inverse()
            sage: f.post_rescaling()
            -1/7

        So here is what it does. If we just plug in the coordinate
        transformation, we get the defining polynomial up to
        scale. This method returns the overall rescaling of the
        equation to bring the result into the standard form::

            sage: cubic(f.defining_polynomials())
            7*x^3 - 7*y^2*z + 1806336*y*z^2 - 155373797376*z^3
            sage: cubic(f.defining_polynomials()) * f.post_rescaling()
            -x^3 + y^2*z - 258048*y*z^2 + 22196256768*z^3
        """
        return self._post


def WeierstrassTransformationWithInverse(domain, codomain,
                                         defining_polynomials, post_multiplication,
                                         inv_defining_polynomials, inv_post_multiplication):
    """
    Construct morphism of a genus-one curve to/from the Weierstrass
    form with its inverse.

    EXAMPLES::

        sage: R.<u,v,w> = QQ[]
        sage: f = EllipticCurve_from_cubic(u^3 + v^3 + w^3, [1,-1,0], morphism=True);  f
        Scheme morphism:
          From: Projective Plane Curve over Rational Field defined by u^3 + v^3 + w^3
          To:   Elliptic Curve defined by y^2 - 9*y = x^3 - 27 over Rational Field
          Defn: Defined on coordinates by sending (u : v : w) to
            (-w : 3*u : 1/3*u + 1/3*v)

        Scheme morphism:
          From: Closed subscheme of Projective Space of dimension 2 over Rational Field defined by:
          u^3 + v^3 + w^3
          To:   Elliptic Curve defined by y^2 + 2*x*y + 1/3*y
                = x^3 - x^2 - 1/3*x - 1/27 over Rational Field
          Defn: Defined on coordinates by sending (u : v : w) to
                (-w : -v + w : 3*u + 3*v)
    """
    fwd = WeierstrassTransformationWithInverse_class(
        domain, codomain, defining_polynomials, post_multiplication)
    inv = WeierstrassTransformationWithInverse_class(
        codomain, domain, inv_defining_polynomials, inv_post_multiplication)
    fwd._inverse = inv
    inv._inverse = fwd
    return fwd


class WeierstrassTransformationWithInverse_class(WeierstrassTransformation):

    def inverse(self):
        """
        Return the inverse.

        OUTPUT:

        A morphism in the opposite direction. This may be a rational
        inverse or an analytic inverse.

        EXAMPLES::

            sage: R.<u,v,w> = QQ[]
            sage: f = EllipticCurve_from_cubic(u^3 + v^3 + w^3, [1,-1,0], morphism=True)
            sage: f.inverse()
            Scheme morphism:
              From: Elliptic Curve defined by y^2 - 9*y = x^3 - 27 over Rational Field
              To:   Projective Plane Curve over Rational Field defined by u^3 + v^3 + w^3
              Defn: Defined on coordinates by sending (x : y : z) to
              (1/3*y : -1/3*y + 3*z : -x)
        """
        return self._inverse
