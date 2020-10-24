"""
Kummer surfaces over a general ring
"""

#*****************************************************************************
#  Copyright (C) 2006 David Kohel <kohel@maths.usyd.edu>
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.schemes.projective.projective_space import ProjectiveSpace
from sage.schemes.projective.projective_subscheme\
     import AlgebraicScheme_subscheme_projective
from sage.categories.homset import Hom
from sage.categories.all import Schemes

# The generic genus 2 curve in Weierstrass form:
#
# y^2 + (c9*x^3 + c6*x^2 + c3*x + c0)*y =
#       a12*x^6 + a10*x^5 + a8*x^4 + a6*x^3 + a4*x^2 + a2*x + a0.
#
# Transforms to:
#
# y^2 = (c9^2 + 4*a12)*x^6 + (2*c6*c9 + 4*a10)*x^5
#     + (2*c3*c9 + c6^2 + 4*a8)*x^4 + (2*c0*c9 + 2*c3*c6 + 4*a6)*x^3
#     + (2*c0*c6 + c3^2 + 4*a4)*x^2 + (2*c0*c3 + 4*a2)*x + c0^2 + 4*a0

class KummerSurface(AlgebraicScheme_subscheme_projective):
    def __init__(self, J):
        """
        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: f = x^5 + x + 1
            sage: X = HyperellipticCurve(f)
            sage: J = Jacobian(X)
            sage: K = KummerSurface(J); K
            Closed subscheme of Projective Space of dimension 3 over Rational Field defined by:
            X0^4 - 4*X0*X1^3 + 4*X0^2*X1*X2 - 4*X0*X1^2*X2 + 2*X0^2*X2^2 + X2^4 - 4*X0^3*X3 - 2*X0^2*X1*X3 - 2*X1*X2^2*X3 + X1^2*X3^2 - 4*X0*X2*X3^2
        """
        R = J.base_ring()
        PP = ProjectiveSpace(3, R, ["X0", "X1", "X2", "X3"])
        X0, X1, X2, X3 = PP.gens()
        C = J.curve()
        f, h = C.hyperelliptic_polynomials()
        a12 = f[0]
        a10 = f[1]
        a8 = f[2]
        a6 = f[3]
        a4 = f[4]
        a2 = f[5]
        a0 = f[6]
        if h != 0:
            c6 = h[0]
            c4 = h[1]
            c2 = h[2]
            c0 = h[3]
            a12, a10, a8, a6, a4, a2, a0 = \
                 (4*a12 + c6**2,
                  4*a10 + 2*c4*c6,
                  4*a8 + 2*c2*c6 + c4**2,
                  4*a6 + 2*c0*c6 + 2*c2*c4,
                  4*a4 + 2*c0*c4 + c2**2,
                  4*a2 + 2*c0*c2,
                  4*a0 + c0**2)
        F = \
          (-4*a8*a12 + a10**2)*X0**4 + \
          -4*a6*a12*X0**3*X1 + \
          -2*a6*a10*X0**3*X2 + \
          -4*a12*X0**3*X3 + \
          -4*a4*a12*X0**2*X1**2 + \
          (4*a2*a12 - 4*a4*a10)*X0**2*X1*X2 + \
          -2*a10*X0**2*X1*X3 + \
          (-4*a0*a12 + 2*a2*a10 - 4*a4*a8 + a6**2)*X0**2*X2**2 + \
          -4*a8*X0**2*X2*X3 + \
          -4*a2*a12*X0*X1**3 + \
          (8*a0*a12 - 4*a2*a10)*X0*X1**2*X2 + \
          (4*a0*a10 - 4*a2*a8)*X0*X1*X2**2 + \
          -2*a6*X0*X1*X2*X3 + \
          -2*a2*a6*X0*X2**3 + \
          -4*a4*X0*X2**2*X3 + \
          -4*X0*X2*X3**2 + \
          -4*a0*a12*X1**4 + \
          -4*a0*a10*X1**3*X2 + \
          -4*a0*a8*X1**2*X2**2 + \
          X1**2*X3**2 + \
          -4*a0*a6*X1*X2**3 + \
          -2*a2*X1*X2**2*X3 + \
          (-4*a0*a4 + a2**2)*X2**4 + \
          -4*a0*X2**3*X3
        AlgebraicScheme_subscheme_projective.__init__(self, PP, F)
        X, Y, Z = C.ambient_space().gens()
        if a0 == 0:
            a0 = a2
        phi = Hom(C, self)([X.parent().zero(), Z**2, X*Z, a0*X**2], Schemes())
        C._kummer_morphism = phi
        J._kummer_surface = self
