# sage_setup: distribution = sagemath-schemes
"""
Jacobian of a hyperelliptic curve of genus 2
"""

#*****************************************************************************
#  Copyright (C) 2006 David Kohel <kohel@maths.usyd.edu>
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from . import jacobian_generic
from . import kummer_surface

# The generic genus 2 curve in Weierstrass form:
#
# y^2 + (c0*x^3 + c2*x^2 + c4*x + c6)*y =
#     a0*x^6 + a2*x^5 + a4*x^4 + a6*x^3 + a8*x^2 + a10*x + a12.
#
# Transforms to:
#
# y^2 = (4*a0 + c0^2)*x^6 + (4*a2 + 2*c0*c2)*x^5
#     + (4*a4 + 2*c0*c4 + c2^2)*x^4 + (4*a6 + 2*c0*c6 + 2*c2*c4)*x^3
#     + (4*a8 + 2*c2*c6 + c4^2)*x^2 + (4*a10 + 2*c4*c6)*x + 4*a12 + c6^2


class HyperellipticJacobian_g2(jacobian_generic.HyperellipticJacobian_generic):
    def kummer_surface(self):
        try:
            return self._kummer_surface
        except AttributeError:
            return kummer_surface.KummerSurface(self)
