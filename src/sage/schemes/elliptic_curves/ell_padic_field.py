"""
Elliptic curves over padic fields
"""
# ****************************************************************************
#       Copyright (C) 2007 Robert Bradshaw <robertwb@math.washington.edu>
#                          William Stein   <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from .ell_field import EllipticCurve_field
from . import ell_point
from sage.rings.all import PolynomialRing

# Elliptic curves are very different than genus > 1 hyperelliptic curves,
# there is an "is a" relationship here, and common implementation with regard
# Coleman integration.
from sage.schemes.hyperelliptic_curves.hyperelliptic_padic_field import HyperellipticCurve_padic_field


class EllipticCurve_padic_field(EllipticCurve_field, HyperellipticCurve_padic_field):
    """
    Elliptic curve over a padic field.

    EXAMPLES::

        sage: Qp=pAdicField(17)
        sage: E=EllipticCurve(Qp,[2,3]); E
        Elliptic Curve defined by y^2  = x^3 + (2+O(17^20))*x + (3+O(17^20)) over 17-adic Field with capped relative precision 20
        sage: E == loads(dumps(E))
        True
    """

    _point = ell_point.EllipticCurvePoint_field

    _genus = 1

    def frobenius(self, P=None):
        """
        Return the Frobenius as a function on the group of points of
        this elliptic curve.

        EXAMPLES::

            sage: Qp = pAdicField(13)
            sage: E = EllipticCurve(Qp,[1,1])
            sage: type(E.frobenius())
            <... 'function'>
            sage: point = E(0,1)
            sage: E.frobenius(point)
            (0 : 1 + O(13^20) : 1 + O(13^20))

        Check that :trac:`29709` is fixed::

            sage: Qp = pAdicField(13)
            sage: E = EllipticCurve(Qp,[0,0,1,0,1])
            sage: E.frobenius(E(1,1))
            Traceback (most recent call last):
            ...
            NotImplementedError: Curve must be in weierstrass normal form.
            sage: E = EllipticCurve(Qp,[0,1,0,0,1])
            sage: E.frobenius(E(0,1))
            (0 : 1 + O(13^20) : 1 + O(13^20))
        """
        try:
            _frob = self._frob
        except AttributeError:
            K = self.base_field()
            p = K.prime()
            x = PolynomialRing(K, 'x').gen(0)

            a1, a2, a3, a4, a6 = self.a_invariants()
            if a1 != 0 or a3 != 0:
                raise NotImplementedError("Curve must be in weierstrass normal form.")

            f = x*x*x + a2*x*x + a4*x + a6
            h = (f(x**p) - f**p)

            # internal function: I don't know how to doctest it...
            def _frob(P):
                x0 = P[0]
                y0 = P[1]
                uN = (1 + h(x0)/y0**(2*p)).sqrt()
                yres = y0**p * uN
                xres = x0**p
                if (yres-y0).valuation() == 0:
                    yres = -yres
                return self.point([xres, yres, K(1)])

            self._frob = _frob

        if P is None:
            return _frob
        else:
            return _frob(P)
