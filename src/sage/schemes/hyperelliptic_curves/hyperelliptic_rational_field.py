"""
Hyperelliptic curves over the rationals
"""
#*****************************************************************************
#  Copyright (C) 2006 David Kohel <kohel@maths.usyd.edu>
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import sage.rings.abc

from sage.rings.padics.all import pAdicField

from sage.schemes.curves.projective_curve import ProjectivePlaneCurve_field

from . import hyperelliptic_generic


class HyperellipticCurve_rational_field(hyperelliptic_generic.HyperellipticCurve_generic,
                                        ProjectivePlaneCurve_field):

    def matrix_of_frobenius(self, p, prec=20):

        # BUG: should get this method from HyperellipticCurve_generic
        def my_chage_ring(self, R):
            from .constructor import HyperellipticCurve
            f, h = self._hyperelliptic_polynomials
            y = self._printing_ring.gen()
            x = self._printing_ring.base_ring().gen()
            return HyperellipticCurve(f.change_ring(R), h, "%s,%s"%(x,y))

        import sage.schemes.hyperelliptic_curves.monsky_washnitzer as monsky_washnitzer
        if isinstance(p, (sage.rings.abc.pAdicField, sage.rings.abc.pAdicRing)):
            K = p
        else:
            K = pAdicField(p, prec)
        frob_p, forms = monsky_washnitzer.matrix_of_frobenius_hyperelliptic(my_chage_ring(self, K))
        return frob_p

