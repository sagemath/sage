"""
Cyclic covers curves over a general ring

EXAMPLES::

    sage: ZZx.<x> = ZZ[]
    sage: C = CyclicCover(5, x^5 + x + 1); C
    Cyclic Cover of P^1 over Integer Ring defined by y^5 = x^5 + x + 1
    sage: C.genus()
    6
    sage: D = C.projective_closure(); D
    Projective Plane Curve over Integer Ring defined by x0^5 + x0^4*x1 + x1^5 - x2^5
    sage: D.change_ring(QQ).genus()
    6
    sage: C.change_ring(GF(5))
    Traceback (most recent call last):
    ...
    ValueError: As the characteristic divides the order of the cover, this model is not smooth.

    sage: GF7x.<x> = GF(7)[]
    sage: C = CyclicCover(3, x^9 + x + 1)
    sage: C
    Cyclic Cover of P^1 over Finite Field of size 7 defined by y^3 = x^9 + x + 1
    sage: C.genus()
    7
    sage: C.projective_closure()
    Traceback (most recent call last):
    ...
    NotImplementedError: Weighted Projective Space is not implemented
"""

# *****************************************************************************
#  Copyright (C) 2018 Edgar Costa <edgarc@mit.edu>
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# *****************************************************************************

from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.structure.category_object import normalize_names
from sage.arith.misc import GCD
from sage.schemes.curves.affine_curve import AffinePlaneCurve


class CyclicCover_generic(AffinePlaneCurve):
    def __init__(self, AA, r, f, names=None):
        """
        Cyclic covers over a general ring.

        INPUT:

        - ``A`` -- ambient affine space

        - ``r`` -- degree of the cover

        - ``f`` -- univariate polynomial

        - ``names`` -- (default: ``["x","y"]``); names for the
          coordinate functions

        TESTS::

            sage: ZZx.<x> = ZZ[]
            sage: C = CyclicCover(5, x^5 + x + 1); C
            Cyclic Cover of P^1 over Integer Ring defined by y^5 = x^5 + x + 1
            sage: C.genus()
            6
            sage: D = C.projective_closure(); D
            Projective Plane Curve over Integer Ring defined by x0^5 + x0^4*x1 + x1^5 - x2^5
            sage: D.change_ring(QQ).genus()
            6
            sage: C.change_ring(GF(5))
            Traceback (most recent call last):
            ...
            ValueError: As the characteristic divides the order of the cover,
            this model is not smooth.


            sage: GF7x.<x> = GF(7)[]
            sage: C = CyclicCover(3, x^9 + x + 1)
            sage: C
            Cyclic Cover of P^1 over Finite Field of size 7 defined by y^3 = x^9 + x + 1
            sage: C.genus()
            7
            sage: C.projective_closure()
            Traceback (most recent call last):
            ...
            NotImplementedError: Weighted Projective Space is not implemented
        """
        x, y = AA.gens()
        self._r = r
        self._d = f.degree()
        self._delta = GCD(self._r, self._d)
        self._genus = ((self._d - 1) * (self._r - 1) - (self._delta - 1)) // 2
        self._f = f

        F = y**r - f(x)
        AffinePlaneCurve.__init__(self, AA, F)
        if names is None:
            names = ("x", "y")
        else:
            names = normalize_names(2, names)
        self._names = names

    def change_ring(self, R):
        """
        Return this CyclicCover over a new base ring R.

        EXAMPLES::

            sage: ZZx.<x> = ZZ[]
            sage: C = CyclicCover(5, x^5 + x + 1)
            sage: C.change_ring(GF(5))
            Traceback (most recent call last):
            ...
            ValueError: As the characteristic divides the order of the cover,
            this model is not smooth.
            sage: C.change_ring(GF(3))
            Traceback (most recent call last):
            ...
            ValueError: Not a smooth Cyclic Cover of P^1: singularity in the
            provided affine patch.
            sage: C.change_ring(GF(17))
            Cyclic Cover of P^1 over Finite Field of size 17 defined by y^5 = x^5 + x + 1
        """
        from .constructor import CyclicCover

        return CyclicCover(self._r, self._f.change_ring(R), self._names)

    base_extend = change_ring

    def _repr_(self):
        """
        String representation of cyclic covers.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: CyclicCover(2, x^5 + x + 1)
            Cyclic Cover of P^1 over Rational Field defined by y^2 = x^5 + x + 1
            sage: CyclicCover(3, x^5 + x + 1)
            Cyclic Cover of P^1 over Rational Field defined by y^3 = x^5 + x + 1
            sage: CyclicCover(5, x^5 + x + 1)
            Cyclic Cover of P^1 over Rational Field defined by y^5 = x^5 + x + 1
            sage: CyclicCover(15, x^9 + x + 1)
            Cyclic Cover of P^1 over Rational Field defined by y^15 = x^9 + x + 1
        """
        R = self.base_ring()
        x, y = self.ambient_space().gens()
        r = self._r
        f = self._f
        return "Cyclic Cover of P^1 over %s defined by %s = %s" % (R, y**r, f(x))

    def __eq__(self, other):
        """
        Test of equality.

        EXAMPLES::

            sage: ZZx.<x> = ZZ[]
            sage: C0 = CyclicCover(5, x^5 + x + 1)
            sage: C1 = C0.change_ring(QQ)
            sage: C1 == C0
            False
            sage: C2 = CyclicCover(3, x^5 + x + 1)
            sage: C2 == C0
            False
            sage: C3 = CyclicCover(5, x^6 + x + 1)
            sage: C3 == C0
            False
            sage: C0 == CyclicCover(5, x^5 + x + 1)
            True
        """
        if not isinstance(other, CyclicCover_generic):
            return False

        return (
            (self.base_ring() == other.base_ring())
            and (self._r == other._r)
            and (self._f == other._f)
        )

    def __ne__(self, other):
        """
        Test of not equality.

        EXAMPLES::

            sage: ZZx.<x> = ZZ[]
            sage: C0 = CyclicCover(5, x^5 + x + 1)
            sage: C1 = C0.change_ring(QQ)
            sage: C1 != C0
            True
            sage: C2 = CyclicCover(3, x^5 + x + 1)
            sage: C2 != C0
            True
            sage: C3 = CyclicCover(5, x^6 + x + 1)
            sage: C3 != C0
            True
            sage: C0 != CyclicCover(5, x^5 + x + 1)
            False
        """
        return not self == other

    def order(self):
        """
        The order of the cover.

        EXAMPLES::

            sage: ZZx.<x> = ZZ[]
            sage: CyclicCover(5, x^5 + x + 1).order()
            5
            sage: CyclicCover(3, x^5 + x + 1).order()
            3
        """
        return self._r

    def genus(self):
        """
        The geometric genus of the curve.

        EXAMPLES::

            sage: ZZx.<x> = ZZ[]
            sage: CyclicCover(5, x^5 + x + 1).genus()
            6
            sage: CyclicCover(3, x^5 + x + 1).genus()
            4
        """
        return self._genus

    def projective_closure(self, **kwds):
        """
        Return the projective closure of this affine curve.

        EXAMPLES::

            sage: GF7x.<x> = GF(7)[]
            sage: CyclicCover(3, x^9 + x + 1).projective_closure()
            Traceback (most recent call last):
            ...
            NotImplementedError: Weighted Projective Space is not implemented

            sage: ZZx.<x> = ZZ[]
            sage: CyclicCover(5, x^5 + x + 1).projective_closure()
            Projective Plane Curve over Integer Ring defined by x0^5 + x0^4*x1 + x1^5 - x2^5
        """
        # test d = 3 and 4
        if self._d == self._r:
            return AffinePlaneCurve.projective_closure(self, **kwds)
        else:
            raise NotImplementedError("Weighted Projective Space is not implemented")

    def cover_polynomial(self, K=None, var='x'):
        """
        Return the polynomial defining the cyclic cover.

        EXAMPLES::

            sage: ZZx.<x> = ZZ[]; CyclicCover(5, x^5 + x + 1).cover_polynomial()
            x^5 + x + 1
        """

        if K is None:
            return self._f
        else:
            P = PolynomialRing(K, var)
            return P(self._f)

    def is_singular(self):
        r"""
        Return if this curve is singular or not.

        This just checks that the characteristic of the ring does not divide the
        order of the cover and that the defining polynomial of the cover is
        square free.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: CyclicCover(3, x^5 + x + 1).is_singular()
            False
            sage: CyclicCover(3, (x^5 + x + 1)^2, check_smooth=False).is_singular()
            True
        """
        P = self._f.parent()
        r = self._r
        if P(r) == 0:
            return True
        else:
            return not self._f.is_squarefree()

    def is_smooth(self):
        r"""
        Return if this curve is smooth or not.

        This just checks that the characteristic of the ring does not divide the
        order of the cover and that the defining polynomial of the cover is
        square free.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: CyclicCover(3, x^5 + x + 1).is_smooth()
            True
            sage: CyclicCover(3, (x^5 + x + 1)^2, check_smooth=False).is_smooth()
            False
        """
        return not self.is_singular()
