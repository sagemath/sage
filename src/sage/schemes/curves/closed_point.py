# sage_setup: distribution = sagemath-schemes
"""
Closed points of integral curves

A rational point of a curve in Sage is represented by its coordinates. If the
curve is defined over finite field and integral, that is reduced and
irreducible, then it is empowered by the global function field machinery of
Sage. Thus closed points of the curve are computable, as represented by maximal
ideals of the coordinate ring of the ambient space.

EXAMPLES::

    sage: F.<a> = GF(2)
    sage: P.<x,y> = AffineSpace(F, 2)
    sage: C = Curve(y^2 + y - x^3)
    sage: C.closed_points()
    [Point (x, y), Point (x, y + 1)]
    sage: C.closed_points(2)
    [Point (y^2 + y + 1, x + 1),
     Point (y^2 + y + 1, x + y),
     Point (y^2 + y + 1, x + y + 1)]
    sage: C.closed_points(3)
    [Point (x^2 + x + y, x*y + 1, y^2 + x + 1),
     Point (x^2 + x + y + 1, x*y + x + 1, y^2 + x)]

Closed points of projective curves are represented by homogeneous maximal
ideals::

    sage: F.<a> = GF(2)
    sage: P.<x,y,z> = ProjectiveSpace(F, 2)
    sage: C = Curve(x^3*y + y^3*z + x*z^3)
    sage: C.closed_points()
    [Point (x, z), Point (x, y), Point (y, z)]
    sage: C.closed_points(2)
    [Point (y^2 + y*z + z^2, x + y + z)]
    sage: C.closed_points(3)
    [Point (y^3 + y^2*z + z^3, x + y),
     Point (y^3 + y*z^2 + z^3, x + z),
     Point (x^2 + x*z + y*z + z^2, x*y + x*z + z^2, y^2 + x*z),
     Point (x^2 + y*z, x*y + x*z + z^2, y^2 + x*z + y*z),
     Point (x^3 + x*z^2 + z^3, y + z),
     Point (x^2 + y*z + z^2, x*y + x*z + y*z, y^2 + x*z + y*z + z^2),
     Point (x^2 + y*z + z^2, x*y + z^2, y^2 + x*z + y*z)]

Rational points are easily converted to closed points and vice versa if the
closed point is of degree one::

    sage: F.<a> = GF(2)
    sage: P.<x,y,z> = ProjectiveSpace(F, 2)
    sage: C = Curve(x^3*y + y^3*z + x*z^3)
    sage: p1, p2, p3 = C.closed_points()
    sage: p1.rational_point()
    (0 : 1 : 0)
    sage: p2.rational_point()
    (0 : 0 : 1)
    sage: p3.rational_point()
    (1 : 0 : 0)
    sage: _.closed_point()
    Point (y, z)
    sage: _ == p3
    True

AUTHORS:

- Kwankyu Lee (2019-03): initial version
"""

# *****************************************************************************
#       Copyright (C) 2019 Kwankyu Lee <kwankyu@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# *****************************************************************************

from sage.structure.richcmp import richcmp
from sage.schemes.generic.point import SchemeTopologicalPoint_prime_ideal


class CurveClosedPoint(SchemeTopologicalPoint_prime_ideal):
    """
    Base class of closed points of curves.
    """
    pass


class IntegralCurveClosedPoint(CurveClosedPoint):
    """
    Closed points of integral curves.

    INPUT:

    - ``curve`` -- the curve to which the closed point belongs

    - ``prime_ideal`` -- a prime ideal

    - ``degree`` -- degree of the closed point

    EXAMPLES::

        sage: # needs sage.rings.finite_rings
        sage: F.<a> = GF(4)
        sage: P.<x,y> = AffineSpace(F, 2)
        sage: C = Curve(y^2 + y - x^3)
        sage: C.closed_points()
        [Point (x, y),
         Point (x, y + 1),
         Point (x + a, y + a),
         Point (x + a, y + (a + 1)),
         Point (x + (a + 1), y + a),
         Point (x + (a + 1), y + (a + 1)),
         Point (x + 1, y + a),
         Point (x + 1, y + (a + 1))]
    """
    def __init__(self, curve, prime_ideal, degree):
        """
        Initialize.

        TESTS::

            sage: # needs sage.rings.finite_rings
            sage: F.<a> = GF(4)
            sage: P.<x,y> = AffineSpace(F, 2)
            sage: C = Curve(y^2 + y - x^3)
            sage: p = C([0,0]); p
            (0, 0)
            sage: loads(dumps(p)) == p
            True
        """
        super().__init__(curve.ambient_space(), prime_ideal)

        self._curve = curve
        self._degree = degree

    def __hash__(self):
        """
        Return the hash of ``self``.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: F.<a> = GF(4)
            sage: P.<x,y> = AffineSpace(F, 2)
            sage: C = Curve(y^2 + y - x^3)
            sage: pts = C.closed_points()
            sage: p = pts[0]
            sage: {p: 1}
            {Point (x, y): 1}
        """
        return hash((self.parent(), self.prime_ideal()))

    def _richcmp_(self, other, op):
        """
        Compare ``self`` and ``other`` with respect to the operator.

        INPUT:

        - ``other`` -- a closed point

        - ``op`` -- a comparison operator

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: F.<a> = GF(4)
            sage: P.<x,y> = AffineSpace(F, 2)
            sage: C = Curve(y^2 + y - x^3)
            sage: pts = C.closed_points()
            sage: pts[0] == pts[1]
            False
        """
        return richcmp((self._curve, self.prime_ideal()), (other._curve, other.prime_ideal()), op)

    def _repr_(self):
        """
        Return the string representation of the closed point.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: F.<a> = GF(4)
            sage: P.<x,y> = AffineSpace(F, 2)
            sage: C = Curve(y^2 + y - x^3)
            sage: pts = C.closed_points()
            sage: pts[0]
            Point (x, y)
        """
        return "Point ({})".format(', '.join(repr(g) for g in self.prime_ideal().gens()))

    def curve(self):
        """
        Return the curve to which this point belongs.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: F.<a> = GF(4)
            sage: P.<x,y> = AffineSpace(F, 2)
            sage: C = Curve(y^2 + y - x^3)
            sage: pts = C.closed_points()
            sage: p = pts[0]
            sage: p.curve()
            Affine Plane Curve over Finite Field in a of size 2^2 defined by x^3 + y^2 + y
        """
        return self._curve

    def degree(self):
        """
        Return the degree of the point.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: F.<a> = GF(4)
            sage: P.<x,y> = AffineSpace(F, 2)
            sage: C = Curve(y^2 + y - x^3)
            sage: pts = C.closed_points()
            sage: p = pts[0]
            sage: p.degree()
            1
        """
        return self._degree

    def places(self):
        """
        Return all places on this closed point.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: F.<a> = GF(4)
            sage: P.<x,y> = AffineSpace(F, 2)
            sage: C = Curve(y^2 + y - x^3)
            sage: pts = C.closed_points()
            sage: p = pts[0]
            sage: p.places()
            [Place (x, y)]
        """
        return self._curve.places_on(self)

    def place(self):
        """
        Return a place on this closed point.

        If there are more than one, arbitrary one is chosen.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: F.<a> = GF(4)
            sage: P.<x,y> = AffineSpace(F, 2)
            sage: C = Curve(y^2 + y - x^3)
            sage: pts = C.closed_points()
            sage: p = pts[0]
            sage: p.place()
            Place (x, y)
        """
        return self._curve.places_on(self)[0]


class IntegralAffineCurveClosedPoint(IntegralCurveClosedPoint):
    """
    Closed points of affine curves.
    """
    def rational_point(self):
        """
        Return the rational point if this closed point is of degree `1`.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: A.<x,y> = AffineSpace(GF(3^2), 2)
            sage: C = Curve(y^2 - x^5 - x^4 - 2*x^3 - 2*x - 2)
            sage: C.closed_points()
            [Point (x, y + (z2 + 1)),
             Point (x, y + (-z2 - 1)),
             Point (x + (z2 + 1), y + (z2 - 1)),
             Point (x + (z2 + 1), y + (-z2 + 1)),
             Point (x - 1, y + (z2 + 1)),
             Point (x - 1, y + (-z2 - 1)),
             Point (x + (-z2 - 1), y + z2),
             Point (x + (-z2 - 1), y + (-z2)),
             Point (x + 1, y + 1),
             Point (x + 1, y - 1)]
            sage: [p.rational_point() for p in _]
            [(0, 2*z2 + 2),
             (0, z2 + 1),
             (2*z2 + 2, 2*z2 + 1),
             (2*z2 + 2, z2 + 2),
             (1, 2*z2 + 2),
             (1, z2 + 1),
             (z2 + 1, 2*z2),
             (z2 + 1, z2),
             (2, 2),
             (2, 1)]
            sage: set(_) == set(C.rational_points())
            True
        """
        if self.degree() != 1:
            raise ValueError("not a rational point")

        G = self.prime_ideal().groebner_basis()
        C = self._curve
        return C([g.reduce(G) for g in C.ambient_space().gens()])

    def projective(self, i=0):
        """
        Return the point in the projective closure of the curve, of which this
        curve is the ``i``-th affine patch.

        INPUT:

        - ``i`` -- integer

        EXAMPLES::

            sage: F.<a> = GF(2)
            sage: A.<x,y> = AffineSpace(F, 2)
            sage: C = Curve(y^2 + y - x^3, A)
            sage: p1, p2 = C.closed_points()
            sage: p1
            Point (x, y)
            sage: p2
            Point (x, y + 1)
            sage: p1.projective()
            Point (x1, x2)
            sage: p2.projective(0)
            Point (x1, x0 + x2)
            sage: p2.projective(1)
            Point (x0, x1 + x2)
            sage: p2.projective(2)
            Point (x0, x1 + x2)
        """
        C = self.curve()
        A = C.ambient_space()
        ideal = self.prime_ideal()

        phi = A.projective_embedding(i)

        gs = list(phi.codomain().gens())
        xi = gs.pop(i)

        # gens of ideal is a groebner basis in degrevlex order
        S = phi.codomain().coordinate_ring()
        prime = S.ideal(ideal.subs(dict(zip(A.gens(), gs))).homogenize(xi))

        Cp = C.projective_closure(i)
        return Cp._closed_point(Cp, prime, self.degree())


class IntegralProjectiveCurveClosedPoint(IntegralCurveClosedPoint):
    """
    Closed points of projective plane curves.
    """
    def rational_point(self):
        """
        Return the rational point if this closed point is of degree `1`.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: F.<a> = GF(4)
            sage: P.<x,y,z> = ProjectiveSpace(F, 2)
            sage: C = Curve(x^3*y + y^3*z + x*z^3)
            sage: C.closed_points()
            [Point (x, z),
             Point (x, y),
             Point (y, z),
             Point (x + a*z, y + (a + 1)*z),
             Point (x + (a + 1)*z, y + a*z)]
            sage: [p.rational_point() for p in _]
            [(0 : 1 : 0), (0 : 0 : 1), (1 : 0 : 0), (a : a + 1 : 1), (a + 1 : a : 1)]
            sage: set(_) == set(C.rational_points())
            True
        """
        if self.degree() != 1:
            raise ValueError("not a rational point")

        C = self.curve()
        A = C.ambient_space().coordinate_ring()
        prime_ideal = self.prime_ideal()
        for i in range(A.ngens()):
            G = (prime_ideal + A.ideal([A.gen(i) - 1])).groebner_basis()
            if 1 not in G:
                break
        return C([A.gen(j).reduce(G) for j in range(A.ngens())])

    def affine(self, i=None):
        """
        Return the point in the ``i``-th affine patch of the curve.

        INPUT:

        - ``i`` -- integer; if not specified, it is chosen automatically

        EXAMPLES::

            sage: F.<a> = GF(2)
            sage: P.<x,y,z> = ProjectiveSpace(F, 2)
            sage: C = Curve(x^3*y + y^3*z + x*z^3)
            sage: p1, p2, p3 = C.closed_points()
            sage: p1.affine()
            Point (x, z)
            sage: p2.affine()
            Point (x, y)
            sage: p3.affine()
            Point (y, z)
            sage: p3.affine(0)
            Point (y, z)
            sage: p3.affine(1)
            Traceback (most recent call last):
            ...
            ValueError: not in the affine patch
        """
        C = self.curve()
        P = C.ambient_space()
        ideal = self.prime_ideal()
        if i is None:
            for j in range(P.ngens()):
                if P.gen(j) not in ideal:
                    i = j
                    break
        else:
            if P.gen(i) in ideal:
                raise ValueError("not in the affine patch")

        A = P.affine_patch(i)
        phi = A.projective_embedding(i, P)

        prime = A.coordinate_ring().ideal(ideal.subs(dict(zip(P.gens(), phi.defining_polynomials()))))

        Ca = C.affine_patch(i)
        return Ca._closed_point(Ca, prime, self.degree())
