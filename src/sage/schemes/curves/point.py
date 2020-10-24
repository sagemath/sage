"""
Rational points of curves

We can create points on projective curves::

    sage: P.<x,y,z,w> = ProjectiveSpace(QQ, 3)
    sage: C = Curve([x^3 - 2*x*z^2 - y^3, z^3 - w^3 - x*y*z], P)
    sage: Q = C([1,1,0,0])
    sage: Q.parent()
    Set of rational points of Projective Curve over Rational Field defined
    by x^3 - y^3 - 2*x*z^2, -x*y*z + z^3 - w^3

or on affine curves::

    sage: A.<x,y> = AffineSpace(GF(23), 2)
    sage: C = Curve([y - y^4 + 17*x^2 - 2*x + 22], A)
    sage: Q = C([22,21])
    sage: Q.parent()
    Set of rational points of Affine Plane Curve over Finite Field of size
    23 defined by -y^4 - 6*x^2 - 2*x + y - 1

AUTHORS:

- Grayson Jorgenson (2016-6): initial version

"""
#*****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from sage.schemes.affine.affine_point import (SchemeMorphism_point_affine_field,
                                              SchemeMorphism_point_affine_finite_field)
from sage.schemes.projective.projective_point import (SchemeMorphism_point_projective_field,
                                                      SchemeMorphism_point_projective_finite_field)

class ProjectiveCurvePoint_field(SchemeMorphism_point_projective_field):
    """
    Point of a projective curve over a field.
    """
    def is_singular(self):
        r"""
        Return whether this point is a singular point of the projective curve it is on.

        EXAMPLES::

            sage: P.<x,y,z,w> = ProjectiveSpace(QQ, 3)
            sage: C = Curve([x^2 - y^2, z - w], P)
            sage: Q1 = C([0,0,1,1])
            sage: Q1.is_singular()
            True
            sage: Q2 = C([1,1,1,1])
            sage: Q2.is_singular()
            False
        """
        return self.codomain().is_singular(self)


class ProjectivePlaneCurvePoint_field(ProjectiveCurvePoint_field):
    """
    Point of a projective plane curve over a field.
    """
    def multiplicity(self):
        r"""
        Return the multiplicity of this point with respect to the projective
        curve it is on.

        EXAMPLES::

            sage: P.<x,y,z> = ProjectiveSpace(GF(17), 2)
            sage: C = Curve([y^3*z - 16*x^4], P)
            sage: Q = C([0,0,1])
            sage: Q.multiplicity()
            3
        """
        return self.codomain().multiplicity(self)

    def tangents(self):
        r"""
        Return the tangents at this point of the projective plane curve this point is on.

        OUTPUT:

        A list of polynomials in the coordinate ring of the ambient space
        of the curve this point is on.

        EXAMPLES::

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: C = Curve([y^2*z^3 - x^5 + 18*y*x*z^3])
            sage: Q = C([0,0,1])
            sage: Q.tangents()
            [y, 18*x + y]
        """
        return self.codomain().tangents(self)

    def is_ordinary_singularity(self):
        r"""
        Return whether this point is an ordinary singularity of the projective
        plane curve it is on.

        EXAMPLES::

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: C = Curve([z^6 - x^6 - x^3*z^3 - x^3*y^3])
            sage: Q = C([0,1,0])
            sage: Q.is_ordinary_singularity()
            False

        ::

            sage: R.<a> = QQ[]
            sage: K.<b> = NumberField(a^2 - 3)
            sage: P.<x,y,z> = ProjectiveSpace(K, 2)
            sage: C = P.curve([x^2*y^3*z^4 - y^6*z^3 - 4*x^2*y^4*z^3 -
            ....: 4*x^4*y^2*z^3 + 3*y^7*z^2 + 10*x^2*y^5*z^2 + 9*x^4*y^3*z^2 +
            ....: 5*x^6*y*z^2 - 3*y^8*z - 9*x^2*y^6*z - 11*x^4*y^4*z - 7*x^6*y^2*z -
            ....: 2*x^8*z + y^9 + 2*x^2*y^7 + 3*x^4*y^5 + 4*x^6*y^3 + 2*x^8*y])
            sage: Q = C([-1/2, 1/2, 1])
            sage: Q.is_ordinary_singularity()
            True
        """
        return self.codomain().is_ordinary_singularity(self)

    def is_transverse(self, D):
        r"""
        Return whether the intersection of the curve ``D`` at this point with
        the curve this point is on is transverse or not.

        INPUT:

        - ``D`` -- a curve in the same ambient space as the curve this point is on

        EXAMPLES::

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: C = Curve([x^2 - 2*y^2 - 2*z^2], P)
            sage: D = Curve([y - z], P)
            sage: Q = C([2,1,1])
            sage: Q.is_transverse(D)
            True

        ::

            sage: P.<x,y,z> = ProjectiveSpace(GF(17), 2)
            sage: C = Curve([x^4 - 16*y^3*z], P)
            sage: D = Curve([y^2 - z*x], P)
            sage: Q = C([0,0,1])
            sage: Q.is_transverse(D)
            False
        """
        return self.codomain().is_transverse(D, self)


class ProjectivePlaneCurvePoint_finite_field(ProjectivePlaneCurvePoint_field,
                                             SchemeMorphism_point_projective_finite_field):
    """
    Point of a projective plane curve over a finite field.
    """
    pass


class IntegralProjectiveCurvePoint(ProjectiveCurvePoint_field):
    def closed_point(self):
        """
        Return the closed point corresponding to this rational point.

        EXAMPLES::

            sage: P.<x,y,z> = ProjectiveSpace(GF(17), 2)
            sage: C = Curve([x^4 - 16*y^3*z], P)
            sage: C.singular_points()
            [(0 : 0 : 1)]
            sage: p = _[0]
            sage: p.closed_point()
            Point (x, y)
        """
        curve = self.codomain()
        A = curve.ambient_space()
        S = A.coordinate_ring()

        hcoords = self._coords
        for i in range(S.ngens()):
            if hcoords[i]:
                break
        ai = hcoords[i]
        xi = S.gen(i)
        hgens = [ai*S.gen(j) - hcoords[j]*xi for j in range(S.ngens()) if j != i]
        return curve._closed_point(curve, S.ideal(hgens), degree=1)

    def places(self):
        """
        Return all places on this point.

        EXAMPLES::

            sage: P.<x,y,z> = ProjectiveSpace(GF(17), 2)
            sage: C = Curve([x^4 - 16*y^3*z], P)
            sage: C.singular_points()
            [(0 : 0 : 1)]
            sage: p = _[0]
            sage: p.places()
            [Place (y)]
        """
        return self.closed_point().places()

    def place(self):
        """
        Return a place on this point.

        EXAMPLES::

            sage: P.<x,y,z> = ProjectiveSpace(GF(17), 2)
            sage: C = Curve([x^4 - 16*y^3*z], P)
            sage: C.singular_points()
            [(0 : 0 : 1)]
            sage: p = _[0]
            sage: p.place()
            Place (y)
        """
        return self.closed_point().place()


class IntegralProjectiveCurvePoint_finite_field(IntegralProjectiveCurvePoint):
    """
    Point of an integral projective curve over a finite field.
    """
    pass


class IntegralProjectivePlaneCurvePoint(IntegralProjectiveCurvePoint, ProjectivePlaneCurvePoint_field):
    """
    Point of an integral projective plane curve over a field.
    """
    pass


class IntegralProjectivePlaneCurvePoint_finite_field(ProjectivePlaneCurvePoint_finite_field, IntegralProjectiveCurvePoint_finite_field):
    """
    Point of an integral projective plane curve over a finite field.
    """
    pass


class AffineCurvePoint_field(SchemeMorphism_point_affine_field):

    def is_singular(self):
        r"""
        Return whether this point is a singular point of the affine curve it is on.

        EXAMPLES::

            sage: K = QuadraticField(-1)
            sage: A.<x,y,z> = AffineSpace(K, 3)
            sage: C = Curve([(x^4 + 2*z + 2)*y, z - y + 1])
            sage: Q1 = C([0,0,-1])
            sage: Q1.is_singular()
            True
            sage: Q2 = C([-K.gen(),0,-1])
            sage: Q2.is_singular()
            False
        """
        return self.codomain().is_singular(self)


class AffinePlaneCurvePoint_field(AffineCurvePoint_field):
    """
    Point of an affine plane curve over a field.
    """
    def multiplicity(self):
        r"""
        Return the multiplicity of this point with respect to the affine curve it is on.

        EXAMPLES::

            sage: A.<x,y> = AffineSpace(QQ, 2)
            sage: C = A.curve([2*x^7 - 3*x^6*y + x^5*y^2 + 31*x^6 - 40*x^5*y +
            ....: 13*x^4*y^2 - x^3*y^3 + 207*x^5 - 228*x^4*y + 70*x^3*y^2 - 7*x^2*y^3
            ....: + 775*x^4 - 713*x^3*y + 193*x^2*y^2 - 19*x*y^3 + y^4 + 1764*x^3 -
            ....: 1293*x^2*y + 277*x*y^2 - 22*y^3 + 2451*x^2 - 1297*x*y + 172*y^2 +
            ....: 1935*x - 570*y + 675])
            sage: Q = C([-2,1])
            sage: Q.multiplicity()
            4
        """
        return self.codomain().multiplicity(self)

    def tangents(self):
        r"""
        Return the tangents at this point of the affine plane curve this point is on.

        OUTPUT: a list of polynomials in the coordinate ring of the ambient
        space of the curve this point is on.

        EXAMPLES::

            sage: A.<x,y> = AffineSpace(QQ, 2)
            sage: C = A.curve([x^5 - x^3*y^2 + 5*x^4 - x^3*y - 3*x^2*y^2 +
            ....: x*y^3 + 10*x^3 - 3*x^2*y - 3*x*y^2 + y^3 + 10*x^2 - 3*x*y - y^2 +
            ....: 5*x - y + 1])
            sage: Q = C([-1,0])
            sage: Q.tangents()
            [y, x + 1, x - y + 1, x + y + 1]
        """
        return self.codomain().tangents(self)

    def is_ordinary_singularity(self):
        r"""
        Return whether this point is an ordinary singularity of the affine
        plane curve it is on.

        EXAMPLES::

            sage: A.<x,y> = AffineSpace(QQ, 2)
            sage: C = A.curve([x^5 - x^3*y^2 + 5*x^4 - x^3*y - 3*x^2*y^2 +
            ....: x*y^3 + 10*x^3 - 3*x^2*y - 3*x*y^2 + y^3 + 10*x^2 - 3*x*y - y^2 +
            ....: 5*x - y + 1])
            sage: Q = C([-1,0])
            sage: Q.is_ordinary_singularity()
            True

        ::

            sage: A.<x,y> = AffineSpace(GF(7), 2)
            sage: C = A.curve([y^2 - x^7 - 6*x^3])
            sage: Q = C([0,0])
            sage: Q.is_ordinary_singularity()
            False
        """
        return self.codomain().is_ordinary_singularity(self)

    def is_transverse(self, D):
        r"""
        Return whether the intersection of the curve ``D`` at this point with
        the curve this point is on is transverse or not.

        INPUT:

        - ``D`` -- a curve in the same ambient space as the curve this point is on.

        EXAMPLES::

            sage: A.<x,y> = AffineSpace(QQ, 2)
            sage: C = Curve([y - x^2], A)
            sage: D = Curve([y], A)
            sage: Q = C([0,0])
            sage: Q.is_transverse(D)
            False

        ::

            sage: R.<a> = QQ[]
            sage: K.<b> = NumberField(a^2 - 2)
            sage: A.<x,y> = AffineSpace(K, 2)
            sage: C = Curve([y^2 + x^2 - 1], A)
            sage: D = Curve([y - x], A)
            sage: Q = C([-1/2*b,-1/2*b])
            sage: Q.is_transverse(D)
            True
        """
        return self.codomain().is_transverse(D, self)


class AffinePlaneCurvePoint_finite_field(AffinePlaneCurvePoint_field, SchemeMorphism_point_affine_finite_field):
    """
    Point of an affine plane curve over a finite field.
    """
    pass


class IntegralAffineCurvePoint(AffineCurvePoint_field):
    """
    Point of an integral affine curve.
    """
    def closed_point(self):
        """
        Return the closed point that corresponds to this rational point.

        EXAMPLES::

            sage: A.<x,y> = AffineSpace(GF(8), 2)
            sage: C = Curve(x^5 + y^5 + x*y + 1)
            sage: p = C([1,1])
            sage: p.closed_point()
            Point (x + 1, y + 1)
        """
        curve = self.codomain()
        A = curve.ambient_space()
        R = A.coordinate_ring()
        coords = self._coords
        gens = [R.gen(i) - coords[i] for i in range(R.ngens())]
        return curve._closed_point(curve, R.ideal(gens), degree=1)

    def places(self):
        """
        Return all places on this point.

        EXAMPLES::

            sage: A.<x,y> = AffineSpace(GF(2), 2)
            sage: C = Curve(x^5 + y^5 + x*y + 1)
            sage: p = C(-1,-1)
            sage: p
            (1, 1)
            sage: p.closed_point()
            Point (x + 1, y + 1)
            sage: _.places()
            [Place (x + 1, (1/(x^5 + 1))*y^4 + ((x^5 + x^4 + 1)/(x^5 + 1))*y^3
            + ((x^5 + x^3 + 1)/(x^5 + 1))*y^2 + (x^2/(x^5 + 1))*y), Place (x +
            1, (1/(x^5 + 1))*y^4 + ((x^5 + x^4 + 1)/(x^5 + 1))*y^3 + (x^3/(x^5
            + 1))*y^2 + (x^2/(x^5 + 1))*y + x + 1)]
        """
        return self.closed_point().places()

    def place(self):
        """
        Return a place on this point.

        EXAMPLES::

            sage: A.<x,y> = AffineSpace(GF(2), 2)
            sage: C = Curve(x^5 + y^5 + x*y + 1)
            sage: p = C(-1,-1)
            sage: p
            (1, 1)
            sage: p.closed_point()
            Point (x + 1, y + 1)
            sage: _.place()
            Place (x + 1, (1/(x^5 + 1))*y^4 + ((x^5 + x^4 + 1)/(x^5 + 1))*y^3 +
            ((x^5 + x^3 + 1)/(x^5 + 1))*y^2 + (x^2/(x^5 + 1))*y)
        """
        return self.closed_point().place()


class IntegralAffineCurvePoint_finite_field(IntegralAffineCurvePoint):
    """
    Point of an integral affine curve over a finite field.
    """
    pass


class IntegralAffinePlaneCurvePoint(IntegralAffineCurvePoint, AffinePlaneCurvePoint_field):
    """
    Point of an integral affine plane curve over a finite field.
    """
    pass

class IntegralAffinePlaneCurvePoint_finite_field(AffinePlaneCurvePoint_finite_field, IntegralAffineCurvePoint_finite_field):
    """
    Point of an integral affine plane curve over a finite field.
    """
    pass
