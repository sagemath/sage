"""
Affine curves

EXAMPLES:

We can construct curves in either an affine plane::

    sage: A.<x,y> = AffineSpace(QQ, 2)
    sage: C = Curve([y - x^2], A); C
    Affine Plane Curve over Rational Field defined by -x^2 + y

or in higher dimensional affine space::

    sage: A.<x,y,z,w> = AffineSpace(QQ, 4)
    sage: C = Curve([y - x^2, z - w^3, w - y^4], A); C
    Affine Curve over Rational Field defined by -x^2 + y, -w^3 + z, -y^4 + w

AUTHORS:

- William Stein (2005-11-13)

- David Joyner (2005-11-13)

- David Kohel (2006-01)

- Grayson Jorgenson (2016-8)

"""
# ****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#  The full text of the GPL is available at:
#
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from __future__ import absolute_import

from sage.arith.misc import binomial
from sage.categories.fields import Fields
from sage.categories.finite_fields import FiniteFields
from sage.categories.homset import Hom, End
from sage.categories.number_fields import NumberFields
from sage.interfaces.all import singular

from sage.misc.all import add

from sage.rings.all import degree_lowest_rational_function

from sage.rings.number_field.number_field import NumberField
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.qqbar import number_field_elements_from_algebraics, QQbar
from sage.rings.rational_field import is_RationalField
from sage.schemes.affine.affine_space import (AffineSpace,
                                              is_AffineSpace)
from . import point

from sage.schemes.affine.affine_subscheme import AlgebraicScheme_subscheme_affine
from .curve import Curve_generic


class AffineCurve(Curve_generic, AlgebraicScheme_subscheme_affine):

    _point = point.AffineCurvePoint_field

    def _repr_type(self):
        r"""
        Return a string representation of the type of this curve.

        EXAMPLES::

            sage: A.<x,y,z,w> = AffineSpace(QQ, 4)
            sage: C = Curve([x - y, z - w, w - x], A)
            sage: C._repr_type()
            'Affine'
        """
        return "Affine"

    def __init__(self, A, X):
        r"""
        Initialization function.

        EXAMPLES::

            sage: R.<v> = QQ[]
            sage: K.<u> = NumberField(v^2 + 3)
            sage: A.<x,y,z> = AffineSpace(K, 3)
            sage: C = Curve([z - u*x^2, y^2], A); C
            Affine Curve over Number Field in u with defining polynomial v^2 + 3
            defined by (-u)*x^2 + z, y^2

        ::

            sage: A.<x,y,z> = AffineSpace(GF(7), 3)
            sage: C = Curve([x^2 - z, z - 8*x], A); C
            Affine Curve over Finite Field of size 7 defined by x^2 - z, -x + z
        """
        if not is_AffineSpace(A):
            raise TypeError("A (=%s) must be an affine space" % A)
        Curve_generic.__init__(self, A, X)
        d = self.dimension()
        if d != 1:
            raise ValueError("defining equations (=%s) define a scheme of dimension %s != 1" % (X, d))

    def projective_closure(self, i=0, PP=None):
        r"""
        Return the projective closure of this affine curve.

        INPUT:

        - ``i`` -- (default: 0) the index of the affine coordinate chart of the projective space that the affine
          ambient space of this curve embeds into.

        - ``PP`` -- (default: None) ambient projective space to compute the projective closure in. This is
          constructed if it is not given.

        OUTPUT:

        - a curve in projective space.

        EXAMPLES::

            sage: A.<x,y,z> = AffineSpace(QQ, 3)
            sage: C = Curve([y-x^2,z-x^3], A)
            sage: C.projective_closure()
            Projective Curve over Rational Field defined by x1^2 - x0*x2,
            x1*x2 - x0*x3, x2^2 - x1*x3

        ::

            sage: A.<x,y,z> = AffineSpace(QQ, 3)
            sage: C = Curve([y - x^2, z - x^3], A)
            sage: C.projective_closure()
            Projective Curve over Rational Field defined by
            x1^2 - x0*x2, x1*x2 - x0*x3, x2^2 - x1*x3

        ::

            sage: A.<x,y> = AffineSpace(CC, 2)
            sage: C = Curve(y - x^3 + x - 1, A)
            sage: C.projective_closure(1)
            Projective Plane Curve over Complex Field with 53 bits of precision defined by
            x0^3 - x0*x1^2 + x1^3 - x1^2*x2

        ::

            sage: A.<x,y> = AffineSpace(QQ, 2)
            sage: P.<u,v,w> = ProjectiveSpace(QQ, 2)
            sage: C = Curve([y - x^2], A)
            sage: C.projective_closure(1, P).ambient_space() == P
            True
        """
        from .constructor import Curve
        return Curve(AlgebraicScheme_subscheme_affine.projective_closure(self, i, PP))

    def projection(self, indices, AS=None):
        r"""
        Return the projection of this curve onto the coordinates specified by ``indices``.

        INPUT:

        - ``indices`` -- a list or tuple of distinct integers specifying the indices of the coordinates to use
          in the projection. Can also be a list or tuple consisting of variables of the coordinate ring of the
          ambient space of this curve. If integers are used to specify the coordinates, 0 denotes the first
          coordinate. The length of ``indices`` must be between two and one less than the dimension of the ambient
          space of this curve, inclusive.

        - ``AS`` -- (default: None) the affine space the projected curve will be defined in. This space must be
          defined over the same base field as this curve, and must have dimension equal to the length of ``indices``.
          This space is constructed if not specified.

        OUTPUT:

        - a tuple consisting of two elements: a scheme morphism from this curve to affine space of dimension
          equal to the number of coordinates specified in ``indices``, and the affine subscheme that is the image
          of that morphism. If the image is a curve, the second element of the tuple will be a curve.

        EXAMPLES::

            sage: A.<x,y,z> = AffineSpace(QQ, 3)
            sage: C = Curve([y^7 - x^2 + x^3 - 2*z, z^2 - x^7 - y^2], A)
            sage: C.projection([0,1])
            (Scheme morphism:
               From: Affine Curve over Rational Field defined by y^7 + x^3 - x^2 -
            2*z, -x^7 - y^2 + z^2
               To:   Affine Space of dimension 2 over Rational Field
               Defn: Defined on coordinates by sending (x, y, z) to
                     (x, y),
             Affine Plane Curve over Rational Field defined by x1^14 + 2*x0^3*x1^7 -
            2*x0^2*x1^7 - 4*x0^7 + x0^6 - 2*x0^5 + x0^4 - 4*x1^2)
            sage: C.projection([0,1,3,4])
            Traceback (most recent call last):
            ...
            ValueError: (=[0, 1, 3, 4]) must be a list or tuple of length between 2
            and (=2), inclusive

        ::

            sage: A.<x,y,z,w> = AffineSpace(QQ, 4)
            sage: C = Curve([x - 2, y - 3, z - 1], A)
            sage: B.<a,b,c> = AffineSpace(QQ, 3)
            sage: C.projection([0,1,2], AS=B)
            (Scheme morphism:
               From: Affine Curve over Rational Field defined by x - 2, y - 3, z - 1
               To:   Affine Space of dimension 3 over Rational Field
               Defn: Defined on coordinates by sending (x, y, z, w) to
                     (x, y, z),
             Closed subscheme of Affine Space of dimension 3 over Rational Field
            defined by:
               c - 1,
               b - 3,
               a - 2)

        ::

            sage: A.<x,y,z,w,u> = AffineSpace(GF(11), 5)
            sage: C = Curve([x^3 - 5*y*z + u^2, x - y^2 + 3*z^2, w^2 + 2*u^3*y, y - u^2 + z*x], A)
            sage: B.<a,b,c> = AffineSpace(GF(11), 3)
            sage: proj1 = C.projection([1,2,4], AS=B)
            sage: proj1
            (Scheme morphism:
               From: Affine Curve over Finite Field of size 11 defined by x^3 -
            5*y*z + u^2, -y^2 + 3*z^2 + x, 2*y*u^3 + w^2, x*z - u^2 + y
               To:   Affine Space of dimension 3 over Finite Field of size 11
               Defn: Defined on coordinates by sending (x, y, z, w, u) to
                     (y, z, u),
             Affine Curve over Finite Field of size 11 defined by a^2*b - 3*b^3 -
            c^2 + a, c^6 - 5*a*b^4 + b^3*c^2 - 3*a*c^4 + 3*a^2*c^2 - a^3, a^2*c^4 -
            3*b^2*c^4 - 2*a^3*c^2 - 5*a*b^2*c^2 + a^4 - 5*a*b^3 + 2*b^4 + b^2*c^2 -
            3*b*c^2 + 3*a*b, a^4*c^2 + 2*b^4*c^2 - a^5 - 2*a*b^4 + 5*b*c^4 + a*b*c^2
            - 5*a*b^2 + 4*b^3 + b*c^2 + 5*c^2 - 5*a, a^6 - 5*b^6 - 5*b^3*c^2 +
            5*a*b^3 + 2*c^4 - 4*a*c^2 + 2*a^2 - 5*a*b + c^2)
            sage: proj1[1].ambient_space() is B
            True
            sage: proj2 = C.projection([1,2,4])
            sage: proj2[1].ambient_space() is B
            False
            sage: C.projection([1,2,3,5], AS=B)
            Traceback (most recent call last):
            ...
            TypeError: (=Affine Space of dimension 3 over Finite Field of size 11)
            must have dimension (=4)

        ::

            sage: A.<x,y,z,w> = AffineSpace(QQ, 4)
            sage: C = A.curve([x*y - z^3, x*z - w^3, w^2 - x^3])
            sage: C.projection([y,z])
            (Scheme morphism:
               From: Affine Curve over Rational Field defined by -z^3 + x*y, -w^3 +
            x*z, -x^3 + w^2
               To:   Affine Space of dimension 2 over Rational Field
               Defn: Defined on coordinates by sending (x, y, z, w) to
                     (y, z),
             Affine Plane Curve over Rational Field defined by x1^23 - x0^7*x1^4)
            sage: B.<x,y,z> = AffineSpace(QQ, 3)
            sage: C.projection([x,y,z], AS=B)
            (Scheme morphism:
               From: Affine Curve over Rational Field defined by -z^3 + x*y, -w^3 +
            x*z, -x^3 + w^2
               To:   Affine Space of dimension 3 over Rational Field
               Defn: Defined on coordinates by sending (x, y, z, w) to
                     (x, y, z),
             Affine Curve over Rational Field defined by z^3 - x*y, x^8 - x*z^2,
            x^7*z^2 - x*y*z)
            sage: C.projection([y,z,z])
            Traceback (most recent call last):
            ...
            ValueError: (=[y, z, z]) must be a list or tuple of distinct indices or
            variables
        """
        AA = self.ambient_space()
        n = AA.dimension_relative()
        if n == 2:
            raise TypeError("this curve is already a plane curve")
        if self.base_ring() not in Fields():
            raise TypeError("this curve must be defined over a field")
        if len(indices) < 2 or len(indices) >= n:
            raise ValueError("(=%s) must be a list or tuple of length between 2 and (=%s), inclusive" % (indices, n - 1))
        if len(set(indices)) < len(indices):
            raise ValueError("(=%s) must be a list or tuple of distinct indices or variables" % indices)
        if not AS is None:
            if not is_AffineSpace(AS):
                raise TypeError("(=%s) must be an affine space" % AS)
            if AS.dimension_relative() != len(indices):
                raise TypeError("(=%s) must have dimension (=%s)" % (AS, len(indices)))
            if AS.base_ring() != AA.base_ring():
                raise TypeError("(=%s) must be defined over the same base field as this curve" % AS)
        indices = list(indices)
        if all(f in AA.gens() for f in indices):
            indices = [AA.gens().index(f) for f in indices]
            indices.sort()
        else:
            indices = [int(i) for i in indices] # type checking
            indices.sort()
            if indices[0] < 0 or indices[-1] > n - 1:
                raise ValueError("index values must be between 0 and one minus the dimension of the ambient space " \
                                 "of this curve")
        # construct the projection map
        if AS is None:
            AA2 = AffineSpace(self.base_ring(), len(indices))
        else:
            AA2 = AS
        H = Hom(self, AA2)
        psi = H([AA.gens()[i] for i in indices])
        # compute the image via elimination
        removecoords = list(AA.gens())
        for i in range(len(indices) - 1, -1, -1):
            removecoords.pop(indices[i])
        J = self.defining_ideal().elimination_ideal(removecoords)
        K = Hom(AA.coordinate_ring(), AA2.coordinate_ring())
        l = [0]*(n)
        for i in range(len(indices)):
            l[indices[i]] = AA2.gens()[i]
        phi = K(l)
        G = [phi(f) for f in J.gens()]
        try:
            C = AA2.curve(G)
        except (TypeError, ValueError):
            C = AA2.subscheme(G)
        return tuple([psi, C])

    def plane_projection(self, AP=None):
        r"""
        Return a projection of this curve into an affine plane so that the image of the projection is
        a plane curve.

        INPUT:

        - ``AP`` -- (default: None) the affine plane to project this curve into. This space must be defined over
          the same base field as this curve, and must have dimension two. This space will be constructed if not
          specified.

        OUTPUT:

        - a tuple consisting of two elements: a scheme morphism from this curve into an affine plane, and the plane
          curve that defines the image of that morphism.

        EXAMPLES::

            sage: A.<x,y,z,w> = AffineSpace(QQ, 4)
            sage: C = Curve([x^2 - y*z*w, z^3 - w, w + x*y - 3*z^3], A)
            sage: C.plane_projection()
            (Scheme morphism:
              From: Affine Curve over Rational Field defined by -y*z*w + x^2, z^3 -
            w, -3*z^3 + x*y + w
              To:   Affine Space of dimension 2 over Rational Field
              Defn: Defined on coordinates by sending (x, y, z, w) to
                    (x, y), Affine Plane Curve over Rational Field defined by
            x0^2*x1^7 - 16*x0^4)

        ::

            sage: R.<a> = QQ[]
            sage: K.<b> = NumberField(a^2 + 2)
            sage: A.<x,y,z> = AffineSpace(K, 3)
            sage: C = A.curve([x - b, y - 2])
            sage: B.<a,b> = AffineSpace(K, 2)
            sage: proj1 = C.plane_projection(AP=B)
            sage: proj1
            (Scheme morphism:
               From: Affine Curve over Number Field in b with defining polynomial
            a^2 + 2 defined by x + (-b), y - 2
               To:   Affine Space of dimension 2 over Number Field in b with
            defining polynomial a^2 + 2
               Defn: Defined on coordinates by sending (x, y, z) to
                     (x, z),
             Affine Plane Curve over Number Field in b with defining polynomial a^2
            + 2 defined by a + (-b))
            sage: proj1[1].ambient_space() is B
            True
            sage: proj2 = C.plane_projection()
            sage: proj2[1].ambient_space() is B
            False
        """
        n = self.ambient_space().dimension_relative()
        # finds a projection that will have a plane curve as its image
        # the following iterates over all pairs (i,j) with j > i to test all
        # possible projections
        for i in range(0, n - 1):
            for j in range(i + 1, n):
                L = self.projection([i,j], AP)
                if isinstance(L[1], Curve_generic):
                    return L

    def blowup(self, P=None):
        r"""
        Return the blow up of this affine curve at the point ``P``.

        The blow up is described by affine charts. This curve must be irreducible.

        INPUT:

        - ``P`` -- (default: None) a point on this curve at which to blow up. If ``None``, then ``P`` is
          taken to be the origin.

        OUTPUT:

        - a tuple consisting of three elements. The first is a tuple of curves in affine space of
          the same dimension as the ambient space of this curve, which define the blow up in
          each affine chart. The second is a tuple of tuples such that the jth element of the ith
          tuple is the transition map from the ith affine patch to the jth affine patch. Lastly,
          the third element is a tuple consisting of the restrictions of the projection map from
          the blow up back to the original curve, restricted to each affine patch. There the
          ith element will be the projection from the ith affine patch.

        EXAMPLES::

            sage: A.<x,y> = AffineSpace(QQ, 2)
            sage: C = Curve([y^2 - x^3], A)
            sage: C.blowup()
            ((Affine Plane Curve over Rational Field defined by s1^2 - x,
             Affine Plane Curve over Rational Field defined by y*s0^3 - 1),
            ([Scheme endomorphism of Affine Plane Curve over Rational Field defined by s1^2 - x
                Defn: Defined on coordinates by sending (x, s1) to
                      (x, s1), Scheme morphism:
                From: Affine Plane Curve over Rational Field defined by s1^2 - x
                To:   Affine Plane Curve over Rational Field defined by y*s0^3 - 1
                Defn: Defined on coordinates by sending (x, s1) to
                      (x*s1, 1/s1)], [Scheme morphism:
                From: Affine Plane Curve over Rational Field defined by y*s0^3 - 1
                To:   Affine Plane Curve over Rational Field defined by s1^2 - x
                Defn: Defined on coordinates by sending (y, s0) to
                      (y*s0, 1/s0),
              Scheme endomorphism of Affine Plane Curve over Rational Field defined by y*s0^3 - 1
                Defn: Defined on coordinates by sending (y, s0) to
                      (y, s0)]),
            (Scheme morphism:
               From: Affine Plane Curve over Rational Field defined by s1^2 - x
               To:   Affine Plane Curve over Rational Field defined by -x^3 + y^2
               Defn: Defined on coordinates by sending (x, s1) to
                     (x, x*s1), Scheme morphism:
               From: Affine Plane Curve over Rational Field defined by y*s0^3 - 1
               To:   Affine Plane Curve over Rational Field defined by -x^3 + y^2
               Defn: Defined on coordinates by sending (y, s0) to
                     (y*s0, y)))

        ::

            sage: K.<a> = QuadraticField(2)
            sage: A.<x,y,z> = AffineSpace(K, 3)
            sage: C = Curve([y^2 - a*x^5, x - z], A)
            sage: B = C.blowup()
            sage: B[0]
            (Affine Curve over Number Field in a with defining polynomial x^2 - 2 with a = 1.414213562373095? defined by s2 - 1, 2*x^3 + (-a)*s1^2,
             Affine Curve over Number Field in a with defining polynomial x^2 - 2 with a = 1.414213562373095? defined by s0 - s2, 2*y^3*s2^5 + (-a),
             Affine Curve over Number Field in a with defining polynomial x^2 - 2 with a = 1.414213562373095? defined by s0 - 1, 2*z^3 + (-a)*s1^2)
            sage: B[1][0][2]
            Scheme morphism:
              From: Affine Curve over Number Field in a with defining polynomial x^2 - 2 with a = 1.414213562373095? defined by s2 - 1, 2*x^3 + (-a)*s1^2
              To:   Affine Curve over Number Field in a with defining polynomial x^2 - 2 with a = 1.414213562373095? defined by s0 - 1, 2*z^3 + (-a)*s1^2
              Defn: Defined on coordinates by sending (x, s1, s2) to
                    (x*s2, 1/s2, s1/s2)
            sage: B[1][2][0]
            Scheme morphism:
              From: Affine Curve over Number Field in a with defining polynomial x^2 - 2 with a = 1.414213562373095? defined by s0 - 1, 2*z^3 + (-a)*s1^2
              To:   Affine Curve over Number Field in a with defining polynomial x^2 - 2 with a = 1.414213562373095? defined by s2 - 1, 2*x^3 + (-a)*s1^2
              Defn: Defined on coordinates by sending (z, s0, s1) to
                    (z*s0, s1/s0, 1/s0)
            sage: B[2]
            (Scheme morphism:
               From: Affine Curve over Number Field in a with defining polynomial x^2 - 2 with a = 1.414213562373095? defined by s2 - 1, 2*x^3 + (-a)*s1^2
               To:   Affine Curve over Number Field in a with defining polynomial x^2 - 2 with a = 1.414213562373095? defined by (-a)*x^5 + y^2, x - z
               Defn: Defined on coordinates by sending (x, s1, s2) to
                     (x, x*s1, x*s2), Scheme morphism:
               From: Affine Curve over Number Field in a with defining polynomial x^2 - 2 with a = 1.414213562373095? defined by s0 - s2, 2*y^3*s2^5 + (-a)
               To:   Affine Curve over Number Field in a with defining polynomial x^2 - 2 with a = 1.414213562373095? defined by (-a)*x^5 + y^2, x - z
               Defn: Defined on coordinates by sending (y, s0, s2) to
                     (y*s0, y, y*s2), Scheme morphism:
               From: Affine Curve over Number Field in a with defining polynomial x^2 - 2 with a = 1.414213562373095? defined by s0 - 1, 2*z^3 + (-a)*s1^2
               To:   Affine Curve over Number Field in a with defining polynomial x^2 - 2 with a = 1.414213562373095? defined by (-a)*x^5 + y^2, x - z
               Defn: Defined on coordinates by sending (z, s0, s1) to
                     (z*s0, z*s1, z))

        ::

            sage: A.<x,y> = AffineSpace(QQ, 2)
            sage: C = A.curve((y - 3/2)^3 - (x + 2)^5 - (x + 2)^6)
            sage: Q = A([-2,3/2])
            sage: C.blowup(Q)
            ((Affine Plane Curve over Rational Field defined by x^3 - s1^3 + 7*x^2 + 16*x + 12,
              Affine Plane Curve over Rational Field defined by 8*y^3*s0^6 - 36*y^2*s0^6 + 8*y^2*s0^5 +
              54*y*s0^6 - 24*y*s0^5 - 27*s0^6 + 18*s0^5 - 8),
             ([Scheme endomorphism of Affine Plane Curve over Rational Field defined by x^3 - s1^3 + 7*x^2 +
             16*x + 12
                 Defn: Defined on coordinates by sending (x, s1) to
                       (x, s1), Scheme morphism:
                 From: Affine Plane Curve over Rational Field defined by x^3 - s1^3 + 7*x^2 + 16*x + 12
                 To:   Affine Plane Curve over Rational Field defined by 8*y^3*s0^6 - 36*y^2*s0^6 + 8*y^2*s0^5 +
                 54*y*s0^6 - 24*y*s0^5 - 27*s0^6 + 18*s0^5 - 8
                 Defn: Defined on coordinates by sending (x, s1) to
                       (x*s1 + 2*s1 + 3/2, 1/s1)], [Scheme morphism:
                 From: Affine Plane Curve over Rational Field defined by 8*y^3*s0^6 - 36*y^2*s0^6 + 8*y^2*s0^5 +
                 54*y*s0^6 - 24*y*s0^5 - 27*s0^6 + 18*s0^5 - 8
                 To:   Affine Plane Curve over Rational Field defined by x^3 - s1^3 + 7*x^2 + 16*x + 12
                 Defn: Defined on coordinates by sending (y, s0) to
                       (y*s0 - 3/2*s0 - 2, 1/s0),
               Scheme endomorphism of Affine Plane Curve over Rational Field defined by 8*y^3*s0^6 - 36*y^2*s0^6 +
               8*y^2*s0^5 + 54*y*s0^6 - 24*y*s0^5 - 27*s0^6 + 18*s0^5 - 8
                 Defn: Defined on coordinates by sending (y, s0) to
                       (y, s0)]),
             (Scheme morphism:
                From: Affine Plane Curve over Rational Field defined by x^3 - s1^3 + 7*x^2 + 16*x + 12
                To:   Affine Plane Curve over Rational Field defined by -x^6 - 13*x^5 - 70*x^4 - 200*x^3 + y^3 -
                320*x^2 - 9/2*y^2 - 272*x + 27/4*y - 795/8
                Defn: Defined on coordinates by sending (x, s1) to
                      (x, x*s1 + 2*s1 + 3/2), Scheme morphism:
                From: Affine Plane Curve over Rational Field defined by 8*y^3*s0^6 - 36*y^2*s0^6 + 8*y^2*s0^5 +
                54*y*s0^6 - 24*y*s0^5 - 27*s0^6 + 18*s0^5 - 8
                To:   Affine Plane Curve over Rational Field defined by -x^6 - 13*x^5 - 70*x^4 - 200*x^3 + y^3 -
                320*x^2 - 9/2*y^2 - 272*x + 27/4*y - 795/8
                Defn: Defined on coordinates by sending (y, s0) to
                      (y*s0 - 3/2*s0 - 2, y)))

        ::

            sage: A.<x,y,z,w> = AffineSpace(QQ, 4)
            sage: C = A.curve([((x + 1)^2 + y^2)^3 - 4*(x + 1)^2*y^2, y - z, w - 4])
            sage: Q = C([-1,0,0,4])
            sage: B = C.blowup(Q)
            sage: B[0]
            (Affine Curve over Rational Field defined by s3, s1 - s2, x^2*s2^6 +
            2*x*s2^6 + 3*x^2*s2^4 + s2^6 + 6*x*s2^4 + 3*x^2*s2^2 + 3*s2^4 + 6*x*s2^2
            + x^2 - s2^2 + 2*x + 1,
             Affine Curve over Rational Field defined by s3, s2 - 1, y^2*s0^6 +
            3*y^2*s0^4 + 3*y^2*s0^2 + y^2 - 4*s0^2,
             Affine Curve over Rational Field defined by s3, s1 - 1, z^2*s0^6 +
            3*z^2*s0^4 + 3*z^2*s0^2 + z^2 - 4*s0^2,
             Closed subscheme of Affine Space of dimension 4 over Rational Field
            defined by:
               1)
            sage: Q = A([6,2,3,1])
            sage: B = C.blowup(Q)
            Traceback (most recent call last):
            ...
            TypeError: (=(6, 2, 3, 1)) must be a point on this curve

        ::

            sage: A.<x,y> = AffineSpace(QuadraticField(-1), 2)
            sage: C = A.curve([y^2 + x^2])
            sage: C.blowup()
            Traceback (most recent call last):
            ...
            TypeError: this curve must be irreducible
        """
        A = self.ambient_space()
        n = A.dimension_relative()
        if P is None:
            P = A([0]*n)
        try:
            self(P)
        except TypeError:
            raise TypeError("(=%s) must be a point on this curve" % P)
        if not self.base_ring() in Fields():
            raise TypeError("the base ring of this curve must be a field")
        if not self.is_irreducible():
            raise TypeError("this curve must be irreducible")
        # attempt to make the variable names more organized
        # the convention used here is to have the homogeneous coordinates for the projective component of the
        # product space the blow up resides in be generated from the letter 's'. The following loop is in place
        # to prevent conflicts in the names from occurring
        rf = 1
        for i in range(n):
            if str(A.gens()[i])[0] == 's' and len(str(A.gens()[i])) > rf:
                rf = len(str(A.gens()[i]))
        var_names = [str(A.gens()[i]) for i in range(n)] + ['s'*rf + str(i) for i in range(n)]
        R = PolynomialRing(A.base_ring(), 2*n, var_names)
        # move the defining polynomials of this curve into R
        H = Hom(A.coordinate_ring(), R)
        psi = H([R.gens()[i] for i in range(n)])
        n_polys = [psi(f) for f in self.defining_polynomials()]
        # the blow up ideal of A at P is the ideal generated by
        # (z_i - p_i)*s_j - (z_j - p_j)*s_i for i != j from 0,...,n-1
        # in the mixed product space of A^n and P^{n-1} where the z_i are the gens
        # of A^n, the s_i are the gens for P^{n-1}, and P = (p_1,...,p_n). We describe the
        # blow up of this curve at P in each affine chart
        patches = []
        for i in range(n):
            # in this chart, s_i is assumed to be 1
            # substitute in z_j = (z_i - p_i)*s_j + p_j for each j != i
            coords = list(R.gens())
            for j in range(n):
                if j != i:
                    coords[j] = (R.gens()[i] - P[i])*R.gens()[j + n] + P[j]
            c_polys = [f(coords) for f in n_polys]
            var_names = list(R.gens())[n:2*n]
            var_names.pop(i)
            var_names.insert(0, R.gens()[i])
            c_A = AffineSpace(R.base_ring(), n, var_names)
            H = Hom(R, c_A.coordinate_ring())
            coords = [0]*(2*n)
            coords[i] = c_A.gens()[0]
            t = 1
            for j in range(n):
                if j != i:
                    coords[j + n] = c_A.gens()[t]
                    t = t + 1
                else:
                    coords[j + n] = 1
            psi = H(coords)
            c_polys = [psi(f) for f in c_polys]
            # choose the component of the subscheme defined by these polynomials
            # that corresponds to the proper transform
            irr_comps = c_A.subscheme(c_polys).irreducible_components()
            for j in range(len(irr_comps)):
                proper_transform = True
                for f in irr_comps[j].defining_polynomials():
                    if (c_A.gens()[0] - P[i]).divides(f):
                        proper_transform = False
                        break
                if proper_transform:
                    patches.append(c_A.curve(irr_comps[j].defining_polynomials()))
                    break
                elif j + 1 == len(irr_comps):
                    # patch of blowup in this chart is empty
                    patches.append(c_A.subscheme(1))
        # create the transition maps between the charts
        t_maps = []
        for i in range(n):
            maps = []
            for j in range(n):
                AA = patches[i].ambient_space()
                H = Hom(patches[i], patches[j])
                vars = AA.gens()
                homvars = list(AA.gens())
                homvars.pop(0)
                homvars.insert(i, 1)
                coords = [(vars[0] - P[i])*homvars[j] + P[j]]
                for t in range(n):
                    if t != j:
                        coords.append(homvars[t]/homvars[j])
                maps.append(H(coords))
            t_maps.append(maps)
        # create the restrictions of the projection map
        proj_maps = []
        for i in range(n):
            p_A = patches[i].ambient_space()
            H = Hom(patches[i], self)
            homvars = list(p_A.gens())[1:n]
            homvars.insert(i, 1)
            coords = [(p_A.gens()[0] - P[i])*homvars[j] + P[j] for j in range(n)]
            proj_maps.append(H(coords))
        return tuple([tuple(patches), tuple(t_maps), tuple(proj_maps)])

    def resolution_of_singularities(self, extend=False):
        r"""
        Return a nonsingular model for this affine curve created by blowing up its singular points.

        The nonsingular model is given as a collection of affine patches that cover it. If ``extend`` is ``False``
        and if the base field is a number field, or if the base field is a finite field, the model returned may have
        singularities with coordinates not contained in the base field. An error is returned if this curve is already
        nonsingular, or if it has no singular points over its base field. This curve must be irreducible, and must be
        defined over a number field or finite field.

        INPUT:

        - ``extend`` -- (default: False) specifies whether to extend the base field when necessary to find all
          singular points when this curve is defined over a number field. If ``extend`` is ``False``, then only
          singularities with coordinates in the base field of this curve will be resolved. However, setting
          ``extend`` to ``True`` will slow down computations.

        OUTPUT:

        - a tuple consisting of three elements. The first is a tuple of curves in affine space of
          the same dimension as the ambient space of this curve, which represent affine patches
          of the resolution of singularities. The second is a tuple of tuples such that the jth
          element of the ith tuple is the transition map from the ith patch to the jth patch. Lastly,
          the third element is a tuple consisting of birational maps from the patches back to the
          original curve that were created by composing the projection maps generated from the blow up
          computations. There the ith element will be a map from the ith patch.

        EXAMPLES::

            sage: A.<x,y> = AffineSpace(QQ, 2)
            sage: C = Curve([y^2 - x^3], A)
            sage: C.resolution_of_singularities()
            ((Affine Plane Curve over Rational Field defined by s1^2 - x,
              Affine Plane Curve over Rational Field defined by y*s0^3 - 1),
             ((Scheme endomorphism of Affine Plane Curve over Rational Field defined by s1^2 - x
                 Defn: Defined on coordinates by sending (x, s1) to
                       (x, s1), Scheme morphism:
                 From: Affine Plane Curve over Rational Field defined by s1^2 - x
                 To:   Affine Plane Curve over Rational Field defined by y*s0^3 - 1
                 Defn: Defined on coordinates by sending (x, s1) to
                       (x*s1, 1/s1)), (Scheme morphism:
                 From: Affine Plane Curve over Rational Field defined by y*s0^3 - 1
                 To:   Affine Plane Curve over Rational Field defined by s1^2 - x
                 Defn: Defined on coordinates by sending (y, s0) to
                       (y*s0, 1/s0),
               Scheme endomorphism of Affine Plane Curve over Rational Field defined by y*s0^3 - 1
                 Defn: Defined on coordinates by sending (y, s0) to
                       (y, s0))),
             (Scheme morphism:
                From: Affine Plane Curve over Rational Field defined by s1^2 - x
                To:   Affine Plane Curve over Rational Field defined by -x^3 + y^2
                Defn: Defined on coordinates by sending (x, s1) to
                      (x, x*s1), Scheme morphism:
                From: Affine Plane Curve over Rational Field defined by y*s0^3 - 1
                To:   Affine Plane Curve over Rational Field defined by -x^3 + y^2
                Defn: Defined on coordinates by sending (y, s0) to
                      (y*s0, y)))

        ::

            sage: set_verbose(-1)
            sage: K.<a> = QuadraticField(3)
            sage: A.<x,y> = AffineSpace(K, 2)
            sage: C = A.curve(x^4 + 2*x^2 + a*y^3 + 1)
            sage: C.resolution_of_singularities(extend=True)[0] # long time (2 seconds)
            (Affine Plane Curve over Number Field in a0 with defining polynomial y^4 - 4*y^2 + 16 defined by
            24*x^2*ss1^3 + 24*ss1^3 + (a0^3 - 8*a0),
             Affine Plane Curve over Number Field in a0 with defining polynomial y^4 - 4*y^2 + 16 defined by
             24*s1^2*ss0 + (a0^3 - 8*a0)*ss0^2 + (-6*a0^3)*s1,
             Affine Plane Curve over Number Field in a0 with defining polynomial y^4 - 4*y^2 + 16 defined by
             8*y^2*s0^4 + (4*a0^3)*y*s0^3 - 32*s0^2 + (a0^3 - 8*a0)*y)

        ::

            sage: A.<x,y,z> = AffineSpace(GF(5), 3)
            sage: C = Curve([y - x^3, (z - 2)^2 - y^3 - x^3], A)
            sage: R = C.resolution_of_singularities()
            sage: R[0]
            (Affine Curve over Finite Field of size 5 defined by x^2 - s1, s1^4 - x*s2^2 + s1, x*s1^3 - s2^2 + x,
             Affine Curve over Finite Field of size 5 defined by y*s2^2 - y^2 - 1, s2^4 - s0^3 - y^2 - 2, y*s0^3
             - s2^2 + y, Affine Curve over Finite Field of size 5 defined by s0^3*s1 + z*s1^3 + s1^4 - 2*s1^3 - 1,
             z*s0^3 + z*s1^3 - 2*s0^3 - 2*s1^3 - 1, z^2*s1^3 + z*s1^3 - s1^3 - z + s1 + 2)

        ::

            sage: A.<x,y,z,w> = AffineSpace(QQ, 4)
            sage: C = A.curve([((x - 2)^2 + y^2)^2 - (x - 2)^2 - y^2 + (x - 2)^3, z - y - 7, w - 4])
            sage: B = C.resolution_of_singularities()
            sage: B[0]
            (Affine Curve over Rational Field defined by s3, s1 - s2, x^2*s2^4 -
            4*x*s2^4 + 2*x^2*s2^2 + 4*s2^4 - 8*x*s2^2 + x^2 + 7*s2^2 - 3*x + 1,
             Affine Curve over Rational Field defined by s3, s2 - 1, y^2*s0^4 +
            2*y^2*s0^2 + y*s0^3 + y^2 - s0^2 - 1,
             Affine Curve over Rational Field defined by s3, s1 - 1, z^2*s0^4 -
            14*z*s0^4 + 2*z^2*s0^2 + z*s0^3 + 49*s0^4 - 28*z*s0^2 - 7*s0^3 + z^2 +
            97*s0^2 - 14*z + 48,
             Closed subscheme of Affine Space of dimension 4 over Rational Field
            defined by:
               1)

        ::

            sage: A.<x,y> = AffineSpace(QQ, 2)
            sage: C = Curve([y - x^2 + 1], A)
            sage: C.resolution_of_singularities()
            Traceback (most recent call last):
            ...
            TypeError: this curve is already nonsingular

        ::

            sage: A.<x,y> = AffineSpace(QQ, 2)
            sage: C = A.curve([(x^2 + y^2 - y - 2)*(y - x^2 + 2) + y^3])
            sage: C.resolution_of_singularities()
            Traceback (most recent call last):
            ...
            TypeError: this curve has no singular points over its base field. If
            working over a number field use extend=True
        """
        # helper function for extending the base field (in the case of working over a number field)
        def extension(self):
            F = self.base_ring()
            pts = self.change_ring(F.embeddings(QQbar)[0]).rational_points()
            L = [t for pt in pts for t in pt]
            K = number_field_elements_from_algebraics(L)[0]
            if is_RationalField(K):
                return F.embeddings(F)[0]
            else:
                if is_RationalField(F):
                    return F.embeddings(K)[0]
                else:
                    # make sure the defining polynomial variable names are the same for K, N
                    N = NumberField(K.defining_polynomial().parent()(F.defining_polynomial()), str(K.gen()))
                    return N.composite_fields(K, both_maps=True)[0][1]*F.embeddings(N)[0]
        # find the set of singular points of this curve
        # in the case that the base field is a number field, extend it as needed (if extend == True)
        C = self
        n = C.ambient_space().dimension_relative()
        if not self.is_irreducible():
            raise TypeError("this curve must be irreducible")
        if not (self.base_ring() in NumberFields() or self.base_ring() in FiniteFields()):
            raise NotImplementedError("this curve must be defined over either a number field or a finite field")
        if C.base_ring() in NumberFields() and extend:
            C = C.change_ring(extension(C.singular_subscheme()))
        H = End(C)
        placeholder = H(C.ambient_space().gens())
        # the list res holds the data for the patches of the resolution of singularities
        # each element is a list consisting of the curve defining the patch, a list
        # of the transition maps from that patch to the other patches, a projection
        # map from the patch to the original curve, and the set of singular points
        # of the patch
        res = [[C, [placeholder], placeholder, C.singular_points()]]
        if len(res[0][3]) == 0:
            if C.is_smooth():
                raise TypeError("this curve is already nonsingular")
            else:
                raise TypeError("this curve has no singular points over its base field. If working over"\
                                " a number field use extend=True")
        not_resolved = True
        t = 0
        # loop through the patches and blow up each until no patch has singular points
        while not_resolved:
            [BC, t_maps, pi, pts] = [res[t][0], res[t][1], res[t][2], res[t][3]]
            # check if there are any singular points in this patch
            if len(pts) == 0:
                t = t + 1
                if t == len(res):
                    not_resolved = False
                continue
            # the identity map should be replaced for each of the charts of the blow up
            t_maps.pop(t)
            # blow up pts[0]
            B = list(BC.blowup(pts[0]))
            B = [list(B[0]), [list(B[1][i]) for i in range(len(B[1]))], list(B[2])]
            # the t-th element of res will be replaced with the new data corresponding to the charts
            # of the blow up
            res.pop(t)
            # take out the transition maps from the other resolution patches to the t-th patch
            old_maps = [res[i][1].pop(t) for i in range(len(res))]
            patches_to_add = []
            # generate the needed data for each patch of the blow up
            for i in range(len(B[0])):
                # check if there are any singular points where this patch meets the exceptional divisor
                AA = AffineSpace(B[0][i].base_ring(), n - 1, 'x')
                coords = [pts[0][i]]
                coords.extend(list(AA.gens()))
                H = Hom(B[0][i].ambient_space().coordinate_ring(), AA.coordinate_ring())
                poly_hom = H(coords)
                X = AA.subscheme([poly_hom(f) for f in B[0][i].defining_polynomials()])
                # in the case of working over a number field, it might be necessary to extend the base
                # field in order to find all intersection points
                n_pts = []
                if B[0][i].base_ring() in NumberFields() and extend:
                    emb = extension(X)
                    X = X.change_ring(emb)
                    tmp_curve = B[0][i].change_ring(emb)
                    for pt in X.rational_points():
                        tmp_pt = tmp_curve([pts[0][i]] + list(pt))
                        if tmp_curve.is_singular(tmp_pt):
                            n_pts.append(tmp_pt)
                    # avoid needlessly extending the base field
                    if len(n_pts) > 0:
                        # coerce everything to the new base field
                        BC = BC.change_ring(emb)
                        t_maps = [t_maps[j].change_ring(emb) for j in range(len(t_maps))]
                        old_maps = [old_maps[j].change_ring(emb) for j in range(len(old_maps))]
                        pi = pi.change_ring(emb)
                        pts = [pt.change_ring(emb) for pt in pts]
                        # coerce the current blow up data
                        for j in range(len(B[0])):
                            B[0][j] = B[0][j].change_ring(emb)
                        for j in range(len(B[1])):
                            for k in range(len(B[1])):
                                B[1][j][k] = B[1][j][k].change_ring(emb)
                        for j in range(len(B[2])):
                            B[2][j] = B[2][j].change_ring(emb)
                        # coerce the other data in res
                        for j in range(len(res)):
                            res[j][0] = res[j][0].change_ring(emb)
                            for k in range(len(res[j][1])):
                                res[j][1][k] = res[j][1][k].change_ring(emb)
                            res[j][2].change_ring(emb)
                            for k in range(len(res[j][3])):
                                res[j][3][k] = res[j][3][k].change_ring(emb)
                else:
                    for pt in X.rational_points():
                        tmp_pt = B[0][i]([pts[0][i]] + list(pt))
                        if B[0][i].is_singular(tmp_pt):
                            n_pts.append(tmp_pt)
                b_data = [B[0][i]]
                # projection map and its inverse
                t_pi = B[2][i]
                coords = [(BC.ambient_space().gens()[j] - pts[0][j])/(BC.ambient_space().gens()[i] - pts[0][i]) for\
                          j in range(n)]
                coords.pop(i)
                coords.insert(0, BC.ambient_space().gens()[i])
                H = Hom(BC, B[0][i])
                t_pi_inv = H(coords)
                # compose the current transition maps from the original curve to the other patches
                # with the projection map
                L = list(t_maps)
                for j in range(len(t_maps)):
                    L[j] = L[j]*t_pi
                for j in range(len(B[1][i])):
                    L.insert(t + j, B[1][i][j])
                b_data.append(L)
                # update transition maps of each other element of res
                for j in range(len(res)):
                    new_t_map = t_pi_inv*old_maps[j]
                    res[j][1].insert(t + i, new_t_map)
                # create the projection map
                b_data.append(pi*t_pi)
                # singular points
                # translate the singular points of the parent patch (other than that which was the center of the
                # blow up) by the inverse of the first projection map
                for j in range(1, len(pts)):
                    # make sure this point is in this chart before attempting to map it
                    try:
                        n_pts.append(t_pi_inv(BC(pts[j])))
                    except (TypeError, ZeroDivisionError):
                        pass
                b_data.append(n_pts)
                patches_to_add.append(b_data)
            for i in range(len(patches_to_add)):
                res.insert(t + i, patches_to_add[i])
            t = 0
        patches = [res[i][0] for i in range(len(res))]
        t_maps = [tuple(res[i][1]) for i in range(len(res))]
        p_maps = [res[i][2] for i in range(len(res))]
        return tuple([tuple(patches), tuple(t_maps), tuple(p_maps)])

class AffinePlaneCurve(AffineCurve):

    _point = point.AffinePlaneCurvePoint_field

    def __init__(self, A, f):
        r"""
        Initialization function.

        EXAMPLES::

            sage: A.<x,y> = AffineSpace(QQ, 2)
            sage: C = Curve([x^3 - y^2], A); C
            Affine Plane Curve over Rational Field defined by x^3 - y^2

        ::

            sage: A.<x,y> = AffineSpace(CC, 2)
            sage: C = Curve([y^2 + x^2], A); C
            Affine Plane Curve over Complex Field with 53 bits of precision defined
            by x^2 + y^2
        """
        if not (is_AffineSpace(A) and A.dimension != 2):
            raise TypeError("Argument A (= %s) must be an affine plane." % A)
        Curve_generic.__init__(self, A, [f])

    def _repr_type(self):
        r"""
        Return a string representation of the type of this curve.

        EXAMPLES::

            sage: A.<x,y> = AffineSpace(QQ, 2)
            sage: C = Curve([y - 7/2*x^5 + x - 3], A)
            sage: C._repr_type()
            'Affine Plane'
        """
        return "Affine Plane"

    def divisor_of_function(self, r):
        """
        Return the divisor of a function on a curve.

        INPUT: r is a rational function on X

        OUTPUT:


        -  ``list`` - The divisor of r represented as a list of
           coefficients and points. (TODO: This will change to a more
           structural output in the future.)


        EXAMPLES::

            sage: F = GF(5)
            sage: P2 = AffineSpace(2, F, names = 'xy')
            sage: R = P2.coordinate_ring()
            sage: x, y = R.gens()
            sage: f = y^2 - x^9 - x
            sage: C = Curve(f)
            sage: K = FractionField(R)
            sage: r = 1/x
            sage: C.divisor_of_function(r)     # todo: not implemented (broken)
                  [[-1, (0, 0, 1)]]
            sage: r = 1/x^3
            sage: C.divisor_of_function(r)     # todo: not implemented (broken)
                  [[-3, (0, 0, 1)]]
        """
        F = self.base_ring()
        f = self.defining_polynomial()
        pts = self.places_on_curve()
        R = f.parent()
        x, y = R.gens()
        R0 = PolynomialRing(F, 3, names=[str(x), str(y), "t"])
        vars0 = R0.gens()
        t = vars0[2]
        divf = []
        for pt0 in pts:
            if pt0[2] != F(0):
                lcs = self.local_coordinates(pt0, 5)
                yt = lcs[1]
                xt = lcs[0]
                ldg = degree_lowest_rational_function(r(xt, yt), t)
                if ldg != 0:
                    divf.append([ldg, pt0])
        return divf

    def local_coordinates(self, pt, n):
        r"""
        Return local coordinates to precision n at the given point.

        Behaviour is flaky - some choices of `n` are worst that
        others.


        INPUT:


        -  ``pt`` - an F-rational point on X which is not a
           point of ramification for the projection (x,y) - x.

        -  ``n`` - the number of terms desired


        OUTPUT: x = x0 + t y = y0 + power series in t

        EXAMPLES::

            sage: F = GF(5)
            sage: pt = (2,3)
            sage: R = PolynomialRing(F,2, names = ['x','y'])
            sage: x,y = R.gens()
            sage: f = y^2-x^9-x
            sage: C = Curve(f)
            sage: C.local_coordinates(pt, 9)
            [t + 2, -2*t^12 - 2*t^11 + 2*t^9 + t^8 - 2*t^7 - 2*t^6 - 2*t^4 + t^3 - 2*t^2 - 2]
        """
        f = self.defining_polynomial()
        R = f.parent()
        F = self.base_ring()
        p = F.characteristic()
        x0 = F(pt[0])
        y0 = F(pt[1])
        astr = ["a"+str(i) for i in range(1,2*n)]
        x,y = R.gens()
        R0 = PolynomialRing(F,2*n+2,names = [str(x),str(y),"t"]+astr)
        vars0 = R0.gens()
        t = vars0[2]
        yt = y0*t**0+add([vars0[i]*t**(i-2) for i in range(3,2*n+2)])
        xt = x0+t
        ft = f(xt,yt)
        S = singular
        S.eval('ring s = '+str(p)+','+str(R0.gens())+',lp;')
        S.eval('poly f = '+str(ft) + ';')
        c = S('coeffs(%s, t)' % ft)
        N = int(c.size())
        b = ','.join("%s[%s,1]" % (c.name(), i) for i in range(2, N//2-4))
        cmd = 'ideal I = ' + b
        S.eval(cmd)
        S.eval('short=0')    # print using *'s and ^'s.
        c = S.eval('slimgb(I)')
        d = c.split("=")
        d = d[1:]
        d[len(d)-1] += "\n"
        e = [xx[:xx.index("\n")] for xx in d]
        vals = []
        for x in e:
            for y in vars0:
                if str(y) in x:
                    if len(x.replace(str(y),"")) != 0:
                        i = x.find("-")
                        if i>0:
                            vals.append([eval(x[1:i]),x[:i],F(eval(x[i+1:]))])
                        i = x.find("+")
                        if i>0:
                            vals.append([eval(x[1:i]),x[:i],-F(eval(x[i+1:]))])
                    else:
                        vals.append([eval(str(y)[1:]),str(y),F(0)])
        vals.sort()
        return [x0 + t, y0 + add(v[2] * t**(j+1) for j, v in enumerate(vals))]

    def plot(self, *args, **kwds):
        r"""
        Plot the real points on this affine plane curve.

        INPUT:


        -  ``self`` - an affine plane curve

        -  ``*args`` - optional tuples (variable, minimum, maximum) for
           plotting dimensions

        -  ``**kwds`` - optional keyword arguments passed on to
           ``implicit_plot``


        EXAMPLES:

        A cuspidal curve::

            sage: R.<x, y> = QQ[]
            sage: C = Curve(x^3 - y^2)
            sage: C.plot()
            Graphics object consisting of 1 graphics primitive

        A 5-nodal curve of degree 11.  This example also illustrates
        some of the optional arguments::

            sage: R.<x, y> = ZZ[]
            sage: C = Curve(32*x^2 - 2097152*y^11 + 1441792*y^9 - 360448*y^7 + 39424*y^5 - 1760*y^3 + 22*y - 1)
            sage: C.plot((x, -1, 1), (y, -1, 1), plot_points=400)
            Graphics object consisting of 1 graphics primitive

        A line over `\mathbf{RR}`::

            sage: R.<x, y> = RR[]
            sage: C = Curve(R(y - sqrt(2)*x))
            sage: C.plot()
            Graphics object consisting of 1 graphics primitive
        """
        I = self.defining_ideal()
        return I.plot(*args, **kwds)

    def is_transverse(self, C, P):
        r"""
        Return whether the intersection of this curve with the curve ``C`` at the point ``P`` is transverse.

        The intersection at ``P`` is transverse if ``P`` is a nonsingular point of both curves, and if the
        tangents of the curves at ``P`` are distinct.

        INPUT:

        - ``C`` -- a curve in the ambient space of this curve.

        - ``P`` -- a point in the intersection of both curves.

        OUTPUT: Boolean.

        EXAMPLES::

            sage: A.<x,y> = AffineSpace(QQ, 2)
            sage: C = Curve([x^2 + y^2 - 1], A)
            sage: D = Curve([x - 1], A)
            sage: Q = A([1,0])
            sage: C.is_transverse(D, Q)
            False

        ::

            sage: R.<a> = QQ[]
            sage: K.<b> = NumberField(a^3 + 2)
            sage: A.<x,y> = AffineSpace(K, 2)
            sage: C = A.curve([x*y])
            sage: D = A.curve([y - b*x])
            sage: Q = A([0,0])
            sage: C.is_transverse(D, Q)
            False

        ::

            sage: A.<x,y> = AffineSpace(QQ, 2)
            sage: C = Curve([y - x^3], A)
            sage: D = Curve([y + x], A)
            sage: Q = A([0,0])
            sage: C.is_transverse(D, Q)
            True
        """
        if not self.intersects_at(C, P):
            raise TypeError("(=%s) must be a point in the intersection of (=%s) and this curve" % (P, C))
        if self.is_singular(P) or C.is_singular(P):
            return False

        # there is only one tangent at a nonsingular point of a plane curve
        return not self.tangents(P)[0] == C.tangents(P)[0]

    def multiplicity(self, P):
        r"""
        Return the multiplicity of this affine plane curve at the point ``P``.

        In the special case of affine plane curves, the multiplicity of an affine
        plane curve at the point (0,0) can be computed as the minimum of the degrees
        of the homogeneous components of its defining polynomial. To compute the
        multiplicity of a different point, a linear change of coordinates is used.

        This curve must be defined over a field. An error if raised if ``P`` is
        not a point on this curve.

        INPUT:

        - ``P`` -- a point in the ambient space of this curve.

        OUTPUT:

        An integer.

        EXAMPLES::

            sage: A.<x,y> = AffineSpace(QQ, 2)
            sage: C = Curve([y^2 - x^3], A)
            sage: Q1 = A([1,1])
            sage: C.multiplicity(Q1)
            1
            sage: Q2 = A([0,0])
            sage: C.multiplicity(Q2)
            2

        ::

            sage: A.<x,y> = AffineSpace(QQbar,2)
            sage: C = Curve([-x^7 + (-7)*x^6 + y^6 + (-21)*x^5 + 12*y^5 + (-35)*x^4 + 60*y^4 +\
            (-35)*x^3 + 160*y^3 + (-21)*x^2 + 240*y^2 + (-7)*x + 192*y + 63], A)
            sage: Q = A([-1,-2])
            sage: C.multiplicity(Q)
            6

        ::

            sage: A.<x,y> = AffineSpace(QQ, 2)
            sage: C = A.curve([y^3 - x^3 + x^6])
            sage: Q = A([1,1])
            sage: C.multiplicity(Q)
            Traceback (most recent call last):
            ...
            TypeError: (=(1, 1)) is not a point on (=Affine Plane Curve over
            Rational Field defined by x^6 - x^3 + y^3)
        """
        if not self.base_ring() in Fields():
            raise TypeError("curve must be defined over a field")

        # Check whether P is a point on this curve
        try:
            P = self(P)
        except TypeError:
            raise TypeError("(=%s) is not a point on (=%s)" % (P, self))

        # Apply a linear change of coordinates to self so that P becomes (0,0)
        AA = self.ambient_space()
        f = self.defining_polynomials()[0](AA.gens()[0] + P[0], AA.gens()[1] + P[1])

        # Compute the multiplicity of the new curve at (0,0), which is the minimum of the degrees of its
        # nonzero terms
        return min([g.degree() for g in f.monomials()])

    def tangents(self, P, factor=True):
        r"""
        Return the tangents of this affine plane curve at the point ``P``.

        The point ``P`` must be a point on this curve.

        INPUT:

        - ``P`` -- a point on this curve.

        - ``factor`` -- (default: True) whether to attempt computing the polynomials of the individual tangent
          lines over the base field of this curve, or to just return the polynomial corresponding to the union
          of the tangent lines (which requires fewer computations).

        OUTPUT:

        - a list of polynomials in the coordinate ring of the ambient space of this curve.

        EXAMPLES::

            sage: set_verbose(-1)
            sage: A.<x,y> = AffineSpace(QQbar, 2)
            sage: C = Curve([x^5*y^3 + 2*x^4*y^4 + x^3*y^5 + 3*x^4*y^3 + 6*x^3*y^4 + 3*x^2*y^5\
            + 3*x^3*y^3 + 6*x^2*y^4 + 3*x*y^5 + x^5 + 10*x^4*y + 40*x^3*y^2 + 81*x^2*y^3 + 82*x*y^4\
            + 33*y^5], A)
            sage: Q = A([0,0])
            sage: C.tangents(Q)
            [x + 3.425299577684700?*y, x + (1.949159013086856? + 1.179307909383728?*I)*y,
            x + (1.949159013086856? - 1.179307909383728?*I)*y, x + (1.338191198070795? + 0.2560234251008043?*I)*y,
            x + (1.338191198070795? - 0.2560234251008043?*I)*y]
            sage: C.tangents(Q, factor=False)
            [120*x^5 + 1200*x^4*y + 4800*x^3*y^2 + 9720*x^2*y^3 + 9840*x*y^4 + 3960*y^5]

        ::

            sage: R.<a> = QQ[]
            sage: K.<b> = NumberField(a^2 - 3)
            sage: A.<x,y> = AffineSpace(K, 2)
            sage: C = Curve([(x^2 + y^2 - 2*x)^2 - x^2 - y^2], A)
            sage: Q = A([0,0])
            sage: C.tangents(Q)
            [x + (-1/3*b)*y, x + (1/3*b)*y]

        ::

            sage: A.<x,y> = AffineSpace(QQ, 2)
            sage: C = A.curve([y^2 - x^3 - x^2])
            sage: Q = A([0,0])
            sage: C.tangents(Q)
            [x - y, x + y]

        ::

            sage: A.<x,y> = AffineSpace(QQ, 2)
            sage: C = A.curve([y*x - x^4 + 2*x^2])
            sage: Q = A([1,1])
            sage: C.tangents(Q)
            Traceback (most recent call last):
            ...
            TypeError: (=(1, 1)) is not a point on (=Affine Plane Curve over
            Rational Field defined by -x^4 + 2*x^2 + x*y)
        """
        r = self.multiplicity(P)
        f = self.defining_polynomial()
        # move P to (0,0)
        vars = self.ambient_space().gens()
        coords = [vars[0] + P[0], vars[1] + P[1]]
        f = f(coords)
        coords = [vars[0] - P[0], vars[1] - P[1]] # coords to change back with
        deriv = [f.derivative(vars[0],i).derivative(vars[1], r-i)([0,0]) for i in range(r+1)]
        T = sum([binomial(r,i)*deriv[i]*(vars[0])**i*(vars[1])**(r-i) for i in range(r+1)])
        if not factor:
            return [T(coords)]
        if self.base_ring() == QQbar:
            fact = []
            # first add tangents corresponding to vars[0], vars[1] if they divide T
            t = min([e[0] for e in T.exponents()])
            # vars[0] divides T
            if t > 0:
                fact.append(vars[0])
                # divide T by that power of vars[0]
                T = self.ambient_space().coordinate_ring()(dict([((v[0] - t,v[1]), h) for (v,h) in T.dict().items()]))
            t = min([e[1] for e in T.exponents()])
            # vars[1] divides T
            if t > 0:
                fact.append(vars[1])
                # divide T by that power of vars[1]
                T = self.ambient_space().coordinate_ring()(dict([((v[0],v[1] - t), h) for (v,h) in T.dict().items()]))
            # T is homogeneous in var[0], var[1] if nonconstant, so dehomogenize
            if not T in self.base_ring():
                if T.degree(vars[0]) > 0:
                    T = T(vars[0], 1)
                    roots = T.univariate_polynomial().roots()
                    fact.extend([vars[0] - roots[i][0]*vars[1] for i in range(len(roots))])
                else:
                    T = T(1, vars[1])
                    roots = T.univariate_polynomial().roots()
                    fact.extend([vars[1] - roots[i][0]*vars[0] for i in range(len(roots))])
            return [ff(coords) for ff in fact]
        else:
            return [l[0](coords) for l in T.factor()]

    def is_ordinary_singularity(self, P):
        r"""
        Return whether the singular point ``P`` of this affine plane curve is an ordinary singularity.

        The point ``P`` is an ordinary singularity of this curve if it is a singular point, and
        if the tangents of this curve at ``P`` are distinct.

        INPUT:

        - ``P`` -- a point on this curve.

        OUTPUT:

        - Boolean. True or False depending on whether ``P`` is or is not an ordinary singularity of this
          curve, respectively. An error is raised if ``P`` is not a singular point of this curve.

        EXAMPLES::

            sage: A.<x,y> = AffineSpace(QQ, 2)
            sage: C = Curve([y^2 - x^3], A)
            sage: Q = A([0,0])
            sage: C.is_ordinary_singularity(Q)
            False

        ::

            sage: R.<a> = QQ[]
            sage: K.<b> = NumberField(a^2 - 3)
            sage: A.<x,y> = AffineSpace(K, 2)
            sage: C = Curve([(x^2 + y^2 - 2*x)^2 - x^2 - y^2], A)
            sage: Q = A([0,0])
            sage: C.is_ordinary_singularity(Q)
            True

        ::

            sage: A.<x,y> = AffineSpace(QQ, 2)
            sage: C = A.curve([x^2*y - y^2*x + y^2 + x^3])
            sage: Q = A([-1,-1])
            sage: C.is_ordinary_singularity(Q)
            Traceback (most recent call last):
            ...
            TypeError: (=(-1, -1)) is not a singular point of (=Affine Plane Curve
            over Rational Field defined by x^3 + x^2*y - x*y^2 + y^2)
        """
        r = self.multiplicity(P)
        if r < 2:
            raise TypeError("(=%s) is not a singular point of (=%s)" % (P,self))

        T = self.tangents(P, factor=False)[0]
        vars = self.ambient_space().gens()

        # use resultants to determine if there is a higher multiplicity tangent
        if T.degree(vars[0]) > 0:
            return T.resultant(T.derivative(vars[0]), vars[0]) != 0
        else:
            return T.resultant(T.derivative(vars[1]), vars[1]) != 0

    def rational_parameterization(self):
        r"""
        Return a rational parameterization of this curve.

        This curve must have rational coefficients and be absolutely irreducible (i.e. irreducible
        over the algebraic closure of the rational field). The curve must also be rational (have
        geometric genus zero).

        The rational parameterization may have coefficients in a quadratic extension of the rational
        field.

        OUTPUT:

        - a birational map between `\mathbb{A}^{1}` and this curve, given as a scheme morphism.

        EXAMPLES::

            sage: A.<x,y> = AffineSpace(QQ, 2)
            sage: C = Curve([y^2 - x], A)
            sage: C.rational_parameterization()
            Scheme morphism:
              From: Affine Space of dimension 1 over Rational Field
              To:   Affine Plane Curve over Rational Field defined by y^2 - x
              Defn: Defined on coordinates by sending (t) to
                    (t^2, t)

        ::

            sage: A.<x,y> = AffineSpace(QQ, 2)
            sage: C = Curve([(x^2 + y^2 - 2*x)^2 - x^2 - y^2], A)
            sage: C.rational_parameterization()
            Scheme morphism:
              From: Affine Space of dimension 1 over Rational Field
              To:   Affine Plane Curve over Rational Field defined by x^4 +
            2*x^2*y^2 + y^4 - 4*x^3 - 4*x*y^2 + 3*x^2 - y^2
              Defn: Defined on coordinates by sending (t) to
                    ((-12*t^4 + 6*t^3 + 4*t^2 - 2*t)/(-25*t^4 + 40*t^3 - 26*t^2 +
            8*t - 1), (-9*t^4 + 12*t^3 - 4*t + 1)/(-25*t^4 + 40*t^3 - 26*t^2 + 8*t - 1))

        ::

            sage: A.<x,y> = AffineSpace(QQ, 2)
            sage: C = Curve([x^2 + y^2 + 7], A)
            sage: C.rational_parameterization()
            Scheme morphism:
              From: Affine Space of dimension 1 over Number Field in a with defining polynomial a^2 + 7
              To:   Affine Plane Curve over Number Field in a with defining
            polynomial a^2 + 7 defined by x^2 + y^2 + 7
              Defn: Defined on coordinates by sending (t) to
                    ((-7*t^2 + 7)/((-a)*t^2 + (-a)), 14*t/((-a)*t^2 + (-a)))
        """
        para = self.projective_closure(i=0).rational_parameterization().defining_polynomials()
        # these polynomials are homogeneous in two indeterminants, so dehomogenize wrt one of the variables
        R = para[0].parent()
        A_line = AffineSpace(R.base_ring(), 1, 't')
        para = [A_line.coordinate_ring()(para[i].substitute({R.gens()[0]: 1})) for i in range(3)]
        C = self.change_ring(R.base_ring())
        # because of the parameter i=0, the projective closure is constructed with respect to the
        # affine patch corresponding to the first coordinate being nonzero. Thus para[0] will not be
        # the zero polynomial, and dehomogenization won't change this
        H = Hom(A_line, C)
        return H([para[1]/para[0], para[2]/para[0]])

    def fundamental_group(self):
        r"""
        Return a presentation of the fundamental group of the complement
        of ``self``.

        .. NOTE::

            The curve must be defined over the rationals or a number field
            with an embedding over `\QQbar`.

        EXAMPLES::

            sage: A.<x,y> = AffineSpace(QQ, 2)
            sage: C = A.curve(y^2 - x^3 - x^2)
            sage: C.fundamental_group() # optional - sirocco
            Finitely presented group < x0 |  >

        In the case of number fields, they need to have an embedding
        to the algebraic field::

            sage: a = QQ[x](x^2+5).roots(QQbar)[0][0]
            sage: F = NumberField(a.minpoly(), 'a', embedding=a)
            sage: F.inject_variables()
            Defining a
            sage: A.<x,y> = AffineSpace(F, 2)
            sage: C = A.curve(y^2 - a*x^3 - x^2)
            sage: C.fundamental_group() # optional - sirocco
            Finitely presented group < x0 |  >

        .. WARNING::

            This functionality requires the sirocco package to be installed.
        """
        from sage.schemes.curves.zariski_vankampen import fundamental_group
        F = self.base_ring()
        from sage.rings.qqbar import QQbar
        if QQbar.coerce_map_from(F) is None:
            raise NotImplementedError("the base field must have an embedding"
                                      " to the algebraic field")
        f = self.defining_polynomial()
        return fundamental_group(f, projective=False)

    def riemann_surface(self,**kwargs):
        r"""Return the complex riemann surface determined by this curve

        OUTPUT:

         - RiemannSurface object

        EXAMPLES::

            sage: R.<x,y>=QQ[]
            sage: C=Curve(x^3+3*y^3+5)
            sage: C.riemann_surface()
            Riemann surface defined by polynomial f = x^3 + 3*y^3 + 5 = 0, with 53 bits of precision

        """
        from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
        return RiemannSurface(self.defining_polynomial(),**kwargs)


class AffinePlaneCurve_finite_field(AffinePlaneCurve):

    _point = point.AffinePlaneCurvePoint_finite_field

    def rational_points(self, algorithm="enum"):
        r"""
        Return sorted list of all rational points on this curve.

        Use *very* naive point enumeration to find all rational points on
        this curve over a finite field.

        EXAMPLES::

            sage: A.<x,y> = AffineSpace(2,GF(9,'a'))
            sage: C = Curve(x^2 + y^2 - 1)
            sage: C
            Affine Plane Curve over Finite Field in a of size 3^2 defined by x^2 + y^2 - 1
            sage: C.rational_points()
            [(0, 1), (0, 2), (1, 0), (2, 0), (a + 1, a + 1), (a + 1, 2*a + 2), (2*a + 2, a + 1), (2*a + 2, 2*a + 2)]
        """
        f = self.defining_polynomial()
        R = f.parent()
        K = R.base_ring()
        points = []
        for x in K:
            for y in K:
                if f(x,y) == 0:
                    points.append(self((x,y)))
        points.sort()
        return points


class AffinePlaneCurve_prime_finite_field(AffinePlaneCurve_finite_field):
    # CHECK WHAT ASSUMPTIONS ARE MADE REGARDING AFFINE VS. PROJECTIVE MODELS!!!
    # THIS IS VERY DIRTY STILL -- NO DATASTRUCTURES FOR DIVISORS.

    def riemann_roch_basis(self, D):
        r"""
        Interfaces with Singular's BrillNoether command.

        INPUT:


        -  ``self`` - a plane curve defined by a polynomial eqn f(x,y)
           = 0 over a prime finite field F = GF(p) in 2 variables x,y
           representing a curve X: f(x,y) = 0 having n F-rational
           points (see the Sage function places_on_curve)

        -  ``D`` - an n-tuple of integers
           `(d1, ..., dn)` representing the divisor
           `Div = d1*P1+...+dn*Pn`, where
           `X(F) = \{P1,...,Pn\}`.
           *The ordering is that dictated by places_on_curve.*


        OUTPUT: basis of L(Div)

        EXAMPLES::

            sage: R = PolynomialRing(GF(5),2,names = ["x","y"])
            sage: x, y = R.gens()
            sage: f = y^2 - x^9 - x
            sage: C = Curve(f)
            sage: D = [6,0,0,0,0,0]
            sage: C.riemann_roch_basis(D)
            [1, (-x*z^5 + y^2*z^4)/x^6, (-x*z^6 + y^2*z^5)/x^7, (-x*z^7 + y^2*z^6)/x^8]
        """
        f = self.defining_polynomial()
        gens = f.parent().gens()
        p = self.base_ring().characteristic()
        G = singular(','.join(str(x) for x in D), type='intvec')

        singular.lib('brnoeth')
        singular.ring(p, gens, 'lp')

        X = singular(f).Adj_div()
        P = singular.NSplaces(1, X)
        T = P[1][2]
        T.set_ring()  # necessary

        return [g[1].sage() / g[2].sage() for g in G.BrillNoether(P)]

    def rational_points(self, algorithm="enum"):
        r"""
        Return sorted list of all rational points on this curve.

        INPUT:


        -  ``algorithm`` - string:

           +  ``'enum'`` - straightforward enumeration

           +  ``'bn'`` - via Singular's Brill-Noether package.

           +  ``'all'`` - use all implemented algorithms and
              verify that they give the same answer, then return it


        .. note::

           The Brill-Noether package does not always work. When it
           fails a RuntimeError exception is raised.

        EXAMPLES::

            sage: x, y = (GF(5)['x,y']).gens()
            sage: f = y^2 - x^9 - x
            sage: C = Curve(f); C
            Affine Plane Curve over Finite Field of size 5 defined by -x^9 + y^2 - x
            sage: C.rational_points(algorithm='bn')
            [(0, 0), (2, 2), (2, 3), (3, 1), (3, 4)]
            sage: C = Curve(x - y + 1)
            sage: C.rational_points()
            [(0, 1), (1, 2), (2, 3), (3, 4), (4, 0)]

        We compare Brill-Noether and enumeration::

            sage: x, y = (GF(17)['x,y']).gens()
            sage: C = Curve(x^2 + y^5 + x*y - 19)
            sage: v = C.rational_points(algorithm='bn')
            sage: w = C.rational_points(algorithm='enum')
            sage: len(v)
            20
            sage: v == w
            True
        """
        if algorithm == "enum":

            return AffinePlaneCurve_finite_field.rational_points(self, algorithm="enum")

        elif algorithm == "bn":
            f = self.defining_polynomial()._singular_()
            singular = f.parent()
            singular.lib('brnoeth')
            try:
                X1 = f.Adj_div()
            except (TypeError, RuntimeError) as s:
                raise RuntimeError(str(s) + "\n\n ** Unable to use the Brill-Noether Singular package to compute all points (see above).")

            X2 = singular.NSplaces(1, X1)
            R = X2[5][1][1]
            singular.set_ring(R)

            # We use sage_flattened_str_list since iterating through
            # the entire list through the sage/singular interface directly
            # would involve hundreds of calls to singular, and timing issues
            # with the expect interface could crop up.  Also, this is vastly
            # faster (and more robust).
            v = singular('POINTS').sage_flattened_str_list()
            pnts = [self(int(v[3*i]), int(v[3*i+1]))
                    for i in range(len(v)//3) if int(v[3*i+2])]
            # remove multiple points
            return sorted(set(pnts))

        elif algorithm == "all":

            S_enum = self.rational_points(algorithm = "enum")
            S_bn = self.rational_points(algorithm = "bn")
            if S_enum != S_bn:
                raise RuntimeError("Bug in rational_points -- different algorithms give different answers for curve %s!" % self)
            return S_enum

        else:
            raise ValueError("No algorithm '%s' known" % algorithm)
