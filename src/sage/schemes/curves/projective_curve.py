"""
Projective curves.

EXAMPLES:

We can construct curves in either a projective plane::

    sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
    sage: C = Curve([y*z^2 - x^3], P); C
    Projective Plane Curve over Rational Field defined by -x^3 + y*z^2

or in higher dimensional projective spaces::

    sage: P.<x,y,z,w> = ProjectiveSpace(QQ, 3)
    sage: C = Curve([y*w^3 - x^4, z*w^3 - x^4], P); C
    Projective Curve over Rational Field defined by -x^4 + y*w^3, -x^4 + z*w^3

AUTHORS:

- William Stein (2005-11-13)

- David Joyner (2005-11-13)

- David Kohel (2006-01)

- Moritz Minzlaff (2010-11)

- Grayson Jorgenson (2016-8)
"""

#*****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import division
from __future__ import absolute_import

from sage.categories.fields import Fields
from sage.categories.number_fields import NumberFields
from sage.categories.homset import Hom, End
from sage.interfaces.all import singular
from sage.matrix.constructor import matrix
from sage.misc.all import add, sage_eval
from sage.rings.all import degree_lowest_rational_function
from sage.rings.number_field.number_field import NumberField
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.qqbar import (number_field_elements_from_algebraics,
                              QQbar)
from sage.rings.rational_field import is_RationalField
from sage.schemes.affine.affine_space import AffineSpace
from sage.schemes.projective.projective_space import ProjectiveSpace, is_ProjectiveSpace

from . import point

from sage.schemes.generic.algebraic_scheme import AlgebraicScheme_subscheme_projective
from sage.schemes.projective.projective_space import (is_ProjectiveSpace,
                                                      ProjectiveSpace)

from .curve import Curve_generic


class ProjectiveCurve(Curve_generic, AlgebraicScheme_subscheme_projective):

    _point = point.ProjectiveCurvePoint_field

    def _repr_type(self):
        r"""
        Return a string representation of the type of this curve.

        EXAMPLES::

            sage: P.<x,y,z,w> = ProjectiveSpace(QQ, 3)
            sage: C = Curve([y^3 - z^3 - w^3, z*x^3 - y^4])
            sage: C._repr_type()
            'Projective'
        """
        return "Projective"

    def __init__(self, A, X):
        r"""
        Initialization function.

        EXAMPLES::

            sage: P.<x,y,z,w,u> = ProjectiveSpace(GF(7), 4)
            sage: C = Curve([y*u^2 - x^3, z*u^2 - x^3, w*u^2 - x^3, y^3 - x^3], P); C
            Projective Curve over Finite Field of size 7 defined by -x^3 + y*u^2,
            -x^3 + z*u^2, -x^3 + w*u^2, -x^3 + y^3

        ::

            sage: K.<u> = CyclotomicField(11)
            sage: P.<x,y,z,w> = ProjectiveSpace(K, 3)
            sage: C = Curve([y*w - u*z^2 - x^2, x*w - 3*u^2*z*w], P); C
            Projective Curve over Cyclotomic Field of order 11 and degree 10 defined
            by -x^2 + (-u)*z^2 + y*w, x*w + (-3*u^2)*z*w
        """
        if not is_ProjectiveSpace(A):
            raise TypeError("A (=%s) must be a projective space"%A)
        Curve_generic.__init__(self, A, X)
        d = self.dimension()
        if d != 1:
            raise ValueError("defining equations (=%s) define a scheme of dimension %s != 1"%(X,d))

    def affine_patch(self, i, AA=None):
        r"""
        Return the i-th affine patch of this projective curve.

        INPUT:

        - ``i`` -- affine coordinate chart of the projective ambient space of this curve to compute affine patch
          with respect to.

        - ``AA`` -- (default: None) ambient affine space, this is constructed if it is not given.

        OUTPUT:

        - a curve in affine space.

        EXAMPLES::

            sage: P.<x,y,z,w> = ProjectiveSpace(CC, 3)
            sage: C = Curve([y*z - x^2, w^2 - x*y], P)
            sage: C.affine_patch(0)
            Affine Curve over Complex Field with 53 bits of precision defined by
            x0*x1 - 1.00000000000000, x2^2 - x0

        ::

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: C = Curve(x^3 - x^2*y + y^3 - x^2*z, P)
            sage: C.affine_patch(1)
            Affine Plane Curve over Rational Field defined by x0^3 - x0^2*x1 - x0^2 + 1

        ::

            sage: A.<x,y> = AffineSpace(QQ, 2)
            sage: P.<u,v,w> = ProjectiveSpace(QQ, 2)
            sage: C = Curve([u^2 - v^2], P)
            sage: C.affine_patch(1, A).ambient_space() == A
            True
        """
        from .constructor import Curve
        return Curve(AlgebraicScheme_subscheme_projective.affine_patch(self, i, AA))

    def projection(self, P=None, PS=None):
        r"""
        Return a projection of this curve into projective space of dimension one less than the dimension of
        the ambient space of this curve.

        This curve must not already be a plane curve. Over finite fields, if this curve contains all points
        in its ambient space, then an error will be returned.

        INPUT:

        - ``P`` -- (default: None) a point not on this curve that will be used to define the projection map;
          this is constructed if not specified.

        - ``PS`` -- (default: None) the projective space the projected curve will be defined in. This space must
          be defined over the same base ring as this curve, and must have dimension one less than that of the
          ambient space of this curve. This space will be constructed if not specified.

        OUTPUT:

        - a tuple consisting of two elements: a scheme morphism from this curve into a projective space of
          dimension one less than that of the ambient space of this curve, and the projective curve that
          is the image of that morphism.

        EXAMPLES::

            sage: K.<a> = CyclotomicField(3)
            sage: P.<x,y,z,w> = ProjectiveSpace(K, 3)
            sage: C = Curve([y*w - x^2, z*w^2 - a*x^3], P)
            sage: L.<a,b,c> = ProjectiveSpace(K, 2)
            sage: proj1 = C.projection(PS=L)
            sage: proj1
            (Scheme morphism:
               From: Projective Curve over Cyclotomic Field of order 3 and degree 2
            defined by -x^2 + y*w, (-a)*x^3 + z*w^2
               To:   Projective Space of dimension 2 over Cyclotomic Field of order
            3 and degree 2
               Defn: Defined on coordinates by sending (x : y : z : w) to
                     (x : y : -z + w),
             Projective Plane Curve over Cyclotomic Field of order 3 and degree 2
            defined by a^6 + (-a)*a^3*b^3 - a^4*b*c)
            sage: proj1[1].ambient_space() is L
            True
            sage: proj2 = C.projection()
            sage: proj2[1].ambient_space() is L
            False

        ::

            sage: P.<x,y,z,w,a,b,c> = ProjectiveSpace(QQ, 6)
            sage: C = Curve([y - x, z - a - b, w^2 - c^2, z - x - a, x^2 - w*z], P)
            sage: C.projection()
            (Scheme morphism:
               From: Projective Curve over Rational Field defined by -x + y, z - a -
            b, w^2 - c^2, -x + z - a, x^2 - z*w
               To:   Projective Space of dimension 5 over Rational Field
               Defn: Defined on coordinates by sending (x : y : z : w : a : b : c)
            to
                     (x : y : -z + w : a : b : c),
             Projective Curve over Rational Field defined by x1 - x4, x0 - x4, x2*x3
            + x3^2 + x2*x4 + 2*x3*x4, x2^2 - x3^2 - 2*x3*x4 + x4^2 - x5^2, x2*x4^2 +
            x3*x4^2 + x4^3 - x3*x5^2 - x4*x5^2, x4^4 - x3^2*x5^2 - 2*x3*x4*x5^2 -
            x4^2*x5^2)

        ::

            sage: P.<x,y,z,w> = ProjectiveSpace(GF(2), 3)
            sage: C = P.curve([(x - y)*(x - z)*(x - w)*(y - z)*(y - w), x*y*z*w*(x+y+z+w)])
            sage: C.projection()
            Traceback (most recent call last):
            ...
            NotImplementedError: this curve contains all points of its ambient space

        ::

            sage: P.<x,y,z,w,u> = ProjectiveSpace(GF(7), 4)
            sage: C = P.curve([x^3 - y*z*u, w^2 - u^2 + 2*x*z, 3*x*w - y^2])
            sage: L.<a,b,c,d> = ProjectiveSpace(GF(7), 3)
            sage: C.projection(PS=L)
            (Scheme morphism:
               From: Projective Curve over Finite Field of size 7 defined by x^3 -
            y*z*u, 2*x*z + w^2 - u^2, -y^2 + 3*x*w
               To:   Projective Space of dimension 3 over Finite Field of size 7
               Defn: Defined on coordinates by sending (x : y : z : w : u) to
                     (x : y : z : w),
             Projective Curve over Finite Field of size 7 defined by b^2 - 3*a*d,
            a^5*b + a*b*c^3*d - 3*b*c^2*d^3, a^6 + a^2*c^3*d - 3*a*c^2*d^3)
            sage: Q.<a,b,c> = ProjectiveSpace(GF(7), 2)
            sage: C.projection(PS=Q)
            Traceback (most recent call last):
            ...
            TypeError: (=Projective Space of dimension 2 over Finite Field of size
            7) must have dimension (=3)
            

        ::

            sage: PP.<x,y,z,w> = ProjectiveSpace(QQ, 3)
            sage: C = PP.curve([x^3 - z^2*y, w^2 - z*x])
            sage: Q = PP([1,0,1,1])
            sage: C.projection(P=Q)
            (Scheme morphism:
               From: Projective Curve over Rational Field defined by x^3 - y*z^2,
            -x*z + w^2
               To:   Projective Space of dimension 2 over Rational Field
               Defn: Defined on coordinates by sending (x : y : z : w) to
                     (y : -x + z : -x + w),
             Projective Plane Curve over Rational Field defined by x0*x1^5 -
            6*x0*x1^4*x2 + 14*x0*x1^3*x2^2 - 16*x0*x1^2*x2^3 + 9*x0*x1*x2^4 -
            2*x0*x2^5 - x2^6)
            sage: LL.<a,b,c> = ProjectiveSpace(QQ, 2)
            sage: Q = PP([0,0,0,1])
            sage: C.projection(PS=LL, P=Q)
            (Scheme morphism:
               From: Projective Curve over Rational Field defined by x^3 - y*z^2,
            -x*z + w^2
               To:   Projective Space of dimension 2 over Rational Field
               Defn: Defined on coordinates by sending (x : y : z : w) to
                     (x : y : z),
             Projective Plane Curve over Rational Field defined by a^3 - b*c^2)
            sage: Q = PP([0,0,1,0])
            sage: C.projection(P=Q)
            Traceback (most recent call last):
            ...
            TypeError: (=(0 : 0 : 1 : 0)) must be a point not on this curve

        ::

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: C = P.curve(y^2 - x^2 + z^2)
            sage: C.projection()
            Traceback (most recent call last):
            ...
            TypeError: this curve is already a plane curve
        """
        PP = self.ambient_space()
        n = PP.dimension_relative()
        if n == 2:
            raise TypeError("this curve is already a plane curve")
        if self.base_ring() not in Fields():
            raise TypeError("this curve must be defined over a field")
        if not PS is None:
            if not is_ProjectiveSpace(PS):
                raise TypeError("(=%s) must be a projective space"%PS)
            if PS.dimension_relative() != n - 1:
                raise TypeError("(=%s) must have dimension (=%s)"%(PS,n - 1))
            if PS.base_ring() != PP.base_ring():
                raise TypeError("(=%s) must be defined over the same base field as this curve"%PS)
        if P is None:
            # find a point not on the curve if not given
            if self.base_ring().characteristic() == 0:
                # when working over a characteristic 0 field, we can construct a point not on the curve.
                # we do this by constructing a point on which at least one nonzero element of the defining ideal of
                # this curve does not vanish
                F = 0
                # find a nonzero element
                for i in range(len(self.defining_polynomials())):
                    if self.defining_polynomials()[i] != 0:
                        F = self.defining_polynomials()[i]
                # find a point on which it doesn't vanish
                l = list(PP.gens())
                for i in range(n + 1):
                    l[i] = 0
                    while(F(l) == 0):
                        l[i] = l[i] + 1
                Q = PP(l) # will be a point not on the curve
            else:
                # if the base ring is a finite field, iterate over all points in the ambient space and check which
                # are on this curve
                Q = None
                for P in PP.rational_points():
                    try:
                        self(P)
                    except TypeError:
                        Q = P
                        break
                if Q is None:
                    raise NotImplementedError("this curve contains all points of its ambient space")
        else:
            # make sure the given point is in the ambient space of the curve, but not on the curve
            Q = None
            try:
                Q = self(P)
            except TypeError:
                pass
            if not Q is None:
                raise TypeError("(=%s) must be a point not on this curve"%P)
            try:
                Q = self.ambient_space()(P)
            except TypeError:
                raise TypeError("(=%s) must be a point in the ambient space of this curve"%P)
        # in order to create the change of coordinates map, need to find a coordinate of Q that is nonzero
        j = 0
        while Q[j] == 0:
            j = j + 1
        # use this Q to project. Apply a change of coordinates to move Q to (0:...:0:1:0:...:0)
        # where 1 is in the jth coordinate
        if PS is None:
            PP2 = ProjectiveSpace(self.base_ring(), n - 1)
        else:
            PP2 = PS
        H = Hom(self, PP2)
        coords = [PP.gens()[i] - Q[i]/Q[j]*PP.gens()[j] for i in range(n + 1)]
        coords.pop(j)
        psi = H(coords)
        # compute image of psi via elimination
        # first construct the image of this curve by the change of coordinates. This can be found by composing the
        # defining polynomials of this curve with the polynomials defining the inverse of the change of coordinates
        invcoords = [Q[i]*PP.gens()[j] + PP.gens()[i] for i in range(n + 1)]
        invcoords[j] = Q[j]*PP.gens()[j]
        I = PP.coordinate_ring().ideal([f(invcoords) for f in self.defining_polynomials()])
        J = I.elimination_ideal(PP.gens()[j])
        K = Hom(PP.coordinate_ring(), PP2.coordinate_ring())
        l = list(PP2.gens())
        l.insert(j, 0)
        phi = K(l)
        G = [phi(f) for f in J.gens()]
        C = PP2.curve(G)
        return tuple([psi, C])

    def plane_projection(self, PP=None):
        r"""
        Return a projection of this curve into a projective plane.

        INPUT:

        - ``PP`` -- (default: None) the projective plane the projected curve will be defined in. This space must
          be defined over the same base field as this curve, and must have dimension two. This space is constructed
          if not specified.

        OUTPUT:

        - a tuple consisting of two elements: a scheme morphism from this curve into a projective plane,
          and the projective curve that is the image of that morphism.

        EXAMPLES::

            sage: P.<x,y,z,w,u,v> = ProjectiveSpace(QQ, 5)
            sage: C = P.curve([x*u - z*v, w - y, w*y - x^2, y^3*u*2*z - w^4*w])
            sage: L.<a,b,c> = ProjectiveSpace(QQ, 2)
            sage: proj1 = C.plane_projection(PP=L)
            sage: proj1
            (Scheme morphism:
               From: Projective Curve over Rational Field defined by x*u - z*v, -y +
            w, -x^2 + y*w, -w^5 + 2*y^3*z*u
               To:   Projective Space of dimension 2 over Rational Field
               Defn: Defined on coordinates by sending (x : y : z : w : u : v) to
                     (x : -z + u : -z + v),
             Projective Plane Curve over Rational Field defined by a^8 + 6*a^7*b +
            4*a^5*b^3 - 4*a^7*c - 2*a^6*b*c - 4*a^5*b^2*c + 2*a^6*c^2)
            sage: proj1[1].ambient_space() is L
            True
            sage: proj2 = C.projection()
            sage: proj2[1].ambient_space() is L
            False

        ::

            sage: P.<x,y,z,w,u> = ProjectiveSpace(GF(7), 4)
            sage: C = P.curve([x^2 - 6*y^2, w*z*u - y^3 + 4*y^2*z, u^2 - x^2])
            sage: C.plane_projection()
            (Scheme morphism:
               From: Projective Curve over Finite Field of size 7 defined by x^2 +
            y^2, -y^3 - 3*y^2*z + z*w*u, -x^2 + u^2
               To:   Projective Space of dimension 2 over Finite Field of size 7
               Defn: Defined on coordinates by sending (x : y : z : w : u) to
                     (y : z : -x + w),
             Projective Plane Curve over Finite Field of size 7 defined by x0^10 -
            2*x0^9*x1 + 3*x0^8*x1^2 - 2*x0^7*x1^3 + x0^6*x1^4 + 2*x0^6*x1^2*x2^2 -
            2*x0^5*x1^3*x2^2 - x0^4*x1^4*x2^2 + x0^2*x1^4*x2^4)

        ::

            sage: P.<x,y,z> = ProjectiveSpace(GF(17), 2)
            sage: C = P.curve(x^2 - y*z - z^2)
            sage: C.plane_projection()
            Traceback (most recent call last):
            ...
            TypeError: this curve is already a plane curve
        """
        PS = self.ambient_space()
        n = PS.dimension_relative()
        if n == 2:
            raise TypeError("this curve is already a plane curve")
        C = self
        H = Hom(PS, PS)
        phi = H([PS.gens()[i] for i in range(n + 1)])
        for i in range(n - 2):
            if i == n - 3:
                L = C.projection(PS=PP)
            else:
                L = C.projection()
            C = L[1]
            # compose the scheme morphisms that are created
            K = Hom(phi.codomain().coordinate_ring(), PS.coordinate_ring())
            psi = K(phi.defining_polynomials())
            H = Hom(self, L[1].ambient_space())
            phi = H([psi(L[0].defining_polynomials()[i]) for i in range(len(L[0].defining_polynomials()))])
        return tuple([phi, C])

    def arithmetic_genus(self):
        r"""
        Return the arithmetic genus of this projective curve.

        This is the arithmetic genus `g_a(C)` as defined in [Hartshorne]_. If `P` is the
        Hilbert polynomial of the defining ideal of this curve, then the arithmetic genus
        of this curve is `1 - P(0)`. This curve must be irreducible.

        OUTPUT: Integer.

        EXAMPLES::

            sage: P.<x,y,z,w> = ProjectiveSpace(QQ, 3)
            sage: C = P.curve([w*z - x^2, w^2 + y^2 + z^2])
            sage: C.arithmetic_genus()
            1

        ::

            sage: P.<x,y,z,w,t> = ProjectiveSpace(GF(7), 4)
            sage: C = P.curve([t^3 - x*y*w, x^3 + y^3 + z^3, z - w])
            sage: C.arithmetic_genus()
            10
        """
        if not self.is_irreducible():
            raise TypeError("this curve must be irreducible")
        return 1 - self.defining_ideal().hilbert_polynomial()(0)

    def is_complete_intersection(self):
        r"""
        Return whether this projective curve is or is not a complete intersection.

        OUTPUT: Boolean.

        EXAMPLES::

            sage: P.<x,y,z,w> = ProjectiveSpace(QQ, 3)
            sage: C = Curve([x*y - z*w, x^2 - y*w, y^2*w - x*z*w], P)
            sage: C.is_complete_intersection()
            False

        ::

            sage: P.<x,y,z,w> = ProjectiveSpace(QQ, 3)
            sage: C = Curve([y*w - x^2, z*w^2 - x^3], P)
            sage: C.is_complete_intersection()
            True

            sage: P.<x,y,z,w> = ProjectiveSpace(QQ, 3)
            sage: C = Curve([z^2 - y*w, y*z - x*w, y^2 - x*z], P)
            sage: C.is_complete_intersection()
            False
        """
        singular.lib("sing.lib")
        I = singular.simplify(self.defining_ideal(), 10)
        L = singular.is_ci(I).sage()
        return len(self.ambient_space().gens()) - len(I.sage().gens()) == L[-1]

class ProjectivePlaneCurve(ProjectiveCurve):

    _point = point.ProjectivePlaneCurvePoint_field

    def __init__(self, A, f):
        r"""
        Initialization function.

        EXAMPLES::

            sage: P.<x,y,z> = ProjectiveSpace(QQbar, 2)
            sage: C = Curve([y*z - x^2 - QQbar.gen()*z^2], P); C
            Projective Plane Curve over Algebraic Field defined by
            -x^2 + y*z + (-I)*z^2

        ::

            sage: P.<x,y,z> = ProjectiveSpace(GF(5^2, 'v'), 2)
            sage: C = Curve([y^2*z - x*z^2 - z^3], P); C
            Projective Plane Curve over Finite Field in v of size 5^2 defined by y^2*z - x*z^2 - z^3
        """
        if not (is_ProjectiveSpace(A) and A.dimension != 2):
            raise TypeError("Argument A (= %s) must be a projective plane."%A)
        Curve_generic.__init__(self, A, [f])

    def _repr_type(self):
        r"""
        Return a string representation of the type of this curve.

        EXAMPLES::

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: C = Curve([y*z^3 - 5/7*x^4 + 4*x^3*z - 9*z^4], P)
            sage: C._repr_type()
            'Projective Plane'
        """
        return "Projective Plane"

    def arithmetic_genus(self):
        r"""
        Return the arithmetic genus of this projective curve.

        This is the arithmetic genus `g_a(C)` as defined in [Hartshorne]_. For a projective
        plane curve of degree `d`, this is simply `(d-1)(d-2)/2`. It need *not* equal
        the geometric genus (the genus of the normalization of the curve). This curve must be
        irreducible.

        OUTPUT: Integer.

        EXAMPLE::

            sage: x,y,z = PolynomialRing(GF(5), 3, 'xyz').gens()
            sage: C = Curve(y^2*z^7 - x^9 - x*z^8); C
            Projective Plane Curve over Finite Field of size 5 defined by -x^9 + y^2*z^7 - x*z^8
            sage: C.arithmetic_genus()
            28
            sage: C.genus()
            4

        ::

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: C = Curve([y^3*x - x^2*y*z - 7*z^4])
            sage: C.arithmetic_genus()
            3

        REFERENCES:

        ..  [Hartshorne] \R. Hartshorne. Algebraic Geometry. Springer-Verlag, New York, 1977.
        """
        if not self.is_irreducible():
            raise TypeError("this curve must be irreducible")
        d = self.defining_polynomial().total_degree()
        return int((d-1)*(d-2)/2)

    def divisor_of_function(self, r):
        """
        Return the divisor of a function on a curve.

        INPUT: r is a rational function on X

        OUTPUT:


        -  ``list`` - The divisor of r represented as a list of
           coefficients and points. (TODO: This will change to a more
           structural output in the future.)


        EXAMPLES::

            sage: FF = FiniteField(5)
            sage: P2 = ProjectiveSpace(2, FF, names = ['x','y','z'])
            sage: R = P2.coordinate_ring()
            sage: x, y, z = R.gens()
            sage: f = y^2*z^7 - x^9 - x*z^8
            sage: C = Curve(f)
            sage: K = FractionField(R)
            sage: r = 1/x
            sage: C.divisor_of_function(r)     # todo: not implemented  !!!!
            [[-1, (0, 0, 1)]]
            sage: r = 1/x^3
            sage: C.divisor_of_function(r)     # todo: not implemented  !!!!
            [[-3, (0, 0, 1)]]
        """
        F = self.base_ring()
        f = self.defining_polynomial()
        x, y, z = f.parent().gens()
        pnts = self.rational_points()
        divf = []
        for P in pnts:
            if P[2] != F(0):
                # What is the '5' in this line and the 'r()' in the next???
                lcs = self.local_coordinates(P,5)
                ldg = degree_lowest_rational_function(r(lcs[0],lcs[1]),z)
                if ldg[0] != 0:
                    divf.append([ldg[0],P])
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

            sage: FF = FiniteField(5)
            sage: P2 = ProjectiveSpace(2, FF, names = ['x','y','z'])
            sage: x, y, z = P2.coordinate_ring().gens()
            sage: C = Curve(y^2*z^7-x^9-x*z^8)
            sage: pt = C([2,3,1])
            sage: C.local_coordinates(pt,9)     # todo: not implemented  !!!!
                  [2 + t, 3 + 3*t^2 + t^3 + 3*t^4 + 3*t^6 + 3*t^7 + t^8 + 2*t^9 + 3*t^11 + 3*t^12]
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
        yt = y0*t**0 + add([vars0[i]*t**(i-2) for i in range(3,2*n+2)])
        xt = x0+t
        ft = f(xt,yt)
        S = singular
        S.eval('ring s = '+str(p)+','+str(R0.gens())+',lp;')
        S.eval('poly f = '+str(ft))
        cmd = 'matrix c = coeffs ('+str(ft)+',t)'
        S.eval(cmd)
        N = int(S.eval('size(c)'))
        b = ["c["+str(i)+",1]," for i in range(2, N//2 - 4)]
        b = ''.join(b)
        b = b[:len(b)-1] #to cut off the trailing comma
        cmd = 'ideal I = '+b
        S.eval(cmd)
        c = S.eval('slimgb(I)')
        d = c.split("=")
        d = d[1:]
        d[len(d)-1] += "\n"
        e = [x[:x.index("\n")] for x in d]
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
        k = len(vals)
        v = [x0+t,y0+add([vals[i][2]*t**(i+1) for i in range(k)])]
        return v

    def plot(self, *args, **kwds):
        """
        Plot the real points of an affine patch of this projective
        plane curve.


        INPUT:


        -  ``self`` - an affine plane curve

        -  ``patch`` - (optional) the affine patch to be plotted; if not
           specified, the patch corresponding to the last projective
           coordinate being nonzero

        -  ``*args`` - optional tuples (variable, minimum, maximum) for
           plotting dimensions

        -  ``**kwds`` - optional keyword arguments passed on to
           ``implicit_plot``


        EXAMPLES:

        A cuspidal curve::

            sage: R.<x, y, z> = QQ[]
            sage: C = Curve(x^3 - y^2*z)
            sage: C.plot()
            Graphics object consisting of 1 graphics primitive

        The other affine patches of the same curve::

            sage: C.plot(patch=0)
            Graphics object consisting of 1 graphics primitive
            sage: C.plot(patch=1)
            Graphics object consisting of 1 graphics primitive

        An elliptic curve::

            sage: E = EllipticCurve('101a')
            sage: C = Curve(E)
            sage: C.plot()
            Graphics object consisting of 1 graphics primitive
            sage: C.plot(patch=0)
            Graphics object consisting of 1 graphics primitive
            sage: C.plot(patch=1)
            Graphics object consisting of 1 graphics primitive

        A hyperelliptic curve::

            sage: P.<x> = QQ[]
            sage: f = 4*x^5 - 30*x^3 + 45*x - 22
            sage: C = HyperellipticCurve(f)
            sage: C.plot()
            Graphics object consisting of 1 graphics primitive
            sage: C.plot(patch=0)
            Graphics object consisting of 1 graphics primitive
            sage: C.plot(patch=1)
            Graphics object consisting of 1 graphics primitive
        """
        # if user hasn't specified a favourite affine patch, take the
        # one avoiding "infinity", i.e. the one corresponding to the
        # last projective coordinate being nonzero
        patch = kwds.pop('patch', self.ngens() - 1)
        from .constructor import Curve
        C = Curve(self.affine_patch(patch))
        return C.plot(*args, **kwds)

    def is_singular(self, P=None):
        r"""
        Return whether this curve is singular or not, or if a point ``P`` is provided,
        whether ``P`` is a singular point of this curve.

        INPUT:

        - ``P`` -- (default: None) a point on this curve.

        OUTPUT:

        - Boolean. If no point ``P`` is provided, returns True of False depending on whether
          this curve is singular or not. If a point ``P`` is provided, returns True or False
          depending on whether ``P`` is or is not a singular point of this curve.

        EXAMPLES:

        Over `\QQ`::

            sage: F = QQ
            sage: P2.<X,Y,Z> = ProjectiveSpace(F,2)
            sage: C = Curve(X^3-Y^2*Z)
            sage: C.is_singular()
            True

        Over a finite field::

            sage: F = GF(19)
            sage: P2.<X,Y,Z> = ProjectiveSpace(F,2)
            sage: C = Curve(X^3+Y^3+Z^3)
            sage: C.is_singular()
            False
            sage: D = Curve(X^4-X*Z^3)
            sage: D.is_singular()
            True
            sage: E = Curve(X^5+19*Y^5+Z^5)
            sage: E.is_singular()
            True
            sage: E = Curve(X^5+9*Y^5+Z^5)
            sage: E.is_singular()
            False

        Over `\CC`::

            sage: F = CC
            sage: P2.<X,Y,Z> = ProjectiveSpace(F,2)
            sage: C = Curve(X)
            sage: C.is_singular()
            False
            sage: D = Curve(Y^2*Z-X^3)
            sage: D.is_singular()
            True
            sage: E = Curve(Y^2*Z-X^3+Z^3)
            sage: E.is_singular()
            False

        Showing that :trac:`12187` is fixed::

            sage: F.<X,Y,Z> = GF(2)[]
            sage: G = Curve(X^2+Y*Z)
            sage: G.is_singular()
            False

        ::

            sage: P.<x,y,z> = ProjectiveSpace(CC, 2)
            sage: C = Curve([y^4 - x^3*z], P)
            sage: Q = P([0,0,1])
            sage: C.is_singular()
            True
        """
        if P is None:
            poly = self.defining_polynomial()
            return poly.parent().ideal(poly.gradient()+[poly]).dimension() > 0
        else:
            return not self.is_smooth(P)

    def degree(self):
        r"""
        Return the degree of this projective curve.

        For a plane curve, this is just the degree of its defining polynomial.

        OUTPUT: integer.

        EXAMPLES::

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: C = P.curve([y^7 - x^2*z^5 + 7*z^7])
            sage: C.degree()
            7
        """
        return self.defining_polynomial().degree()

    def tangents(self, P, factor=True):
            r"""
            Return the tangents of this projective plane curve at the point ``P``.

            These are found by homogenizing the tangents of an affine patch of this curve
            containing ``P``. The point ``P`` must be a point on this curve.

            INPUT:

            - ``P`` -- a point on this curve.

            - ``factor`` -- (default: True) whether to attempt computing the polynomials of the individual tangent
              lines over the base field of this curve, or to just return the polynomial corresponding to the union
              of the tangent lines (which requires fewer computations).

            OUTPUT:

            - a list of polynomials in the coordinate ring of the ambient space of this curve.

            EXAMPLES::

                sage: set_verbose(-1)
                sage: P.<x,y,z> = ProjectiveSpace(QQbar, 2)
                sage: C = Curve([x^3*y + 2*x^2*y^2 + x*y^3 + x^3*z + 7*x^2*y*z + 14*x*y^2*z + 9*y^3*z], P)
                sage: Q = P([0,0,1])
                sage: C.tangents(Q)
                [x + 4.147899035704788?*y, x + (1.426050482147607? + 0.3689894074818041?*I)*y,
                x + (1.426050482147607? - 0.3689894074818041?*I)*y]
                sage: C.tangents(Q, factor=False)
                [6*x^3 + 42*x^2*y + 84*x*y^2 + 54*y^3]

            ::

                sage: P.<x,y,z> = ProjectiveSpace(QQ,2)
                sage: C = P.curve([x^2*y^3*z^4 - y^6*z^3 - 4*x^2*y^4*z^3 - 4*x^4*y^2*z^3 + 3*y^7*z^2 +\
                10*x^2*y^5*z^2 + 9*x^4*y^3*z^2 + 5*x^6*y*z^2 - 3*y^8*z - 9*x^2*y^6*z - 11*x^4*y^4*z -\
                7*x^6*y^2*z - 2*x^8*z + y^9 + 2*x^2*y^7 + 3*x^4*y^5 + 4*x^6*y^3 + 2*x^8*y])
                sage: Q = P([0,1,1])
                sage: C.tangents(Q)
                [-y + z, 3*x^2 - y^2 + 2*y*z - z^2]

            ::

                sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
                sage: C = P.curve([z^3*x + y^4 - x^2*z^2])
                sage: Q = P([1,1,1])
                sage: C.tangents(Q)
                Traceback (most recent call last):
                ...
                TypeError: (=(1 : 1 : 1)) is not a point on (=Projective Plane Curve
                over Rational Field defined by y^4 - x^2*z^2 + x*z^3)
            """
            PP = self.ambient_space()
            # Check whether P is a point on this curve
            try:
                P = self(P)
            except TypeError:
                raise TypeError("(=%s) is not a point on (=%s)"%(P,self))

            # Find an affine chart of the ambient space of self that contains P
            i = 0
            while(P[i] == 0):
                i = i + 1
            C = self.affine_patch(i)
            L = C.tangents(C(P.dehomogenize(i)), factor)
            R = PP.coordinate_ring()
            H = Hom(C.ambient_space().coordinate_ring(), R)
            G = list(R.gens())
            x = G.pop(i)
            phi = H(G)
            return [phi(g).homogenize(x) for g in L]

    def is_ordinary_singularity(self, P):
        r"""
        Return whether the singular point ``P`` of this projective plane curve is an ordinary singularity.

        The point ``P`` is an ordinary singularity of this curve if it is a singular point, and
        if the tangents of this curve at ``P`` are distinct.

        INPUT:

        - ``P`` -- a point on this curve.

        OUTPUT:

        - Boolean. True or False depending on whether ``P`` is or is not an ordinary singularity of this
          curve, respectively. An error is raised if ``P`` is not a singular point of this curve.

        EXAMPLES::

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: C = Curve([y^2*z^3 - x^5], P)
            sage: Q = P([0,0,1])
            sage: C.is_ordinary_singularity(Q)
            False

        ::

            sage: R.<a> = QQ[]
            sage: K.<b> = NumberField(a^2 - 3)
            sage: P.<x,y,z> = ProjectiveSpace(K, 2)
            sage: C = P.curve([x^2*y^3*z^4 - y^6*z^3 - 4*x^2*y^4*z^3 - 4*x^4*y^2*z^3 + 3*y^7*z^2 + 10*x^2*y^5*z^2\
            + 9*x^4*y^3*z^2 + 5*x^6*y*z^2 - 3*y^8*z - 9*x^2*y^6*z - 11*x^4*y^4*z - 7*x^6*y^2*z - 2*x^8*z + y^9 +\
            2*x^2*y^7 + 3*x^4*y^5 + 4*x^6*y^3 + 2*x^8*y])
            sage: Q = P([0,1,1])
            sage: C.is_ordinary_singularity(Q)
            True

        ::

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: C = P.curve([z^5 - y^5 + x^5 + x*y^2*z^2])
            sage: Q = P([0,1,1])
            sage: C.is_ordinary_singularity(Q)
            Traceback (most recent call last):
            ...
            TypeError: (=(0 : 1 : 1)) is not a singular point of (=Projective Plane
            Curve over Rational Field defined by x^5 - y^5 + x*y^2*z^2 + z^5)
        """
        r = self.multiplicity(P)
        if r < 2:
            raise TypeError("(=%s) is not a singular point of (=%s)"%(P,self))

        # Find an affine chart of the ambient space of self that contains P
        i = 0
        while(P[i] == 0):
            i = i + 1
        C = self.affine_patch(i)
        return C.is_ordinary_singularity(C(P.dehomogenize(i)))

    def quadratic_transform(self):
        r"""
        Return a birational map from this curve to the proper transform of this curve with respect to the standard
        Cremona transformation.

        The standard Cremona transformation is the birational automorphism of `\mathbb{P}^{2}` defined
        `(x : y : z)\mapsto (yz : xz : xy)`.

        OUTPUT:

        - a scheme morphism representing the restriction of the standard Cremona transformation from this curve
          to the proper transform.

        EXAMPLES::

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: C = Curve(x^3*y - z^4 - z^2*x^2, P)
            sage: C.quadratic_transform()
            Scheme morphism:
              From: Projective Plane Curve over Rational Field defined by x^3*y -
            x^2*z^2 - z^4
              To:   Projective Plane Curve over Rational Field defined by -x^3*y -
            x*y*z^2 + z^4
              Defn: Defined on coordinates by sending (x : y : z) to
                    (y*z : x*z : x*y)

        ::

            sage: P.<x,y,z> = ProjectiveSpace(GF(17), 2)
            sage: C = P.curve([y^7*z^2 - 16*x^9 + x*y*z^7 + 2*z^9])
            sage: C.quadratic_transform()
            Scheme morphism:
              From: Projective Plane Curve over Finite Field of size 17 defined by
            x^9 + y^7*z^2 + x*y*z^7 + 2*z^9
              To:   Projective Plane Curve over Finite Field of size 17 defined by
            2*x^9*y^7 + x^8*y^6*z^2 + x^9*z^7 + y^7*z^9
              Defn: Defined on coordinates by sending (x : y : z) to
                    (y*z : x*z : x*y)
        """
        PP = self.ambient_space()
        R = PP.coordinate_ring()
        L = R.gens()
        coords = [L[1]*L[2], L[0]*L[2], L[0]*L[1]]
        G = self.defining_polynomial()(coords)
        # remove the component of the curve corresponding to the exceptional divisor
        degs = [G.degree()]*len(L)
        for F in G.monomials():
            for i in range(len(L)):
                if F.degree(L[i]) < degs[i]:
                    degs[i] = F.degree(L[i])
        T = []
        for item in G.dict().items():
            tup = tuple([item[0][i] - degs[i] for i in range(len(L))])
            T.append(tuple([tup, item[1]]))
        G = R(dict(T))
        H = Hom(self, PP.curve(G))
        phi = H(coords)
        return phi

    def excellent_position(self, Q):
        r"""
        Return a transformation of this curve into one in excellent position with respect to the point ``Q``.

        Here excellent position is defined as in [Fulton89]_. A curve `C` of degree `d` containing the point
        `(0 : 0 : 1)` with multiplicity `r` is said to be in excellent position if none of the coordinate lines
        are tangent to `C` at any of the fundamental points `(1 : 0 : 0)`, `(0 : 1 : 0)`, and `(0 : 0 : 1)`, and
        if the two coordinate lines containing `(0 : 0 : 1)` intersect `C` transversally in `d - r` distinct
        non-fundamental points, and if the other coordinate line intersects `C` transversally at `d` distinct,
        non-fundamental points.

        INPUT:

        - ``Q`` -- a point on this curve.

        OUTPUT:

        - a scheme morphism from this curve to a curve in excellent position that is a restriction of a change
          of coordinates map of the projective plane.

        EXAMPLES::

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: C = Curve([x*y - z^2], P)
            sage: Q = P([1,1,1])
            sage: C.excellent_position(Q)
            Scheme morphism:
              From: Projective Plane Curve over Rational Field defined by x*y - z^2
              To:   Projective Plane Curve over Rational Field defined by -x^2 -
            3*x*y - 4*y^2 - x*z - 3*y*z
              Defn: Defined on coordinates by sending (x : y : z) to
                    (-x + 1/2*y + 1/2*z : -1/2*y + 1/2*z : x + 1/2*y - 1/2*z)

        ::

            sage: R.<a> = QQ[]
            sage: K.<b> = NumberField(a^2 - 3)
            sage: P.<x,y,z> = ProjectiveSpace(K, 2)
            sage: C = P.curve([z^2*y^3*x^4 - y^6*x^3 - 4*z^2*y^4*x^3 - 4*z^4*y^2*x^3 + 3*y^7*x^2 + 10*z^2*y^5*x^2\
            + 9*z^4*y^3*x^2 + 5*z^6*y*x^2 - 3*y^8*x - 9*z^2*y^6*x - 11*z^4*y^4*x - 7*z^6*y^2*x - 2*z^8*x + y^9 +\
            2*z^2*y^7 + 3*z^4*y^5 + 4*z^6*y^3 + 2*z^8*y])
            sage: Q = P([1,0,0])
            sage: C.excellent_position(Q)
            Scheme morphism:
              From: Projective Plane Curve over Number Field in b with defining
            polynomial a^2 - 3 defined by -x^3*y^6 + 3*x^2*y^7 - 3*x*y^8 + y^9 +
            x^4*y^3*z^2 - 4*x^3*y^4*z^2 + 10*x^2*y^5*z^2 - 9*x*y^6*z^2 + 2*y^7*z^2 -
            4*x^3*y^2*z^4 + 9*x^2*y^3*z^4 - 11*x*y^4*z^4 + 3*y^5*z^4 + 5*x^2*y*z^6 -
            7*x*y^2*z^6 + 4*y^3*z^6 - 2*x*z^8 + 2*y*z^8
              To:   Projective Plane Curve over Number Field in b with defining
            polynomial a^2 - 3 defined by 900*x^9 - 7410*x^8*y + 29282*x^7*y^2 -
            69710*x^6*y^3 + 110818*x^5*y^4 - 123178*x^4*y^5 + 96550*x^3*y^6 -
            52570*x^2*y^7 + 18194*x*y^8 - 3388*y^9 - 1550*x^8*z + 9892*x^7*y*z -
            30756*x^6*y^2*z + 58692*x^5*y^3*z - 75600*x^4*y^4*z + 67916*x^3*y^5*z -
            42364*x^2*y^6*z + 16844*x*y^7*z - 3586*y^8*z + 786*x^7*z^2 -
            3958*x^6*y*z^2 + 9746*x^5*y^2*z^2 - 14694*x^4*y^3*z^2 +
            15174*x^3*y^4*z^2 - 10802*x^2*y^5*z^2 + 5014*x*y^6*z^2 - 1266*y^7*z^2 -
            144*x^6*z^3 + 512*x^5*y*z^3 - 912*x^4*y^2*z^3 + 1024*x^3*y^3*z^3 -
            816*x^2*y^4*z^3 + 512*x*y^5*z^3 - 176*y^6*z^3 + 8*x^5*z^4 - 8*x^4*y*z^4
            - 16*x^3*y^2*z^4 + 16*x^2*y^3*z^4 + 8*x*y^4*z^4 - 8*y^5*z^4
              Defn: Defined on coordinates by sending (x : y : z) to
                    (1/4*y + 1/2*z : -1/4*y + 1/2*z : x + 1/4*y - 1/2*z)

        ::

            sage: set_verbose(-1)
            sage: a = QQbar(sqrt(2))
            sage: P.<x,y,z> = ProjectiveSpace(QQbar, 2)
            sage: C = Curve([(-1/4*a)*x^3 + (-3/4*a)*x^2*y + (-3/4*a)*x*y^2 + (-1/4*a)*y^3 - 2*x*y*z], P)
            sage: Q = P([0,0,1])
            sage: C.excellent_position(Q)
            Scheme morphism:
              From: Projective Plane Curve over Algebraic Field defined by
            (-0.3535533905932738?)*x^3 + (-1.060660171779822?)*x^2*y +
            (-1.060660171779822?)*x*y^2 + (-0.3535533905932738?)*y^3 + (-2)*x*y*z
              To:   Projective Plane Curve over Algebraic Field defined by
            (-2.828427124746190?)*x^3 + (-2)*x^2*y + 2*y^3 + (-2)*x^2*z + 2*y^2*z
              Defn: Defined on coordinates by sending (x : y : z) to
                    (1/2*x + 1/2*y : (-1/2)*x + 1/2*y : 1/2*x + (-1/2)*y + z)

        REFERENCES:

        ..  [Fulton89] \W. Fulton. Algebraic curves: an introduction to algebraic geometry. Addison-Wesley,
            Redwood City CA (1989).
        """
        PP = self.ambient_space()
        # check that Q is on this curve
        try:
            Q = self(Q)
        except TypeError:
            raise TypeError("(=%s) must be a point on this curve"%Q)
        r = self.multiplicity(Q)
        d = self.degree()
        # first move Q to (0 : 0 : 1), (1 : 0 : 0), or (0 : 1 : 0)
        # this makes it easier to construct the main transformation
        i = 0
        while Q[i] == 0:
            i = i + 1
        coords = [PP.gens()[j] + Q[j]/Q[i]*PP.gens()[i] for j in range(3)]
        coords[i] = PP.gens()[i]
        accoords = [PP.gens()[j] - Q[j]/Q[i]*PP.gens()[i] for j in range(3)] # coords used in map construction
        accoords[i] = PP.gens()[i]
        baseC = PP.curve(self.defining_polynomial()(coords))
        P = [0]*3
        P[i] = 1
        P = PP(P)
        l = [0,1,2]
        l.pop(i)
        # choose points forming a triangle with one vertex at P to map to the coordinate triangle
        good = False
        a = 0
        while not good:
            a = a + 1
            # find points to map to (1 : 0 : 0) and (0 : 1 : 0), not on the curve
            Px = [0]*3
            Px[l[0]] = a
            Px[l[1]] = 1
            Py = [0]*3
            Py[l[0]] = -a
            Py[l[1]] = 1
            Py[i] = 1
            try:
                Px = baseC(Px)
                Py = baseC(Py)
                continue
            except TypeError:
                pass
            # by construction, P, Px, Py are linearly independent so the following matrix is invertible
            M = matrix([[Px[j], Py[j], P[j]] for j in range(3)])
            # M defines a change of coordinates sending (1 : 0 : 0) to Py, (0 : 1 : 0) to Px, (0 : 0 : 1) to P; the
            # inverse of the transformation we want, used to create the new defining polynomial
            coords = [sum([M.row(j)[k]*PP.gens()[k] for k in range(3)]) for j in range(3)]
            C = PP.curve(baseC.defining_polynomial()(coords))
            # check tangents at (0 : 0 : 1)
            T = C.tangents(PP([0,0,1]), factor=False)[0]
            if all([e[0] > 0 for e in T.exponents()]) or all([e[1] > 0 for e in T.exponents()]):
                continue
            # check that the other intersections of C with the exceptional lines are correct
            need_continue = False
            for j in range(3):
                poly = C.defining_polynomial().subs({PP.gens()[j]: 0})
                # this is a homogeneous polynomial in the other two variables
                # and so should factor completely into homogeneous linear factors
                # each corresponding to an intersection point where the jth coord is 0.
                # check if there are enough roots, up to multiplicity (that is, that PP.gens()[j]
                # doesn't divide the defining polynomial of C)
                if poly.degree() != d:
                    need_continue = True
                    break
                # if j != 2, then there should be d - r multiplicity 1 roots,
                # besides the root corresponding to (0 : 0 : 1)
                # if j == 2, then all roots should have multiplicity 1
                npoly = poly
                if j != 2:
                    # since (0 : 0 : 1) has multiplicity r, divide out by the highest
                    # shared power of the corresponding variable before doing the resultant computations
                    if j == 0:
                        div_pow = min([e[1] for e in npoly.exponents()])
                        npoly = PP.coordinate_ring()(dict([((v[0],v[1] - div_pow,v[2]),g) for (v,g) in\
                                                         npoly.dict().items()]))
                    else:
                        div_pow = min([e[0] for e in npoly.exponents()])
                        npoly = PP.coordinate_ring()(dict([((v[0] - div_pow,v[1],v[2]),g) for (v,g) in\
                                                         npoly.dict().items()]))
                    # check the degree again
                    if npoly.degree() != d - r:
                        need_continue = True
                        break
                    # check that npoly isn't a constant now
                    if npoly.degree() > 0:
                        t = 0
                        while npoly.degree(PP.gens()[t]) == 0:
                            t = t + 1
                        if npoly.resultant(npoly.derivative(PP.gens()[t]), PP.gens()[t]) == 0:
                            need_continue = True
                            break
                else:
                    t = 0
                    while npoly.degree(PP.gens()[t]) == 0:
                        t = t + 1
                    if poly.resultant(poly.derivative(PP.gens()[t]), PP.gens()[t]) == 0:
                        need_continue = True
                        break
                # check that intersections with the line PP.gens()[j] are transverse.
                # at a simple point P of the curve, the tangent at that point is
                # given by F_x(P)*x + F_y(P)*y + F_z(P)*z where F is the defining polynomial
                # of the curve
                tmp_l = [0,1,2]
                tmp_l.pop(j)
                poly1 = npoly.derivative(PP.gens()[tmp_l[0]])
                poly2 = npoly.derivative(PP.gens()[tmp_l[1]])
                if poly1.degree() > 0 or poly2.degree() > 0:
                    t = 0
                    while poly1.degree(PP.gens()[t]) == 0 and poly2.degree(PP.gens()[t]) == 0:
                        t = t + 1
                    # maybe a stricter check than necessary
                    if poly1.resultant(poly2, PP.gens()[t]) == 0:
                        need_continue = True
                        break
            if need_continue:
                continue
            good = True
            # coords for map
            M = M.inverse()
            accoords2 = [sum([M.row(j)[k]*PP.gens()[k] for k in range(3)]) for j in range(3)]
            H = Hom(self, C)
            phi = H([f(accoords) for f in accoords2])
        return phi

    def ordinary_model(self):
        r"""
        Return a birational map from this curve to a plane curve with only ordinary singularities.

        Currently only implemented over number fields. If not all of the coordinates of the non-ordinary
        singularities of this curve are contained in its base field, then the domain and codomain of the
        map returned will be defined over an extension. This curve must be irreducible.

        OUTPUT:

        - a scheme morphism from this curve to a curve with only ordinary singularities that defines a
          birational map between the two curves.

        EXAMPLES::

            sage: set_verbose(-1)
            sage: K = QuadraticField(3)
            sage: P.<x,y,z> = ProjectiveSpace(K, 2)
            sage: C = Curve([x^5 - K.0*y*z^4], P)
            sage: C.ordinary_model()
            Scheme morphism:
              From: Projective Plane Curve over Number Field in a with defining
            polynomial x^2 - 3 defined by x^5 + (-a)*y*z^4
              To:   Projective Plane Curve over Number Field in a with defining
            polynomial x^2 - 3 defined by (-a)*x^5*y + (-4*a)*x^4*y^2 +
            (-6*a)*x^3*y^3 + (-4*a)*x^2*y^4 + (-a)*x*y^5 + (-a - 1)*x^5*z + (-4*a +
            5)*x^4*y*z + (-6*a - 10)*x^3*y^2*z + (-4*a + 10)*x^2*y^3*z + (-a -
            5)*x*y^4*z + y^5*z
              Defn: Defined on coordinates by sending (x : y : z) to
                    (-1/4*x^2 - 1/2*x*y + 1/2*x*z + 1/2*y*z - 1/4*z^2 : 1/4*x^2 +
            1/2*x*y + 1/2*y*z - 1/4*z^2 : -1/4*x^2 + 1/4*z^2)

        ::

            sage: set_verbose(-1)
            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: C = Curve([y^2*z^2 - x^4 - x^3*z], P)
            sage: D = C.ordinary_model(); D # long time (2 seconds)
            Scheme morphism:
              From: Projective Plane Curve over Rational Field defined by -x^4 -
            x^3*z + y^2*z^2
              To:   Projective Plane Curve over Rational Field defined by 4*x^6*y^3
            - 24*x^5*y^4 + 36*x^4*y^5 + 8*x^6*y^2*z - 40*x^5*y^3*z + 24*x^4*y^4*z +
            72*x^3*y^5*z - 4*x^6*y*z^2 + 8*x^5*y^2*z^2 - 56*x^4*y^3*z^2 +
            104*x^3*y^4*z^2 + 44*x^2*y^5*z^2 + 8*x^6*z^3 - 16*x^5*y*z^3 -
            24*x^4*y^2*z^3 + 40*x^3*y^3*z^3 + 48*x^2*y^4*z^3 + 8*x*y^5*z^3 -
            8*x^5*z^4 + 36*x^4*y*z^4 - 56*x^3*y^2*z^4 + 20*x^2*y^3*z^4 +
            40*x*y^4*z^4 - 16*y^5*z^4
              Defn: Defined on coordinates by sending (x : y : z) to
                    (-3/64*x^4 + 9/64*x^2*y^2 - 3/32*x*y^3 - 1/16*x^3*z -
            1/8*x^2*y*z + 1/4*x*y^2*z - 1/16*y^3*z - 1/8*x*y*z^2 + 1/16*y^2*z^2 :
            -1/64*x^4 + 3/64*x^2*y^2 - 1/32*x*y^3 + 1/16*x*y^2*z - 1/16*y^3*z +
            1/16*y^2*z^2 : 3/64*x^4 - 3/32*x^3*y + 3/64*x^2*y^2 + 1/16*x^3*z -
            3/16*x^2*y*z + 1/8*x*y^2*z - 1/8*x*y*z^2 + 1/16*y^2*z^2)
            sage: all([D.codomain().is_ordinary_singularity(Q) for Q in D.codomain().singular_points()]) # long time
            True

        ::

            sage: set_verbose(-1)
            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: C = Curve([(x^2 + y^2 - y*z - 2*z^2)*(y*z - x^2 + 2*z^2)*z + y^5], P)
            sage: C.ordinary_model() # long time (5 seconds)
            Scheme morphism:
              From: Projective Plane Curve over Number Field in a with defining
            polynomial y^2 - 2 defined by y^5 - x^4*z - x^2*y^2*z + 2*x^2*y*z^2 +
            y^3*z^2 + 4*x^2*z^3 + y^2*z^3 - 4*y*z^4 - 4*z^5
              To:   Projective Plane Curve over Number Field in a with defining
            polynomial y^2 - 2 defined by (-29*a + 1)*x^8*y^6 + (10*a + 158)*x^7*y^7
            + (-109*a - 31)*x^6*y^8 + (-80*a - 198)*x^8*y^5*z + (531*a +
            272)*x^7*y^6*z + (170*a - 718)*x^6*y^7*z + (19*a - 636)*x^5*y^8*z +
            (-200*a - 628)*x^8*y^4*z^2 + (1557*a - 114)*x^7*y^5*z^2 + (2197*a -
            2449)*x^6*y^6*z^2 + (1223*a - 3800)*x^5*y^7*z^2 + (343*a -
            1329)*x^4*y^8*z^2 + (-323*a - 809)*x^8*y^3*z^3 + (1630*a -
            631)*x^7*y^4*z^3 + (4190*a - 3126)*x^6*y^5*z^3 + (3904*a -
            7110)*x^5*y^6*z^3 + (1789*a - 5161)*x^4*y^7*z^3 + (330*a -
            1083)*x^3*y^8*z^3 + (-259*a - 524)*x^8*y^2*z^4 + (720*a -
            605)*x^7*y^3*z^4 + (3082*a - 2011)*x^6*y^4*z^4 + (4548*a -
            5462)*x^5*y^5*z^4 + (2958*a - 6611)*x^4*y^6*z^4 + (994*a -
            2931)*x^3*y^7*z^4 + (117*a - 416)*x^2*y^8*z^4 + (-108*a - 184)*x^8*y*z^5
            + (169*a - 168)*x^7*y^2*z^5 + (831*a - 835)*x^6*y^3*z^5 + (2225*a -
            1725)*x^5*y^4*z^5 + (1970*a - 3316)*x^4*y^5*z^5 + (952*a -
            2442)*x^3*y^6*z^5 + (217*a - 725)*x^2*y^7*z^5 + (16*a - 77)*x*y^8*z^5 +
            (-23*a - 35)*x^8*z^6 + (43*a + 24)*x^7*y*z^6 + (21*a - 198)*x^6*y^2*z^6
            + (377*a - 179)*x^5*y^3*z^6 + (458*a - 537)*x^4*y^4*z^6 + (288*a -
            624)*x^3*y^5*z^6 + (100*a - 299)*x^2*y^6*z^6 + (16*a - 67)*x*y^7*z^6 -
            5*y^8*z^6
              Defn: Defined on coordinates by sending (x : y : z) to
                    ((-5/128*a - 5/128)*x^4 + (-5/32*a + 5/32)*x^3*y + (-1/16*a +
            3/32)*x^2*y^2 + (1/16*a - 1/16)*x*y^3 + (1/32*a - 1/32)*y^4 - 1/32*x^3*z
            + (3/16*a - 5/8)*x^2*y*z + (1/8*a - 5/16)*x*y^2*z + (1/8*a +
            5/32)*x^2*z^2 + (-3/16*a + 5/16)*x*y*z^2 + (-3/16*a - 1/16)*y^2*z^2 +
            1/16*x*z^3 + (1/4*a + 1/4)*y*z^3 + (-3/32*a - 5/32)*z^4 : (-5/128*a -
            5/128)*x^4 + (5/32*a)*x^3*y + (3/32*a + 3/32)*x^2*y^2 + (-1/16*a)*x*y^3
            + (-1/32*a - 1/32)*y^4 - 1/32*x^3*z + (-11/32*a)*x^2*y*z + (1/8*a +
            5/16)*x*y^2*z + (3/16*a + 1/4)*y^3*z + (1/8*a + 5/32)*x^2*z^2 + (-1/16*a
            - 3/8)*x*y*z^2 + (-3/8*a - 9/16)*y^2*z^2 + 1/16*x*z^3 + (5/16*a +
            1/2)*y*z^3 + (-3/32*a - 5/32)*z^4 : (1/64*a + 3/128)*x^4 + (-1/32*a -
            1/32)*x^3*y + (3/32*a - 9/32)*x^2*y^2 + (1/16*a - 3/16)*x*y^3 - 1/32*y^4
            + (3/32*a + 1/8)*x^2*y*z + (-1/8*a + 1/8)*x*y^2*z + (-1/16*a)*y^3*z +
            (-1/16*a - 3/32)*x^2*z^2 + (1/16*a + 1/16)*x*y*z^2 + (3/16*a +
            3/16)*y^2*z^2 + (-3/16*a - 1/4)*y*z^3 + (1/16*a + 3/32)*z^4)
        """
        # helper function for extending the base field
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
        if not self.base_ring() in NumberFields():
            raise NotImplementedError("the base ring of this curve must be a number field")
        if not self.is_irreducible():
            raise TypeError("this curve must be irreducible")
        C_orig = self
        C = self
        PP = C.ambient_space()
        # extend the base field if necessary to find all singular points
        emb = extension(C.singular_subscheme())
        PP = PP.change_ring(emb)
        C = C.change_ring(emb)
        C_orig = C_orig.change_ring(emb)
        pts = C.singular_points()
        H = End(C)
        phi = H(list(C.ambient_space().gens()))
        while len(pts) > 0:
            for i in range(len(pts) - 1, -1, -1):
                try:
                    if C.is_ordinary_singularity(pts[i]):
                        pts.pop(i)
                except TypeError:
                    pts.pop(i)
            if len(pts) > 0:
                temp_exc = C.excellent_position(pts[0])
                temp_qua = temp_exc.codomain().quadratic_transform()
                C = temp_qua.codomain()
                phi = temp_qua*temp_exc*phi
                # transform the old points
                for i in range(len(pts) - 1, -1, -1):
                    # find image if it is a point the composition map is defined on
                    try:
                        temp_pt = (temp_qua*temp_exc)(temp_exc.domain()(pts[i]))
                        pts.pop(i)
                        if not PP(list(temp_pt)) in [PP(list(tpt)) for tpt in pts]:
                            pts.append(temp_pt)
                    except (TypeError, ValueError):
                        pass
                # add points from the intersection of C and the line z
                PPline = ProjectiveSpace(PP.base_ring(), 1)
                # make sure the conversion happens in the right order
                ringH = Hom(PP.coordinate_ring(), PPline.coordinate_ring())
                psi = ringH(list(PPline.gens()) + [0])
                X = PPline.subscheme([psi(f) for f in C.singular_subscheme().defining_polynomials()])
                emb = extension(X)
                PP = PP.change_ring(emb)
                phi = phi.change_ring(emb)
                C = C.change_ring(emb)
                C_orig = C_orig.change_ring(emb)
                X = X.change_ring(emb)
                pts = [PP(pt.change_ring(emb)) for pt in pts]
                newpts = [PP(list(pt) + [0]) for pt in X.rational_points()]
                # avoid duplicates
                for pt in newpts:
                    if not PP(list(pt)) in [PP(list(tpt)) for tpt in pts]:
                        pts.append(pt)
        return phi

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

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: C = Curve([x^2 - y^2], P)
            sage: D = Curve([x - y], P)
            sage: Q = P([1,1,0])
            sage: C.is_transverse(D, Q)
            False

        ::

            sage: K = QuadraticField(-1)
            sage: P.<x,y,z> = ProjectiveSpace(K, 2)
            sage: C = Curve([y^2*z - K.0*x^3], P)
            sage: D = Curve([z*x + y^2], P)
            sage: Q = P([0,0,1])
            sage: C.is_transverse(D, Q)
            False

        ::

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: C = Curve([x^2 - 2*y^2 - 2*z^2], P)
            sage: D = Curve([y - z], P)
            sage: Q = P([2,1,1])
            sage: C.is_transverse(D, Q)
            True
        """
        if not self.intersects_at(C, P):
            raise TypeError("(=%s) must be a point in the intersection of (=%s) and this curve"%(P,C))
        if self.is_singular(P) or C.is_singular(P):
            return False

        # there is only one tangent at a nonsingular point of a plane curve
        return not self.tangents(P)[0] == C.tangents(P)[0]

    def rational_parameterization(self):
        r"""
        Return a rational parameterization of this curve.

        This curve must have rational coefficients and be absolutely irreducible (i.e. irreducible
        over the algebraic closure of the rational field). The curve must also be rational (have
        geometric genus zero).

        The rational parameterization may have coefficients in a quadratic extension of the rational
        field.

        OUTPUT:

        - a birational map between `\mathbb{P}^{1}` and this curve, given as a scheme morphism.

        EXAMPLES::

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: C = Curve([y^2*z - x^3], P)
            sage: C.rational_parameterization()
            Scheme morphism:
              From: Projective Space of dimension 1 over Rational Field
              To:   Projective Plane Curve over Rational Field defined by -x^3 + y^2*z
              Defn: Defined on coordinates by sending (s : t) to
                    (s^2*t : s^3 : t^3)

        ::

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: C = Curve([x^3 - 4*y*z^2 + x*z^2 - x*y*z], P)
            sage: C.rational_parameterization()
            Scheme morphism:
              From: Projective Space of dimension 1 over Rational Field
              To:   Projective Plane Curve over Rational Field defined by x^3 - x*y*z + x*z^2 - 4*y*z^2
              Defn: Defined on coordinates by sending (s : t) to
                    (4*s^2*t + s*t^2 : s^2*t + t^3 : 4*s^3 + s^2*t)

        ::

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: C = Curve([x^2 + y^2 + z^2], P)
            sage: C.rational_parameterization()
            Scheme morphism:
              From: Projective Space of dimension 1 over Number Field in a with defining polynomial a^2 + 1
              To:   Projective Plane Curve over Number Field in a with defining
              polynomial a^2 + 1 defined by x^2 + y^2 + z^2
              Defn: Defined on coordinates by sending (s : t) to
                    ((-a)*s^2 + (-a)*t^2 : s^2 - t^2 : 2*s*t)
        """
        if self.genus() != 0:
            raise TypeError("this curve must have geometric genus zero")
        if not is_RationalField(self.base_ring()):
            raise TypeError("this curve must be defined over the rational field")
        singular.lib("paraplanecurves.lib")
        R = singular.paraPlaneCurve(self.defining_polynomial())
        singular.setring(R)
        param = singular('PARA').sage().gens()
        R = R.sage()
        C = self.change_ring(R.base_ring())
        H = Hom(ProjectiveSpace(R.base_ring(), 1, R.gens()), C)
        return H(param)

class ProjectivePlaneCurve_finite_field(ProjectivePlaneCurve):

    _point = point.ProjectivePlaneCurvePoint_finite_field

    def rational_points_iterator(self):
        r"""
        Return a generator object for the rational points on this curve.

        INPUT:

        - ``self`` -- a projective curve

        OUTPUT:

        A generator of all the rational points on the curve defined over its base field.

        EXAMPLE::

            sage: F = GF(37)
            sage: P2.<X,Y,Z> = ProjectiveSpace(F,2)
            sage: C = Curve(X^7+Y*X*Z^5*55+Y^7*12)
            sage: len(list(C.rational_points_iterator()))
            37

        ::

            sage: F = GF(2)
            sage: P2.<X,Y,Z> = ProjectiveSpace(F,2)
            sage: C = Curve(X*Y*Z)
            sage: a = C.rational_points_iterator()
            sage: next(a)
            (1 : 0 : 0)
            sage: next(a)
            (0 : 1 : 0)
            sage: next(a)
            (1 : 1 : 0)
            sage: next(a)
            (0 : 0 : 1)
            sage: next(a)
            (1 : 0 : 1)
            sage: next(a)
            (0 : 1 : 1)
            sage: next(a)
            Traceback (most recent call last):
            ...
            StopIteration

        ::

            sage: F = GF(3^2,'a')
            sage: P2.<X,Y,Z> = ProjectiveSpace(F,2)
            sage: C = Curve(X^3+5*Y^2*Z-33*X*Y*X)
            sage: b = C.rational_points_iterator()
            sage: next(b)
            (0 : 1 : 0)
            sage: next(b)
            (0 : 0 : 1)
            sage: next(b)
            (2*a + 2 : a : 1)
            sage: next(b)
            (2 : a + 1 : 1)
            sage: next(b)
            (a + 1 : 2*a + 1 : 1)
            sage: next(b)
            (1 : 2 : 1)
            sage: next(b)
            (2*a + 2 : 2*a : 1)
            sage: next(b)
            (2 : 2*a + 2 : 1)
            sage: next(b)
            (a + 1 : a + 2 : 1)
            sage: next(b)
            (1 : 1 : 1)
            sage: next(b)
            Traceback (most recent call last):
            ...
            StopIteration

        """
        g = self.defining_polynomial()
        K = g.parent().base_ring()
        from sage.rings.polynomial.all import PolynomialRing
        R = PolynomialRing(K,'X')
        X = R.gen()
        one = K.one()
        zero = K.zero()

        # the point with  Z = 0 = Y
        try:
            t = self.point([one,zero,zero])
            yield(t)
        except TypeError:
            pass

        # points with Z = 0, Y = 1
        g10 = R(g(X,one,zero))
        if g10.is_zero():
            for x in K:
                yield(self.point([x,one,zero]))
        else:
            for x in g10.roots(multiplicities=False):
                yield(self.point([x,one,zero]))

        # points with Z = 1
        for y in K:
            gy1 = R(g(X,y,one))
            if gy1.is_zero():
                for x in K:
                    yield(self.point([x,y,one]))
            else:
                for x in gy1.roots(multiplicities=False):
                    yield(self.point([x,y,one]))

    def rational_points(self, algorithm="enum", sort=True):
        r"""
        Return the rational points on this curve computed via enumeration.


        INPUT:

        - ``algorithm`` (string, default: 'enum') -- the algorithm to
          use.  Currently this is ignored.

        - ``sort`` (boolean, default ``True``) -- whether the output
          points should be sorted.  If False, the order of the output
          is non-deterministic.

        OUTPUT:

        A list of all the rational points on the curve defined over
        its base field, possibly sorted.

        .. note::

           This is a slow Python-level implementation.


        EXAMPLES::

            sage: F = GF(7)
            sage: P2.<X,Y,Z> = ProjectiveSpace(F,2)
            sage: C = Curve(X^3+Y^3-Z^3)
            sage: C.rational_points()
            [(0 : 1 : 1), (0 : 2 : 1), (0 : 4 : 1), (1 : 0 : 1), (2 : 0 : 1), (3 : 1 : 0), (4 : 0 : 1), (5 : 1 : 0), (6 : 1 : 0)]


        ::

            sage: F = GF(1237)
            sage: P2.<X,Y,Z> = ProjectiveSpace(F,2)
            sage: C = Curve(X^7+7*Y^6*Z+Z^4*X^2*Y*89)
            sage: len(C.rational_points())
            1237

        ::

            sage: F = GF(2^6,'a')
            sage: P2.<X,Y,Z> = ProjectiveSpace(F,2)
            sage: C = Curve(X^5+11*X*Y*Z^3 + X^2*Y^3 - 13*Y^2*Z^3)
            sage: len(C.rational_points())
            104

        ::

            sage: R.<x,y,z> = GF(2)[]
            sage: f = x^3*y + y^3*z + x*z^3
            sage: C = Curve(f); pts = C.rational_points()
            sage: pts
            [(0 : 0 : 1), (0 : 1 : 0), (1 : 0 : 0)]

        """
        points = list(self.rational_points_iterator())
        if sort:
            points.sort()
        return points

class ProjectivePlaneCurve_prime_finite_field(ProjectivePlaneCurve_finite_field):
    def _points_via_singular(self, sort=True):
        r"""
        Return all rational points on this curve, computed using Singular's
        Brill-Noether implementation.

        INPUT:


        -  ``sort`` - bool (default: True), if True return the
           point list sorted. If False, returns the points in the order
           computed by Singular.


        EXAMPLE::

            sage: x, y, z = PolynomialRing(GF(5), 3, 'xyz').gens()
            sage: f = y^2*z^7 - x^9 - x*z^8
            sage: C = Curve(f); C
            Projective Plane Curve over Finite Field of size 5 defined by
            -x^9 + y^2*z^7 - x*z^8
            sage: C._points_via_singular()
            [(0 : 0 : 1), (0 : 1 : 0), (2 : 2 : 1), (2 : 3 : 1),
             (3 : 1 : 1), (3 : 4 : 1)]
            sage: C._points_via_singular(sort=False)     #random
            [(0 : 1 : 0), (3 : 1 : 1), (3 : 4 : 1), (2 : 2 : 1),
             (0 : 0 : 1), (2 : 3 : 1)]


        .. note::

           The Brill-Noether package does not always work (i.e., the
           'bn' algorithm. When it fails a RuntimeError exception is
           raised.
        """
        f = self.defining_polynomial()._singular_()
        singular = f.parent()
        singular.lib('brnoeth')
        try:
            X1 = f.Adj_div()
        except (TypeError, RuntimeError) as s:
            raise RuntimeError(str(s) + "\n\n ** Unable to use the\
                                          Brill-Noether Singular package to\
                                          compute all points (see above).")

        X2 = singular.NSplaces(1, X1)
        R = X2[5][1][1]
        singular.set_ring(R)

        # We use sage_flattened_str_list since iterating through
        # the entire list through the sage/singular interface directly
        # would involve hundreds of calls to singular, and timing issues with
        # the expect interface could crop up.  Also, this is vastly
        # faster (and more robust).
        v = singular('POINTS').sage_flattened_str_list()
        pnts = [self(int(v[3*i]), int(v[3*i+1]), int(v[3*i+2]))
                for i in range(len(v)//3)]
        # singular always dehomogenizes with respect to the last variable
        # so if this variable divides the curve equation, we need to add
        # points at infinity
        F = self.defining_polynomial()
        z = F.parent().gens()[-1]
        if z.divides(F):
            pnts += [self(1,a,0) for a in self.base_ring()]
            pnts += [self(0,1,0)]
        # remove multiple points
        pnts = list(set(pnts))
        if sort:
            pnts.sort()
        return pnts

    def riemann_roch_basis(self, D):
        r"""
        Return a basis for the Riemann-Roch space corresponding to
        `D`.

        This uses Singular's Brill-Noether implementation.

        INPUT:

        -  ``D`` - a divisor

        OUTPUT:

        A list of function field elements that form a basis of the Riemann-Roch space

        EXAMPLE::

            sage: R.<x,y,z> = GF(2)[]
            sage: f = x^3*y + y^3*z + x*z^3
            sage: C = Curve(f); pts = C.rational_points()
            sage: D = C.divisor([ (4, pts[0]), (4, pts[2]) ])
            sage: C.riemann_roch_basis(D)
            [x/y, 1, z/y, z^2/y^2, z/x, z^2/(x*y)]

        ::

            sage: R.<x,y,z> = GF(5)[]
            sage: f = x^7 + y^7 + z^7
            sage: C = Curve(f); pts = C.rational_points()
            sage: D = C.divisor([ (3, pts[0]), (-1,pts[1]), (10, pts[5]) ])
            sage: C.riemann_roch_basis(D)
            [(-2*x + y)/(x + y), (-x + z)/(x + y)]


        .. NOTE::

            Currently this only works over prime field and divisors
            supported on rational points.
        """
        f = self.defining_polynomial()._singular_()
        singular = f.parent()
        singular.lib('brnoeth')
        try:
            X1 = f.Adj_div()
        except (TypeError, RuntimeError) as s:
            raise RuntimeError(str(s) + "\n\n ** Unable to use the Brill-Noether Singular package to compute all points (see above).")
        X2 = singular.NSplaces(1, X1)
        # retrieve list of all computed closed points (possibly of degree >1)
        v = X2[3].sage_flattened_str_list()    # We use sage_flattened_str_list since iterating through
                                               # the entire list through the sage/singular interface directly
                                               # would involve hundreds of calls to singular, and timing issues with
                                               # the expect interface could crop up.  Also, this is vastly
                                               # faster (and more robust).
        v = [ v[i].partition(',') for i in range(len(v)) ]
        pnts = [ ( int(v[i][0]), int(v[i][2])-1 ) for i in range(len(v))]
        # retrieve coordinates of rational points
        R = X2[5][1][1]
        singular.set_ring(R)
        v = singular('POINTS').sage_flattened_str_list()
        coords = [self(int(v[3*i]), int(v[3*i+1]), int(v[3*i+2])) for i in range(len(v)//3)]
        # build correct representation of D for singular
        Dsupport = D.support()
        Dcoeffs = []
        for x in pnts:
            if x[0] == 1:
                Dcoeffs.append(D.coefficient(coords[x[1]]))
            else:
                Dcoeffs.append(0)
        Dstr = str(tuple(Dcoeffs))
        G = singular(','.join([str(x) for x in Dcoeffs]), type='intvec')
        # call singular's brill noether routine and return
        T = X2[1][2]
        T.set_ring()
        LG = G.BrillNoether(X2)
        LG = [X.split(',\n') for X in LG.sage_structured_str_list()]
        x,y,z = self.ambient_space().coordinate_ring().gens()
        vars = {'x':x, 'y':y, 'z':z}
        V = [(sage_eval(a, vars)/sage_eval(b, vars)) for a, b in LG]
        return V

    def rational_points(self, algorithm="enum", sort=True):
        r"""
        INPUT:


        -  ``algorithm`` - string:

        -  ``'enum'`` - straightforward enumeration

        -  ``'bn'`` - via Singular's brnoeth package.


        EXAMPLE::

            sage: x, y, z = PolynomialRing(GF(5), 3, 'xyz').gens()
            sage: f = y^2*z^7 - x^9 - x*z^8
            sage: C = Curve(f); C
            Projective Plane Curve over Finite Field of size 5 defined by
            -x^9 + y^2*z^7 - x*z^8
            sage: C.rational_points()
            [(0 : 0 : 1), (0 : 1 : 0), (2 : 2 : 1), (2 : 3 : 1),
             (3 : 1 : 1), (3 : 4 : 1)]
            sage: C = Curve(x - y + z)
            sage: C.rational_points()
            [(0 : 1 : 1), (1 : 1 : 0), (1 : 2 : 1), (2 : 3 : 1),
             (3 : 4 : 1), (4 : 0 : 1)]
            sage: C = Curve(x*z+z^2)
            sage: C.rational_points('all')
            [(0 : 1 : 0), (1 : 0 : 0), (1 : 1 : 0), (2 : 1 : 0),
             (3 : 1 : 0), (4 : 0 : 1), (4 : 1 : 0), (4 : 1 : 1),
             (4 : 2 : 1), (4 : 3 : 1), (4 : 4 : 1)]

        .. note::

           The Brill-Noether package does not always work (i.e., the
           'bn' algorithm. When it fails a RuntimeError exception is
           raised.
        """
        if algorithm == "enum":

            return ProjectivePlaneCurve_finite_field.rational_points(self,
                                                                algorithm="enum",
                                                                sort=sort)

        elif algorithm == "bn":

            return self._points_via_singular(sort=sort)

        elif algorithm == "all":

            S_enum = self.rational_points(algorithm = "enum")
            S_bn = self.rational_points(algorithm = "bn")
            if S_enum != S_bn:
                raise RuntimeError("Bug in rational_points -- different\
                                     algorithms give different answers for\
                                     curve %s!"%self)
            return S_enum

        else:

            raise ValueError("No algorithm '%s' known"%algorithm)

def Hasse_bounds(q, genus=1):
    r"""
    Return the Hasse-Weil bounds for the cardinality of a nonsingular
    curve defined over `\GF{q}` of given ``genus``.

    INPUT:

    - ``q`` (int) -- a prime power

    - ``genus`` (int, default 1) -- a non-negative integer,

    OUTPUT:

    (tuple)  The Hasse bounds (lb,ub) for the cardinality of a curve of
    genus ``genus`` defined over `\GF{q}`.

    EXAMPLES::

        sage: Hasse_bounds(2)
        (1, 5)
        sage: Hasse_bounds(next_prime(10^30))
        (999999999999998000000000000058, 1000000000000002000000000000058)
    """
    if genus==1:
        rq = (4*q).isqrt()
    else:
        rq = (4*(genus**2)*q).isqrt()
    return (q+1-rq,q+1+rq)

# Fix pickles from changing class names and plane_curves folder name
from sage.structure.sage_object import register_unpickle_override
register_unpickle_override('sage.schemes.plane_curves.projective_curve',
                           'ProjectiveCurve_generic', ProjectivePlaneCurve)
