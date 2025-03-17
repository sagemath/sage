# sage.doctest: needs sage.libs.singular
r"""
Affine curves

Affine curves in Sage are curves in an affine space or an affine plane.

EXAMPLES:

We can construct curves in either an affine plane::

    sage: A.<x,y> = AffineSpace(QQ, 2)
    sage: C = Curve([y - x^2], A); C
    Affine Plane Curve over Rational Field defined by -x^2 + y

or in higher dimensional affine space::

    sage: A.<x,y,z,w> = AffineSpace(QQ, 4)
    sage: C = Curve([y - x^2, z - w^3, w - y^4], A); C
    Affine Curve over Rational Field defined by -x^2 + y, -w^3 + z, -y^4 + w

Integral affine curves over finite fields
-----------------------------------------

If the curve is defined over a finite field and integral, that is reduced and
irreducible, its function field is tightly coupled with the curve so that
advanced computations based on Sage's global function field machinery are
available.

EXAMPLES::

    sage: k.<a> = GF(2)
    sage: A.<x,y,z> = AffineSpace(k, 3)
    sage: C = Curve([x^2 + x - y^3, y^4 - y - z^3], A)
    sage: C.genus()
    10
    sage: C.function_field()
    Function field in z defined by z^9 + x^8 + x^6 + x^5 + x^4 + x^3 + x

Closed points of arbitrary degree can be computed::

    sage: # long time
    sage: C.closed_points()
    [Point (x, y, z), Point (x + 1, y, z)]
    sage: C.closed_points(2)
    [Point (x^2 + x + 1, y + 1, z),
     Point (y^2 + y + 1, x + y, z),
     Point (y^2 + y + 1, x + y + 1, z)]
    sage: p = _[0]
    sage: p.places()
    [Place (x^2 + x + 1, (1/(x^4 + x^2 + 1))*z^7 + (1/(x^4 + x^2 + 1))*z^6 + 1)]

The places at infinity correspond to the extra closed points of the curve's
projective closure::

    sage: C.places_at_infinity()                # long time
    [Place (1/x, 1/x*z)]

It is easy to transit to and from the function field of the curve::

    sage: fx = C(x)
    sage: fy = C(y)
    sage: fx^2 + fx - fy^3
    0
    sage: fx.divisor()
    -9*Place (1/x, 1/x*z)
     + 9*Place (x, z)
    sage: p, = fx.zeros()
    sage: C.place_to_closed_point(p)
    Point (x, y, z)
    sage: _.rational_point()
    (0, 0, 0)
    sage: _.closed_point()
    Point (x, y, z)
    sage: _.place()
    Place (x, z)

Integral affine curves over `\QQ`
---------------------------------

An integral curve over `\QQ` is equipped also with the function field. Unlike
over finite fields, it is not possible to enumerate closed points.

EXAMPLES::

    sage: A.<x,y> = AffineSpace(QQ, 2)
    sage: C = Curve(x^2 + y^2 -1)
    sage: p = C(0,1)
    sage: p
    (0, 1)
    sage: p.closed_point()
    Point (x, y - 1)
    sage: pl = _.place()
    sage: C.parametric_representation(pl)
    (s + ..., 1 - 1/2*s^2 - 1/8*s^4 - 1/16*s^6 + ...)
    sage: sx, sy = _
    sage: sx = sx.polynomial(10); sx
    s
    sage: sy = sy.polynomial(10); sy
    -7/256*s^10 - 5/128*s^8 - 1/16*s^6 - 1/8*s^4 - 1/2*s^2 + 1
    sage: s = var('s')                                                                  # needs sage.symbolic
    sage: P1 = parametric_plot([sx, sy], (s, -1, 1), color='red')                       # needs sage.plot sage.symbolic
    sage: P2 = C.plot((x, -1, 1), (y, 0, 2))  # half circle                             # needs sage.plot sage.symbolic
    sage: P1 + P2                                                                       # needs sage.plot sage.symbolic
    Graphics object consisting of 2 graphics primitives

AUTHORS:

- William Stein (2005-11-13)

- David Joyner (2005-11-13)

- David Kohel (2006-01)

- Grayson Jorgenson (2016-08)

- Kwankyu Lee (2019-05): added integral affine curves
"""
# ****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from builtins import sum as add

import sage.rings.abc

from sage.arith.misc import binomial
from sage.categories.fields import Fields
from sage.categories.finite_fields import FiniteFields
from sage.categories.homset import Hom, End, hom
from sage.categories.number_fields import NumberFields
from sage.matrix.constructor import matrix
from sage.misc.cachefunc import cached_method
from sage.misc.lazy_attribute import lazy_attribute
from sage.misc.lazy_import import lazy_import
from sage.rings.infinity import infinity
from sage.rings.polynomial.multi_polynomial_element import degree_lowest_rational_function
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.rational_field import RationalField

from sage.schemes.affine.affine_space import AffineSpace, AffineSpace_generic
from sage.schemes.affine.affine_subscheme import (AlgebraicScheme_subscheme_affine,
                                                  AlgebraicScheme_subscheme_affine_field)

lazy_import('sage.interfaces.singular', 'singular')
lazy_import('sage.rings.number_field.number_field', 'NumberField')
lazy_import('sage.rings.qqbar', 'number_field_elements_from_algebraics')

from .curve import Curve_generic

from .point import (AffineCurvePoint_field,
                    AffinePlaneCurvePoint_field,
                    AffinePlaneCurvePoint_finite_field,
                    IntegralAffineCurvePoint,
                    IntegralAffineCurvePoint_finite_field,
                    IntegralAffinePlaneCurvePoint,
                    IntegralAffinePlaneCurvePoint_finite_field)

from .closed_point import IntegralAffineCurveClosedPoint


class AffineCurve(Curve_generic, AlgebraicScheme_subscheme_affine):
    """
    Affine curves.

    EXAMPLES::

        sage: # needs sage.rings.number_field
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

    def __init__(self, A, X):
        r"""
        Initialize.

        EXAMPLES::

            sage: # needs sage.rings.number_field
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
        if not isinstance(A, AffineSpace_generic):
            raise TypeError("A (={}) must be an affine space".format(A))

        Curve_generic.__init__(self, A, X)

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

    def projective_closure(self, i=0, PP=None):
        r"""
        Return the projective closure of this affine curve.

        INPUT:

        - ``i`` -- (default: 0) the index of the affine coordinate chart of the
          projective space that the affine ambient space of this curve embeds into

        - ``PP`` -- (default: ``None``) ambient projective space to compute the
          projective closure in; this is constructed if it is not given

        OUTPUT: a curve in projective space

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


class AffinePlaneCurve(AffineCurve):
    """
    Affine plane curves.
    """

    def __init__(self, A, f):
        r"""
        Initialize.

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
        if not (isinstance(A, AffineSpace_generic) and A.dimension != 2):
            raise TypeError("Argument A (= %s) must be an affine plane." % A)

        super().__init__(A, [f])

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

        INPUT:

        - ``r`` -- a rational function on X

        OUTPUT:

        - ``list`` -- the divisor of r represented as a list of coefficients
          and points. (TODO: This will change to a more structural output in
          the future.)

        EXAMPLES::

            sage: F = GF(5)
            sage: P2 = AffineSpace(2, F, names='xy')
            sage: R = P2.coordinate_ring()
            sage: x, y = R.gens()
            sage: f = y^2 - x^9 - x
            sage: C = Curve(f)
            sage: K = FractionField(R)
            sage: r = 1/x
            sage: C.divisor_of_function(r)      # not implemented (broken)
                  [[-1, (0, 0, 1)]]
            sage: r = 1/x^3
            sage: C.divisor_of_function(r)      # not implemented (broken)
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

        - ``pt`` -- an F-rational point on X which is not a
          point of ramification for the projection (x,y) - x

        - ``n`` -- the number of terms desired

        OUTPUT: x = x0 + t y = y0 + power series in t

        EXAMPLES::

            sage: F = GF(5)
            sage: pt = (2,3)
            sage: R = PolynomialRing(F, 2, names = ['x','y'])
            sage: x,y = R.gens()
            sage: f = y^2 - x^9 - x
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
        astr = ["a" + str(i) for i in range(1, 2 * n)]
        x, y = R.gens()
        R0 = PolynomialRing(F, 2 * n + 2, names=[str(x), str(y), "t"] + astr)
        vars0 = R0.gens()
        t = vars0[2]
        yt = y0*t**0+add([vars0[i]*t**(i-2) for i in range(3, 2*n+2)])
        xt = x0+t
        ft = f(xt, yt)
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
                    if x.replace(str(y), ""):
                        i = x.find("-")
                        if i > 0:
                            vals.append([eval(x[1:i]), x[:i], F(eval(x[i+1:]))])
                        i = x.find("+")
                        if i > 0:
                            vals.append([eval(x[1:i]), x[:i], -F(eval(x[i+1:]))])
                    else:
                        vals.append([eval(str(y)[1:]), str(y), F(0)])
        vals.sort()
        return [x0 + t, y0 + add(v[2] * t**(j + 1) for j, v in enumerate(vals))]

    def plot(self, *args, **kwds):
        r"""
        Plot the real points on this affine plane curve.

        INPUT:

        - ``*args`` -- (optional) tuples (variable, minimum, maximum) for
          plotting dimensions

        - ``**kwds`` -- optional keyword arguments passed on to ``implicit_plot``

        EXAMPLES:

        A cuspidal curve::

            sage: R.<x, y> = QQ[]
            sage: C = Curve(x^3 - y^2)
            sage: C.plot()                                      # needs sage.plot
            Graphics object consisting of 1 graphics primitive

        A 5-nodal curve of degree 11.  This example also illustrates
        some of the optional arguments::

            sage: # needs sage.plot
            sage: R.<x, y> = ZZ[]
            sage: C = Curve(32*x^2 - 2097152*y^11 + 1441792*y^9
            ....:            - 360448*y^7 + 39424*y^5 - 1760*y^3 + 22*y - 1)
            sage: C.plot((x, -1, 1), (y, -1, 1), plot_points=400)
            Graphics object consisting of 1 graphics primitive

        A line over `\mathbf{RR}`::

            sage: # needs sage.symbolic sage.plot
            sage: R.<x, y> = RR[]
            sage: C = Curve(R(y - sqrt(2)*x))
            sage: C.plot()
            Graphics object consisting of 1 graphics primitive
        """
        Id = self.defining_ideal()
        return Id.plot(*args, **kwds)

    def is_transverse(self, C, P):
        r"""
        Return whether the intersection of this curve with the curve ``C`` at the point ``P`` is transverse.

        The intersection at ``P`` is transverse if ``P`` is a nonsingular point of both curves, and if the
        tangents of the curves at ``P`` are distinct.

        INPUT:

        - ``C`` -- a curve in the ambient space of this curve

        - ``P`` -- a point in the intersection of both curves

        OUTPUT: boolean

        EXAMPLES::

            sage: A.<x,y> = AffineSpace(QQ, 2)
            sage: C = Curve([x^2 + y^2 - 1], A)
            sage: D = Curve([x - 1], A)
            sage: Q = A([1,0])
            sage: C.is_transverse(D, Q)
            False

        ::

            sage: # needs sage.rings.number_field
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

        - ``P`` -- a point in the ambient space of this curve

        OUTPUT: integer

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

            sage: # needs sage.rings.number_field
            sage: A.<x,y> = AffineSpace(QQbar,2)
            sage: C = Curve([-x^7 + (-7)*x^6 + y^6 + (-21)*x^5 + 12*y^5
            ....:            + (-35)*x^4 + 60*y^4 + (-35)*x^3 + 160*y^3
            ....:            + (-21)*x^2 + 240*y^2 + (-7)*x + 192*y + 63], A)
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
        if self.base_ring() not in Fields():
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

        - ``P`` -- a point on this curve

        - ``factor`` -- boolean (default: ``True``); whether to attempt computing the
          polynomials of the individual tangent lines over the base field of this
          curve, or to just return the polynomial corresponding to the union of
          the tangent lines (which requires fewer computations)

        OUTPUT: list of polynomials in the coordinate ring of the ambient space

        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: set_verbose(-1)
            sage: A.<x,y> = AffineSpace(QQbar, 2)
            sage: C = Curve([x^5*y^3 + 2*x^4*y^4 + x^3*y^5 + 3*x^4*y^3
            ....:            + 6*x^3*y^4 + 3*x^2*y^5 + 3*x^3*y^3
            ....:            + 6*x^2*y^4 + 3*x*y^5 + x^5 + 10*x^4*y
            ....:            + 40*x^3*y^2 + 81*x^2*y^3 + 82*x*y^4 + 33*y^5], A)
            sage: Q = A([0,0])
            sage: C.tangents(Q)
            [x + 3.425299577684700?*y,
             x + (1.949159013086856? + 1.179307909383728?*I)*y,
             x + (1.949159013086856? - 1.179307909383728?*I)*y,
             x + (1.338191198070795? + 0.2560234251008043?*I)*y,
             x + (1.338191198070795? - 0.2560234251008043?*I)*y]
            sage: C.tangents(Q, factor=False)
            [120*x^5 + 1200*x^4*y + 4800*x^3*y^2 + 9720*x^2*y^3 + 9840*x*y^4 + 3960*y^5]

        ::

            sage: # needs sage.rings.number_field
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
        coords = [vars[0] - P[0], vars[1] - P[1]]  # coords to change back with
        deriv = [f.derivative(vars[0], i).derivative(vars[1], r-i)([0, 0]) for i in range(r+1)]
        T = sum([binomial(r, i)*deriv[i]*(vars[0])**i*(vars[1])**(r-i) for i in range(r+1)])
        if not factor:
            return [T(coords)]
        if isinstance(self.base_ring(), sage.rings.abc.AlgebraicField):
            fact = []
            # first add tangents corresponding to vars[0], vars[1] if they divide T
            t = min([e[0] for e in T.exponents()])
            # vars[0] divides T
            if t > 0:
                fact.append(vars[0])
                # divide T by that power of vars[0]
                T = self.ambient_space().coordinate_ring()({(v[0] - t, v[1]): h for (v, h) in T.dict().items()})
            t = min([e[1] for e in T.exponents()])
            # vars[1] divides T
            if t > 0:
                fact.append(vars[1])
                # divide T by that power of vars[1]
                T = self.ambient_space().coordinate_ring()({(v[0], v[1] - t): h for (v, h) in T.dict().items()})
            # T is homogeneous in var[0], var[1] if nonconstant, so dehomogenize
            if T not in self.base_ring():
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
            return [ll[0](coords) for ll in T.factor()]

    def is_ordinary_singularity(self, P):
        r"""
        Return whether the singular point ``P`` of this affine plane curve is
        an ordinary singularity.

        The point ``P`` is an ordinary singularity of this curve if it is a
        singular point, and if the tangents of this curve at ``P`` are
        distinct.

        INPUT:

        - ``P`` -- a point on this curve

        OUTPUT:

        ``True`` or ``False`` depending on whether ``P`` is or is not an ordinary
        singularity of this curve, respectively. An error is raised if ``P`` is
        not a singular point of this curve.

        EXAMPLES::

            sage: A.<x,y> = AffineSpace(QQ, 2)
            sage: C = Curve([y^2 - x^3], A)
            sage: Q = A([0,0])
            sage: C.is_ordinary_singularity(Q)
            False

        ::

            sage: # needs sage.rings.number_field
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
            raise TypeError("(=%s) is not a singular point of (=%s)" % (P, self))

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

        OUTPUT: a birational map between `\mathbb{A}^{1}` and this curve, given
        as a scheme morphism

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


class AffineCurve_field(AffineCurve, AlgebraicScheme_subscheme_affine_field):
    """
    Affine curves over fields.
    """
    _point = AffineCurvePoint_field

    def __init__(self, A, X):
        r"""
        Initialize.

        EXAMPLES::

            sage: # needs sage.rings.number_field
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

        TESTS::

            sage: K.<x,y,z,t> = QQ[]
            sage: t1 = x^2*z^2 + y*t
            sage: t2 = y*z^2 + x^2*t
            sage: C = Curve([x^4 - y^2 - 19, z^4 - t^2 - 23, t1^2 - t2^2 - 19*23])
            Traceback (most recent call last):
            ...
            ValueError: defining equations (=[x^4 - y^2 - 19, z^4 - t^2 - 23,
            x^4*z^4 - y^2*z^4 - x^4*t^2 + y^2*t^2 - 437]) define a scheme of dimension 2 != 1
        """
        super().__init__(A, X)

        if A.base_ring() not in Fields():
            raise TypeError("curve not defined over a field")

        d = super(Curve_generic, self).dimension()
        if d != 1:
            raise ValueError("defining equations (={}) define a scheme of dimension {} != 1".format(X, d))

    def projection(self, indices, AS=None):
        r"""
        Return the projection of this curve onto the coordinates specified by
        ``indices``.

        INPUT:

        - ``indices`` -- list or tuple of distinct integers specifying the
          indices of the coordinates to use in the projection. Can also be a list
          or tuple consisting of variables of the coordinate ring of the ambient
          space of this curve. If integers are used to specify the coordinates, 0
          denotes the first coordinate. The length of ``indices`` must be between
          two and one less than the dimension of the ambient space of this curve,
          inclusive.

        - ``AS`` -- (default: ``None``) the affine space the projected curve will
          be defined in. This space must be defined over the same base field as
          this curve, and must have dimension equal to the length of ``indices``.
          This space is constructed if not specified.

        OUTPUT: a tuple of

        - a scheme morphism from this curve to affine space of dimension equal
          to the number of coordinates specified in ``indices``

        - the affine subscheme that is the image of that morphism. If the image
          is a curve, the second element of the tuple will be a curve.

        EXAMPLES::

            sage: A.<x,y,z> = AffineSpace(QQ, 3)
            sage: C = Curve([y^7 - x^2 + x^3 - 2*z, z^2 - x^7 - y^2], A)
            sage: C.projection([0,1])
            (Scheme morphism:
               From: Affine Curve over Rational Field
                     defined by y^7 + x^3 - x^2 - 2*z, -x^7 - y^2 + z^2
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
             Closed subscheme of Affine Space of dimension 3 over Rational Field defined by:
               c - 1,
               b - 3,
               a - 2)

        ::

            sage: A.<x,y,z,w,u> = AffineSpace(GF(11), 5)
            sage: C = Curve([x^3 - 5*y*z + u^2, x - y^2 + 3*z^2,
            ....:            w^2 + 2*u^3*y, y - u^2 + z*x], A)
            sage: B.<a,b,c> = AffineSpace(GF(11), 3)
            sage: proj1 = C.projection([1,2,4], AS=B); proj1
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
               From: Affine Curve over Rational Field defined by
                     -z^3 + x*y, -w^3 + x*z, -x^3 + w^2
               To:   Affine Space of dimension 2 over Rational Field
               Defn: Defined on coordinates by sending (x, y, z, w) to (y, z),
             Affine Plane Curve over Rational Field defined by x1^23 - x0^7*x1^4)
            sage: B.<x,y,z> = AffineSpace(QQ, 3)
            sage: C.projection([x,y,z], AS=B)
            (Scheme morphism:
               From: Affine Curve over Rational Field defined by
                     -z^3 + x*y, -w^3 + x*z, -x^3 + w^2
               To:   Affine Space of dimension 3 over Rational Field
               Defn: Defined on coordinates by sending (x, y, z, w) to
                     (x, y, z),
             Affine Curve over Rational Field defined by
             z^3 - x*y, x^8 - x*z^2, x^7*z^2 - x*y*z)
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
        if AS is not None:
            if not isinstance(AS, AffineSpace_generic):
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
            indices = [int(i) for i in indices]  # type checking
            indices.sort()
            if indices[0] < 0 or indices[-1] > n - 1:
                raise ValueError("index values must be between 0 and one "
                                 "minus the dimension of the ambient space "
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
        ll = [0]*(n)
        for i in range(len(indices)):
            ll[indices[i]] = AA2.gens()[i]
        phi = K(ll)
        G = [phi(f) for f in J.gens()]
        try:
            C = AA2.curve(G)
        except (TypeError, ValueError):
            C = AA2.subscheme(G)
        return (psi, C)

    def plane_projection(self, AP=None):
        r"""
        Return a projection of this curve into an affine plane so that the
        image of the projection is a plane curve.

        INPUT:

        - ``AP`` -- (default: ``None``) the affine plane to project this curve
          into. This space must be defined over the same base field as this
          curve, and must have dimension two. This space will be constructed if
          not specified.

        OUTPUT: a tuple of

        - a scheme morphism from this curve into an affine plane

        - the plane curve that defines the image of that morphism

        EXAMPLES::

            sage: A.<x,y,z,w> = AffineSpace(QQ, 4)
            sage: C = Curve([x^2 - y*z*w, z^3 - w, w + x*y - 3*z^3], A)
            sage: C.plane_projection()
            (Scheme morphism:
              From: Affine Curve over Rational Field defined by
                    -y*z*w + x^2, z^3 - w, -3*z^3 + x*y + w
              To:   Affine Space of dimension 2 over Rational Field
              Defn: Defined on coordinates by sending (x, y, z, w) to (x, y),
             Affine Plane Curve over Rational Field defined by
             x0^2*x1^7 - 16*x0^4)

        ::

            sage: # needs sage.rings.number_field
            sage: R.<a> = QQ[]
            sage: K.<b> = NumberField(a^2 + 2)
            sage: A.<x,y,z> = AffineSpace(K, 3)
            sage: C = A.curve([x - b, y - 2])
            sage: B.<a,b> = AffineSpace(K, 2)
            sage: proj1 = C.plane_projection(AP=B)
            sage: proj1
            (Scheme morphism:
               From: Affine Curve over Number Field in b
                     with defining polynomial a^2 + 2 defined by x + (-b), y - 2
               To:   Affine Space of dimension 2 over Number Field in b
                     with defining polynomial a^2 + 2
               Defn: Defined on coordinates by sending (x, y, z) to
                     (x, z),
             Affine Plane Curve over Number Field in b
             with defining polynomial a^2 + 2 defined by a + (-b))
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
        for i in range(n - 1):
            for j in range(i + 1, n):
                L = self.projection([i, j], AP)
                if isinstance(L[1], Curve_generic):
                    return L

    def blowup(self, P=None):
        r"""
        Return the blow up of this affine curve at the point ``P``.

        The blow up is described by affine charts. This curve must be irreducible.

        INPUT:

        - ``P`` -- (default: ``None``) a point on this curve at which to blow up;
          if ``None``, then ``P`` is taken to be the origin

        OUTPUT: a tuple of

        - a tuple of curves in affine space of the same dimension as the
          ambient space of this curve, which define the blow up in each affine
          chart.

        - a tuple of tuples such that the j-th element of the i-th tuple is the
          transition map from the i-th affine patch to the j-th affine patch.

        - a tuple consisting of the restrictions of the projection map from the
          blow up back to the original curve, restricted to each affine patch.
          There the i-th element will be the projection from the i-th affine patch.

        EXAMPLES::

            sage: A.<x,y> = AffineSpace(QQ, 2)
            sage: C = Curve([y^2 - x^3], A)
            sage: C.blowup()
            ((Affine Plane Curve over Rational Field defined by s1^2 - x,
              Affine Plane Curve over Rational Field defined by y*s0^3 - 1),
            ([Scheme endomorphism of Affine Plane Curve over Rational Field
               defined by s1^2 - x
                Defn: Defined on coordinates by sending (x, s1) to (x, s1),
              Scheme morphism:
                From: Affine Plane Curve over Rational Field defined by s1^2 - x
                To:   Affine Plane Curve over Rational Field defined by y*s0^3 - 1
                Defn: Defined on coordinates by sending (x, s1) to (x*s1, 1/s1)],
             [Scheme morphism:
                From: Affine Plane Curve over Rational Field defined by y*s0^3 - 1
                To:   Affine Plane Curve over Rational Field defined by s1^2 - x
                Defn: Defined on coordinates by sending (y, s0) to (y*s0, 1/s0),
              Scheme endomorphism of Affine Plane Curve over Rational Field
               defined by y*s0^3 - 1
                Defn: Defined on coordinates by sending (y, s0) to (y, s0)]),
            (Scheme morphism:
               From: Affine Plane Curve over Rational Field defined by s1^2 - x
               To:   Affine Plane Curve over Rational Field defined by -x^3 + y^2
               Defn: Defined on coordinates by sending (x, s1) to (x, x*s1),
             Scheme morphism:
               From: Affine Plane Curve over Rational Field defined by y*s0^3 - 1
               To:   Affine Plane Curve over Rational Field defined by -x^3 + y^2
               Defn: Defined on coordinates by sending (y, s0) to (y*s0, y)))

        ::

            sage: # needs sage.rings.number_field
            sage: K.<a> = QuadraticField(2)
            sage: A.<x,y,z> = AffineSpace(K, 3)
            sage: C = Curve([y^2 - a*x^5, x - z], A)
            sage: B = C.blowup()
            sage: B[0]
            (Affine Curve over Number Field in a with defining polynomial x^2 - 2
              with a = 1.414213562373095? defined by s2 - 1, 2*x^3 + (-a)*s1^2,
             Affine Curve over Number Field in a with defining polynomial x^2 - 2
              with a = 1.414213562373095? defined by s0 - s2, 2*y^3*s2^5 + (-a),
             Affine Curve over Number Field in a with defining polynomial x^2 - 2
              with a = 1.414213562373095? defined by s0 - 1, 2*z^3 + (-a)*s1^2)
            sage: B[1][0][2]
            Scheme morphism:
              From: Affine Curve over Number Field in a
                    with defining polynomial x^2 - 2 with a = 1.414213562373095?
                    defined by s2 - 1, 2*x^3 + (-a)*s1^2
              To:   Affine Curve over Number Field in a
                    with defining polynomial x^2 - 2 with a = 1.414213562373095?
                    defined by s0 - 1, 2*z^3 + (-a)*s1^2
              Defn: Defined on coordinates by sending (x, s1, s2) to
                    (x*s2, 1/s2, s1/s2)
            sage: B[1][2][0]
            Scheme morphism:
              From: Affine Curve over Number Field in a
                    with defining polynomial x^2 - 2 with a = 1.414213562373095?
                    defined by s0 - 1, 2*z^3 + (-a)*s1^2
              To:   Affine Curve over Number Field in a
                    with defining polynomial x^2 - 2 with a = 1.414213562373095?
                    defined by s2 - 1, 2*x^3 + (-a)*s1^2
              Defn: Defined on coordinates by sending (z, s0, s1) to
                    (z*s0, s1/s0, 1/s0)
            sage: B[2]
            (Scheme morphism:
               From: Affine Curve over Number Field in a
                     with defining polynomial x^2 - 2 with a = 1.414213562373095?
                     defined by s2 - 1, 2*x^3 + (-a)*s1^2
               To:   Affine Curve over Number Field in a
                     with defining polynomial x^2 - 2 with a = 1.414213562373095?
                     defined by (-a)*x^5 + y^2, x - z
               Defn: Defined on coordinates by sending (x, s1, s2) to
                     (x, x*s1, x*s2),
             Scheme morphism:
               From: Affine Curve over Number Field in a
                     with defining polynomial x^2 - 2 with a = 1.414213562373095?
                     defined by s0 - s2, 2*y^3*s2^5 + (-a)
               To:   Affine Curve over Number Field in a
                     with defining polynomial x^2 - 2 with a = 1.414213562373095?
                     defined by (-a)*x^5 + y^2, x - z
               Defn: Defined on coordinates by sending (y, s0, s2) to
                     (y*s0, y, y*s2),
             Scheme morphism:
               From: Affine Curve over Number Field in a
                     with defining polynomial x^2 - 2 with a = 1.414213562373095?
                     defined by s0 - 1, 2*z^3 + (-a)*s1^2
               To:   Affine Curve over Number Field in a
                     with defining polynomial x^2 - 2 with a = 1.414213562373095?
                     defined by (-a)*x^5 + y^2, x - z
               Defn: Defined on coordinates by sending (z, s0, s1) to
                     (z*s0, z*s1, z))

        ::

            sage: A.<x,y> = AffineSpace(QQ, 2)
            sage: C = A.curve((y - 3/2)^3 - (x + 2)^5 - (x + 2)^6)
            sage: Q = A([-2,3/2])
            sage: C.blowup(Q)
            ((Affine Plane Curve over Rational Field
               defined by x^3 - s1^3 + 7*x^2 + 16*x + 12,
              Affine Plane Curve over Rational Field
               defined by 8*y^3*s0^6 - 36*y^2*s0^6 + 8*y^2*s0^5
                          + 54*y*s0^6 - 24*y*s0^5 - 27*s0^6 + 18*s0^5 - 8),
             ([Scheme endomorphism of Affine Plane Curve over Rational Field
                defined by x^3 - s1^3 + 7*x^2 + 16*x + 12
                 Defn: Defined on coordinates by sending (x, s1) to (x, s1),
               Scheme morphism:
                 From: Affine Plane Curve over Rational Field
                       defined by x^3 - s1^3 + 7*x^2 + 16*x + 12
                 To:   Affine Plane Curve over Rational Field
                       defined by 8*y^3*s0^6 - 36*y^2*s0^6 + 8*y^2*s0^5
                                  + 54*y*s0^6 - 24*y*s0^5 - 27*s0^6 + 18*s0^5 - 8
                 Defn: Defined on coordinates by sending (x, s1) to
                       (x*s1 + 2*s1 + 3/2, 1/s1)],
              [Scheme morphism:
                 From: Affine Plane Curve over Rational Field
                       defined by 8*y^3*s0^6 - 36*y^2*s0^6 + 8*y^2*s0^5
                                  + 54*y*s0^6 - 24*y*s0^5 - 27*s0^6 + 18*s0^5 - 8
                 To:   Affine Plane Curve over Rational Field
                       defined by x^3 - s1^3 + 7*x^2 + 16*x + 12
                 Defn: Defined on coordinates by sending (y, s0) to
                       (y*s0 - 3/2*s0 - 2, 1/s0),
               Scheme endomorphism of Affine Plane Curve over Rational Field
                defined by 8*y^3*s0^6 - 36*y^2*s0^6 + 8*y^2*s0^5 + 54*y*s0^6
                           - 24*y*s0^5 - 27*s0^6 + 18*s0^5 - 8
                 Defn: Defined on coordinates by sending (y, s0) to (y, s0)]),
             (Scheme morphism:
                From: Affine Plane Curve over Rational Field
                      defined by x^3 - s1^3 + 7*x^2 + 16*x + 12
                To:   Affine Plane Curve over Rational Field
                      defined by -x^6 - 13*x^5 - 70*x^4 - 200*x^3 + y^3
                                 - 320*x^2 - 9/2*y^2 - 272*x + 27/4*y - 795/8
                Defn: Defined on coordinates by sending (x, s1) to
                      (x, x*s1 + 2*s1 + 3/2),
              Scheme morphism:
                From: Affine Plane Curve over Rational Field
                      defined by 8*y^3*s0^6 - 36*y^2*s0^6 + 8*y^2*s0^5
                                 + 54*y*s0^6 - 24*y*s0^5 - 27*s0^6 + 18*s0^5 - 8
                To:   Affine Plane Curve over Rational Field
                      defined by -x^6 - 13*x^5 - 70*x^4 - 200*x^3 + y^3
                                 - 320*x^2 - 9/2*y^2 - 272*x + 27/4*y - 795/8
                Defn: Defined on coordinates by sending (y, s0) to
                      (y*s0 - 3/2*s0 - 2, y)))

        ::

            sage: A.<x,y,z,w> = AffineSpace(QQ, 4)
            sage: C = A.curve([((x + 1)^2 + y^2)^3 - 4*(x + 1)^2*y^2, y - z, w - 4])
            sage: Q = C([-1,0,0,4])
            sage: B = C.blowup(Q)
            sage: B[0]
            (Affine Curve over Rational Field defined by s3, s1 - s2,
              x^2*s2^6 + 2*x*s2^6 + 3*x^2*s2^4 + s2^6 + 6*x*s2^4
              + 3*x^2*s2^2 + 3*s2^4 + 6*x*s2^2 + x^2 - s2^2 + 2*x + 1,
             Affine Curve over Rational Field defined by s3, s2 - 1,
              y^2*s0^6 + 3*y^2*s0^4 + 3*y^2*s0^2 + y^2 - 4*s0^2,
             Affine Curve over Rational Field defined by s3, s1 - 1,
              z^2*s0^6 + 3*z^2*s0^4 + 3*z^2*s0^2 + z^2 - 4*s0^2,
             Closed subscheme of Affine Space of dimension 4 over Rational Field
              defined by: 1)
            sage: Q = A([6,2,3,1])
            sage: B = C.blowup(Q)
            Traceback (most recent call last):
            ...
            TypeError: (=(6, 2, 3, 1)) must be a point on this curve

        ::

            sage: # needs sage.rings.number_field
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
        if self.base_ring() not in Fields():
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
        return (tuple(patches), tuple(t_maps), tuple(proj_maps))

    def resolution_of_singularities(self, extend=False):
        r"""
        Return a nonsingular model for this affine curve created by blowing up
        its singular points.

        The nonsingular model is given as a collection of affine patches that
        cover it. If ``extend`` is ``False`` and if the base field is a number
        field, or if the base field is a finite field, the model returned may
        have singularities with coordinates not contained in the base field. An
        error is returned if this curve is already nonsingular, or if it has no
        singular points over its base field. This curve must be irreducible,
        and must be defined over a number field or finite field.

        INPUT:

        - ``extend`` -- boolean (default: ``False``); specifies whether to
          extend the base field when necessary to find all singular points when
          this curve is defined over a number field. If ``extend`` is
          ``False``, then only singularities with coordinates in the base field
          of this curve will be resolved. However, setting ``extend`` to
          ``True`` will slow down computations.

        OUTPUT: a tuple of

        - a tuple of curves in affine space of the same dimension as the
          ambient space of this curve, which represent affine patches of the
          resolution of singularities.

        - a tuple of tuples such that the j-th element of the i-th tuple is the
          transition map from the i-th patch to the j-th patch.

        - a tuple consisting of birational maps from the patches back to the
          original curve that were created by composing the projection maps
          generated from the blow up computations. There the i-th element will
          be a map from the i-th patch.

        EXAMPLES::

            sage: A.<x,y> = AffineSpace(QQ, 2)
            sage: C = Curve([y^2 - x^3], A)
            sage: C.resolution_of_singularities()
            ((Affine Plane Curve over Rational Field defined by s1^2 - x,
              Affine Plane Curve over Rational Field defined by y*s0^3 - 1),
             ((Scheme endomorphism of Affine Plane Curve over Rational Field
                defined by s1^2 - x
                 Defn: Defined on coordinates by sending (x, s1) to (x, s1),
               Scheme morphism:
                 From: Affine Plane Curve over Rational Field defined by s1^2 - x
                 To:   Affine Plane Curve over Rational Field defined by y*s0^3 - 1
                 Defn: Defined on coordinates by sending (x, s1) to (x*s1, 1/s1)),
              (Scheme morphism:
                 From: Affine Plane Curve over Rational Field defined by y*s0^3 - 1
                 To:   Affine Plane Curve over Rational Field defined by s1^2 - x
                 Defn: Defined on coordinates by sending (y, s0) to (y*s0, 1/s0),
               Scheme endomorphism of Affine Plane Curve over Rational Field
                defined by y*s0^3 - 1
                 Defn: Defined on coordinates by sending (y, s0) to (y, s0))),
             (Scheme morphism:
                From: Affine Plane Curve over Rational Field defined by s1^2 - x
                To:   Affine Plane Curve over Rational Field defined by -x^3 + y^2
                Defn: Defined on coordinates by sending (x, s1) to (x, x*s1),
              Scheme morphism:
                From: Affine Plane Curve over Rational Field defined by y*s0^3 - 1
                To:   Affine Plane Curve over Rational Field defined by -x^3 + y^2
                Defn: Defined on coordinates by sending (y, s0) to (y*s0, y)))

        ::

            sage: # needs sage.rings.number_field
            sage: set_verbose(-1)
            sage: K.<a> = QuadraticField(3)
            sage: A.<x,y> = AffineSpace(K, 2)
            sage: C = A.curve(x^4 + 2*x^2 + a*y^3 + 1)
            sage: C.resolution_of_singularities(extend=True)[0]         # long time (2 s)
            (Affine Plane Curve over Number Field in a0
              with defining polynomial y^4 - 4*y^2 + 16
              defined by 24*x^2*ss1^3 + 24*ss1^3 + (a0^3 - 8*a0),
             Affine Plane Curve over Number Field in a0
              with defining polynomial y^4 - 4*y^2 + 16
              defined by 24*s1^2*ss0 + (a0^3 - 8*a0)*ss0^2 + (-6*a0^3)*s1,
             Affine Plane Curve over Number Field in a0
              with defining polynomial y^4 - 4*y^2 + 16
              defined by 8*y^2*s0^4 + (4*a0^3)*y*s0^3 - 32*s0^2 + (a0^3 - 8*a0)*y)

        ::

            sage: A.<x,y,z> = AffineSpace(GF(5), 3)
            sage: C = Curve([y - x^3, (z - 2)^2 - y^3 - x^3], A)
            sage: R = C.resolution_of_singularities()
            sage: R[0]
            (Affine Curve over Finite Field of size 5
              defined by x^2 - s1, s1^4 - x*s2^2 + s1, x*s1^3 - s2^2 + x,
             Affine Curve over Finite Field of size 5
              defined by y*s2^2 - y^2 - 1, s2^4 - s0^3 - y^2 - 2, y*s0^3 - s2^2 + y,
             Affine Curve over Finite Field of size 5
              defined by s0^3*s1 + z*s1^3 + s1^4 - 2*s1^3 - 1,
                         z*s0^3 + z*s1^3 - 2*s0^3 - 2*s1^3 - 1,
                         z^2*s1^3 + z*s1^3 - s1^3 - z + s1 + 2)

        ::

            sage: A.<x,y,z,w> = AffineSpace(QQ, 4)
            sage: C = A.curve([((x - 2)^2 + y^2)^2 - (x - 2)^2 - y^2 + (x - 2)^3,
            ....:              z - y - 7, w - 4])
            sage: B = C.resolution_of_singularities()
            sage: B[0]
            (Affine Curve over Rational Field defined by s3, s1 - s2,
              x^2*s2^4 - 4*x*s2^4 + 2*x^2*s2^2 + 4*s2^4 - 8*x*s2^2
              + x^2 + 7*s2^2 - 3*x + 1,
             Affine Curve over Rational Field defined by s3, s2 - 1,
              y^2*s0^4 + 2*y^2*s0^2 + y*s0^3 + y^2 - s0^2 - 1,
             Affine Curve over Rational Field defined by s3, s1 - 1,
              z^2*s0^4 - 14*z*s0^4 + 2*z^2*s0^2 + z*s0^3 + 49*s0^4
              - 28*z*s0^2 - 7*s0^3 + z^2 + 97*s0^2 - 14*z + 48,
             Closed subscheme of Affine Space of dimension 4 over Rational Field
              defined by: 1)

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
            from sage.rings.qqbar import QQbar
            pts = self.change_ring(F.embeddings(QQbar)[0]).rational_points()
            L = [t for pt in pts for t in pt]
            K = number_field_elements_from_algebraics(L)[0]
            if isinstance(K, RationalField):
                return F.embeddings(F)[0]
            else:
                if isinstance(F, RationalField):
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
        if not res[0][3]:
            if C.is_smooth():
                raise TypeError("this curve is already nonsingular")
            else:
                raise TypeError("this curve has no singular points over "
                                "its base field. If working over "
                                "a number field use extend=True")
        not_resolved = True
        t = 0
        # loop through the patches and blow up each until no patch has singular points
        while not_resolved:
            [BC, t_maps, pi, pts] = [res[t][0], res[t][1], res[t][2], res[t][3]]
            # check if there are any singular points in this patch
            if not pts:
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
                    if n_pts:
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
                coords = [(BC.ambient_space().gens()[j] - pts[0][j]) / (BC.ambient_space().gens()[i] - pts[0][i])
                          for j in range(n)]
                coords.pop(i)
                coords.insert(0, BC.ambient_space().gens()[i])
                H = Hom(BC, B[0][i])
                t_pi_inv = H(coords)
                # compose the current transition maps from the original curve to the other patches
                # with the projection map
                L = list(t_maps)
                for j in range(len(t_maps)):
                    L[j] = L[j] * t_pi
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
        return (tuple(patches), tuple(t_maps), tuple(p_maps))

    def tangent_line(self, p):
        """
        Return the tangent line at the point ``p``.

        INPUT:

        - ``p`` -- a rational point of the curve

        EXAMPLES::

            sage: A3.<x,y,z> = AffineSpace(3, QQ)
            sage: C = Curve([x + y + z, x^2 - y^2*z^2 + z^3])
            sage: p = C(0,0,0)
            sage: C.tangent_line(p)
            Traceback (most recent call last):
            ...
            ValueError: the curve is not smooth at (0, 0, 0)
            sage: p = C(1,0,-1)
            sage: C.tangent_line(p)
            Affine Curve over Rational Field defined by x + y + z, 2*x + 3*z + 1

        We check that the tangent line at ``p`` is the tangent space at ``p``,
        translated to ``p``. ::

            sage: Tp = C.tangent_space(p)
            sage: Tp
            Closed subscheme of Affine Space of dimension 3 over Rational Field
             defined by: x + y + z, 2*x + 3*z
            sage: phi = A3.translation(A3.origin(), p)
            sage: T = phi * Tp.embedding_morphism()
            sage: T.image()
            Closed subscheme of Affine Space of dimension 3 over Rational Field
             defined by: -2*y + z + 1, x + y + z
            sage: _ == C.tangent_line(p)
            True
        """
        A = self.ambient_space()
        R = A.coordinate_ring()
        gens = R.gens()

        Tp = self.tangent_space(p)

        if Tp.dimension() > 1:
            raise ValueError("the curve is not smooth at {}".format(p))

        from sage.schemes.curves.constructor import Curve

        # translate to p
        I0 = []
        for poly in Tp.defining_polynomials():
            I0.append(poly.subs({x: x - c for x, c in zip(gens, p)}))

        return Curve(I0, A)


class AffinePlaneCurve_field(AffinePlaneCurve, AffineCurve_field):
    """
    Affine plane curves over fields.
    """
    _point = AffinePlaneCurvePoint_field

    def has_vertical_asymptote(self):
        """
        Check if the curve is not a line and has vertical asymptotes.

        EXAMPLES::

            sage: A2.<x,y> = AffineSpace(2, QQ)
            sage: Curve(x).has_vertical_asymptote()
            False
            sage: Curve(y^2 * x + x + y).has_vertical_asymptote()
            True
        """
        A = self.ambient_space()
        R = A.coordinate_ring()
        x, y = R.gens()
        f = self.defining_polynomial().radical()
        dy = f.degree(y)
        dxy = f.coefficient({y: dy}).degree()
        return dxy > 0 and f.degree() > 1

    def is_vertical_line(self):
        """
        Check if the curve is a vertical line.

        EXAMPLES::

            sage: A2.<x, y> = AffineSpace(2, QQ)
            sage: Curve(x - 1).is_vertical_line()
            True
            sage: Curve(x - y).is_vertical_line()
            False
            sage: Curve(y^2 * x + x + y).is_vertical_line()
            False
        """
        A = self.ambient_space()
        R = A.coordinate_ring()
        x, y = R.gens()
        f = self.defining_polynomial().radical()
        return f.degree(y) == 0 and f.degree() == 1

    @cached_method
    def fundamental_group(self, simplified=True, puiseux=True):
        r"""
        Return a presentation of the fundamental group of the complement
        of ``self``.

        INPUT:

        - ``simplified`` -- boolean (default: ``True``); to simplify the presentation

        - ``puiseux`` -- boolean (default: ``True``); to decide if the
          presentation is constructed in the classical way or using Puiseux
          shortcut

        OUTPUT:

        A presentation with generators `x_1, \dots, x_d` and relations. If ``puiseux``
        is ``False`` the relations are `(x_j\cdot \tau)\cdot x_j^{-1}` for `1\leq j<d`
        and `tau` a braid in the braid monodromy; finally the presentation
        is simplified. If ``puiseux`` is ``True``, each
        `tau` is decomposed as `\alpha^{-1}\cdot\beta\cdot\alpha`, where `\beta` is
        a positive braid; the relations are `((x_j\cdot \beta)\cdot x_j^{-1})\cdot \alpha`
        where `j` is an integer of the ``Tietze`` word of `\beta`. This presentation
        is not simplified by default since it represents the homotopy type of
        the complement of the curve.

        .. NOTE::

            The curve must be defined over the rationals or a number field
            with an embedding over `\QQbar`. This functionality requires
            the ``sirocco`` package to be installed.

        EXAMPLES::

            sage: # needs sirocco
            sage: A.<x,y> = AffineSpace(QQ, 2)
            sage: C = A.curve(y^2 - x^3 - x^2)
            sage: C.fundamental_group(puiseux=False)
            Finitely presented group < x0 |  >
            sage: bm = C.braid_monodromy()
            sage: g = C.fundamental_group(simplified=False)
            sage: g.sorted_presentation()
            Finitely presented group < x0, x1 | x1^-1*x0^-1*x1*x0, x1^-1*x0 >
            sage: g.simplified()
            Finitely presented group < x0 |  >

        In the case of number fields, they need to have an embedding
        to the algebraic field::

            sage: # needs sage.rings.number_field
            sage: x = polygen(ZZ)
            sage: a = QQ[x](x^2 + 5).roots(QQbar)[0][0]
            sage: F = NumberField(a.minpoly(), 'a', embedding=a)
            sage: F.inject_variables()
            Defining a
            sage: A.<x,y> = AffineSpace(F, 2)
            sage: C = A.curve(y^2 - a*x^3 - x^2)
            sage: C.fundamental_group()                     # needs sirocco
            Finitely presented group < x0 |  >
            sage: C = A.curve(x * (x - 1))
            sage: C.fundamental_group()                     # needs sirocco
            Finitely presented group < x0, x1 |  >
        """
        from sage.schemes.curves.zariski_vankampen import fundamental_group_from_braid_mon
        bm = self.braid_monodromy()
        if not bm:
            f = self.defining_polynomial()
            x, y = f.parent().gens()
            d0 = f.degree(y)
            f0 = f.coefficient({y: d0})
            d = d0 + f0.degree(x)
        else:
            d = bm[0].parent().strands()
        G = fundamental_group_from_braid_mon(bm, degree=d,
                                             simplified=simplified,
                                             puiseux=puiseux)
        if simplified:
            G = G.simplified()
        return G

    @cached_method
    def braid_monodromy(self):
        r"""
        Compute the braid monodromy of a projection of the curve.

        OUTPUT:

        A list of braids. The braids correspond to paths based in the same point;
        each of this paths is the conjugated of a loop around one of the points
        in the discriminant of the projection of ``self``.

        .. NOTE::

            The projection over the `x` axis is used if there are no vertical asymptotes.
            Otherwise, a linear change of variables is done to fall into the previous case.

        .. NOTE::

            This functionality requires the ``sirocco`` package to be installed.

        EXAMPLES::

            sage: A.<x,y> = AffineSpace(QQ, 2)
            sage: C = A.curve((x^2-y^3)*(x+3*y-5))
            sage: C.braid_monodromy()                                   # needs sirocco
            [s1*s0*(s1*s2)^2*s0*s2^2*s0^-1*(s2^-1*s1^-1)^2*s0^-1*s1^-1,
             s1*s0*(s1*s2)^2*(s0*s2^-1*s1*s2*s1*s2^-1)^2*(s2^-1*s1^-1)^2*s0^-1*s1^-1,
             s1*s0*(s1*s2)^2*s2*s1^-1*s2^-1*s1^-1*s0^-1*s1^-1,
             s1*s0*s2*s0^-1*s2*s1^-1]
            sage: T.<t> = QQ[]
            sage: K.<a> = NumberField(t^3 + 2, 'a')
            sage: A.<x, y> = AffineSpace(K, 2)
            sage: Curve(y^2 + a * x).braid_monodromy()
            Traceback (most recent call last):
            ...
            NotImplementedError: the base field must have an embedding to the algebraic field
        """
        from sage.schemes.curves.zariski_vankampen import braid_monodromy
        F = self.base_ring()
        from sage.rings.qqbar import QQbar
        if QQbar.coerce_map_from(F) is None:
            raise NotImplementedError("the base field must have an embedding"
                                      " to the algebraic field")
        f = self.defining_polynomial()
        return braid_monodromy(f)[0]

    def riemann_surface(self, **kwargs):
        r"""
        Return the complex Riemann surface determined by this curve.

        OUTPUT: a :class:`~sage.schemes.riemann_surfaces.riemann_surface.RiemannSurface` object

        EXAMPLES::

            sage: R.<x,y> = QQ[]
            sage: C = Curve(x^3 + 3*y^3 + 5)
            sage: C.riemann_surface()
            Riemann surface defined by polynomial f = x^3 + 3*y^3 + 5 = 0,
             with 53 bits of precision
        """
        from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
        S = RiemannSurface(self.defining_polynomial(), **kwargs)
        S._curve = self
        return S


class AffinePlaneCurve_finite_field(AffinePlaneCurve_field):
    """
    Affine plane curves over finite fields.
    """
    _point = AffinePlaneCurvePoint_finite_field

    # CHECK WHAT ASSUMPTIONS ARE MADE REGARDING AFFINE VS. PROJECTIVE MODELS!!!
    # THIS IS VERY DIRTY STILL -- NO DATASTRUCTURES FOR DIVISORS.
    def riemann_roch_basis(self, D):
        r"""
        Return a basis of the Riemann-Roch space of the divisor ``D``.

        This interfaces with Singular's Brill-Noether command.

        This curve is assumed to be a plane curve defined by a polynomial
        equation `f(x,y) = 0` over a prime finite field `F = GF(p)` in 2
        variables `x,y` representing a curve `X: f(x,y) = 0` having `n`
        `F`-rational points (see the Sage function ``places_on_curve``)

        INPUT:

        - ``D`` -- an `n`-tuple of integers `(d_1, ..., d_n)` representing the
          divisor `Div = d_1P_1 + \dots + d_nP_n`, where `X(F) = \{P_1, \dots,
          P_n\}`.  The ordering is that dictated by ``places_on_curve``.

        OUTPUT: a basis of `L(Div)`.

        EXAMPLES::

            sage: R = PolynomialRing(GF(5), 2, names=["x","y"])
            sage: x, y = R.gens()
            sage: f = y^2 - x^9 - x
            sage: C = Curve(f)
            sage: D = [6,0,0,0,0,0]
            sage: C.riemann_roch_basis(D)
            [1, (-x*z^5 + y^2*z^4)/x^6, (-x*z^6 + y^2*z^5)/x^7, (-x*z^7 + y^2*z^6)/x^8]
        """
        F = self.base_ring()
        if not F.is_prime_field():
            raise TypeError("only works for curves over prime finite fields")

        p = F.characteristic()
        f = self.defining_polynomial()
        gens = f.parent().gens()

        G = singular(','.join(str(x) for x in D), type='intvec')

        singular.lib('brnoeth')
        singular.ring(p, gens, 'lp')

        X = singular(f).Adj_div()
        P = singular.NSplaces(1, X)
        T = P[1][2]
        T.set_ring()  # necessary

        return [g[1].sage() / g[2].sage() for g in G.BrillNoether(P)]

    def rational_points(self, algorithm='enum'):
        r"""
        Return sorted list of all rational points on this curve.

        INPUT:

        - ``algorithm`` -- possible choices:

           +  ``'enum'`` -- use *very* naive point enumeration to find all
              rational points on this curve over a finite field.

           +  ``'bn'`` -- via Singular's Brill-Noether package.

           +  ``'all'`` -- use all implemented algorithms and verify that they
              give the same answer, then return it

        .. NOTE::

           The Brill-Noether package does not always work. When it fails, a
           :exc:`RuntimeError` exception is raised.

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

            sage: # needs sage.rings.finite_rings
            sage: A.<x,y> = AffineSpace(2, GF(9,'a'))
            sage: C = Curve(x^2 + y^2 - 1); C
            Affine Plane Curve over Finite Field in a of size 3^2
             defined by x^2 + y^2 - 1
            sage: C.rational_points()
            [(0, 1), (0, 2), (1, 0), (2, 0), (a + 1, a + 1),
             (a + 1, 2*a + 2), (2*a + 2, a + 1), (2*a + 2, 2*a + 2)]
        """
        if algorithm == "enum":
            f = self.defining_polynomial()
            K = f.parent().base_ring()
            return sorted(self((x, y)) for x in K for y in K if f(x, y) == 0)

        F = self.base_ring()
        if not F.is_prime_field():
            raise TypeError("other algorithms only work for curves over prime finite fields")

        if algorithm == "bn":
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
            S_enum = self.rational_points(algorithm='enum')
            S_bn = self.rational_points(algorithm='bn')
            if S_enum != S_bn:
                raise RuntimeError("Bug in rational_points -- different algorithms give different answers for curve %s!" % self)
            return S_enum
        else:
            raise ValueError("No algorithm '%s' known" % algorithm)


class IntegralAffineCurve(AffineCurve_field):
    """
    Base class for integral affine curves.
    """
    _point = IntegralAffineCurvePoint
    _closed_point = IntegralAffineCurveClosedPoint

    def function_field(self):
        """
        Return the function field of the curve.

        EXAMPLES::

            sage: A.<x,y> = AffineSpace(QQ, 2)
            sage: C = Curve(x^3 - y^2 - x^4 - y^4)
            sage: C.function_field()
            Function field in y defined by y^4 + y^2 + x^4 - x^3

        ::

            sage: # needs sage.rings.finite_rings
            sage: A.<x,y> = AffineSpace(GF(8), 2)
            sage: C = Curve(x^5 + y^5 + x*y + 1)
            sage: C.function_field()
            Function field in y defined by y^5 + x*y + x^5 + 1
        """
        return self._function_field

    @lazy_attribute
    def _genus(self):
        """
        The geometric genus of the curve.

        EXAMPLES::

            sage: A.<x,y> = AffineSpace(GF(2), 2)
            sage: C = Curve(x^5 + y^5 + x*y + 1)
            sage: C.genus()   # indirect doctest
            1

        TESTS::

            sage: # needs sage.rings.number_field
            sage: R.<T> = QQ[]
            sage: N.<a> = NumberField(T^2 + 1)
            sage: A2.<x,y> = AffineSpace(N, 2)
            sage: C = Curve(y^2 - x^3 + x, A2)
            sage: C.genus()
            1
        """
        # Singular's genus command is usually much faster than the genus method
        # of function fields in Sage. But unfortunately Singular's genus
        # command does not work over extension fields.
        if self.base_ring().degree() > 1:
            return self._function_field.genus()

        # call Singular's genus command
        return self.defining_ideal().genus()

    def __call__(self, *args):
        """
        Return a rational point, a pointset or a function depending on ``args``.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: A.<x,y> = AffineSpace(GF(8), 2)
            sage: C = Curve(x^5 + y^5 + x*y + 1)
            sage: C(1,1)
            (1, 1)
            sage: C(x/y)
            (x/(x^5 + 1))*y^4 + x^2/(x^5 + 1)
            sage: C(GF(8^2))
            Set of rational points of Closed subscheme of Affine Space of dimension 2
            over Finite Field in z6 of size 2^6 defined by: x^5 + y^5 + x*y + 1

        ::

            sage: A.<x,y,z> = AffineSpace(GF(11), 3)
            sage: C = Curve([x*z - y^2, y - z^2, x - y*z], A)
            sage: C([0,0,0])
            (0, 0, 0)
            sage: C(y)
            z^2
            sage: C(A.coordinate_ring()(y))
            z^2
        """
        try:
            return super().__call__(*args)
        except TypeError as e:
            try:
                return self.function(*args)
            except AttributeError:
                raise e

    def function(self, f):
        """
        Return the function field element coerced from ``f``.

        INPUT:

        - ``f`` -- an element of the fraction field of the coordinate ring of
          the ambient space or the coordinate ring of the curve

        OUTPUT: an element of the function field of this curve

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: A.<x,y> = AffineSpace(GF(8), 2)
            sage: C = Curve(x^5 + y^5 + x*y + 1)
            sage: f = C.function(x/y)
            sage: f
            (x/(x^5 + 1))*y^4 + x^2/(x^5 + 1)
            sage: df = f.differential(); df
            ((1/(x^10 + 1))*y^4 + x^6/(x^10 + 1)) d(x)
            sage: df.divisor()                  # long time
            2*Place (1/x, 1/x^4*y^4 + 1/x^3*y^3 + 1/x^2*y^2 + 1/x*y + 1)
             + 2*Place (1/x, 1/x*y + 1)
             - 2*Place (x + 1, y)
             - 2*Place (x^4 + x^3 + x^2 + x + 1, y)
        """
        R = self.ambient_space().coordinate_ring()
        if f not in R and f.parent() is self.coordinate_ring():
            f = f.lift()

        phi = self._map_to_function_field
        num = R(f.numerator())
        den = R(f.denominator())
        return phi(num) / phi(den)

    def coordinate_functions(self):
        """
        Return the coordinate functions.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: A.<x,y> = AffineSpace(GF(8), 2)
            sage: C = Curve(x^5 + y^5 + x*y + 1)
            sage: x, y = C.coordinate_functions()
            sage: x^5 + y^5 + x*y + 1
            0
        """
        return self._coordinate_functions

    def pull_from_function_field(self, f):
        """
        Return the fraction corresponding to ``f``.

        INPUT:

        - ``f`` -- an element of the function field

        OUTPUT:

        A fraction of polynomials in the coordinate ring of the ambient space
        of the curve.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: A.<x,y> = AffineSpace(GF(8), 2)
            sage: C = Curve(x^5 + y^5 + x*y + 1)
            sage: F = C.function_field()
            sage: C.pull_from_function_field(F.gen())
            y
            sage: C.pull_from_function_field(F.one())
            1
            sage: C.pull_from_function_field(F.zero())
            0
            sage: f1 = F.gen()
            sage: f2 = F.base_ring().gen()
            sage: C.function(C.pull_from_function_field(f1)) == f1
            True
            sage: C.function(C.pull_from_function_field(f2)) == f2
            True
        """
        return self._map_from_function_field(f)

    @lazy_attribute
    def _nonsingular_model(self):
        """
        Return the data of a nonsingular model of the curve.

        The data consists of an abstract function field `M`, a map from the
        fraction field of the coordinate ring `R` of the ambient space of the
        curve to the function field, and the inverse map.

        The coordinate ring of the curve is the quotient of `R` by the kernel
        of the map restricted to `R`.

        TESTS::

            sage: A.<x,y,z> = AffineSpace(GF(11), 3)
            sage: C = Curve([x*z - y^2, y - z^2, x - y*z], A)
            sage: C._nonsingular_model
            (Function field in z defined by z^3 + 10*x,
             Ring morphism:
               From: Fraction Field of Multivariate Polynomial Ring in x, y, z
                     over Finite Field of size 11
               To:   Function field in z defined by z^3 + 10*x
               Defn: x |--> x
                     y |--> z^2
                     z |--> z,
             Ring morphism:
               From: Function field in z defined by z^3 + 10*x
               To:   Fraction Field of Multivariate Polynomial Ring in x, y, z
                     over Finite Field of size 11)
        """
        from sage.structure.sequence import Sequence
        from sage.rings.fraction_field import FractionField
        from sage.rings.function_field.constructor import FunctionField
        from sage.rings.function_field.maps import FunctionFieldRingMorphism

        k = self.base_ring()
        I = self.defining_ideal()

        # invlex is the lex order with x < y < z for R = k[x,y,z] for instance
        R = I.parent().ring().change_ring(order='invlex')
        I0 = I.change_ring(R)
        n = R.ngens()

        names = R.variable_names()

        gbasis = I0.groebner_basis()

        if not I0.is_prime():
            raise TypeError("the curve is not integral")

        # Suppose the generators of the defining ideal I of the curve is
        #
        #       -y^2 + x*z, -z^2 + y, -y*z + x.
        #
        # Then the Groebner basis of the ideal with respect to the elimination
        # order invlex is
        #
        #        f0 = z^2 - y,
        #        f1 = y*z - x,
        #        f2 = x*z - y^2,
        #        f3 = y^3 - x^2.
        #
        # Now the task is to find f that has minimal degree as a polynomial
        # in the i-th variable. The result is
        #
        #        f0 z^2
        #        f1 y*z
        #        f2 x*z                       o
        #        f3 y^3             o
        #        ------------------------------
        #                k[x]  k[x,y]  k[x,y,z]
        #
        # Hence x is an independent variable; f3 is the syzygy for y; f2 is the
        # syzygy for z. Now x is the generator of a rational function field F0;
        # y is the generator of the extension F1 of F0 by f3; z is the
        # generator of the extension F2 of F1 by f2.

        basis = list(gbasis)
        syzygy = {}
        for i in range(n):
            S = k[R._first_ngens(i + 1)]
            while basis:
                f = basis.pop()
                if f in S:
                    if i not in syzygy and f:
                        syzygy[i] = f
                else:
                    basis.append(f)
                    break

        # sanity check
        indeps = [i for i in range(n) if i not in syzygy]
        if len(indeps) != 1:
            raise TypeError("not a curve")
        indep = indeps[0]

        F = FunctionField(k, names[indep])
        coords = {indep: F.gen()}

        for i in range(n):
            if i == indep:
                continue
            P = PolynomialRing(F, 'T')
            f = P([R(c).subs(coords) for c in syzygy[i].polynomial(R.gen(i))])
            F = F.extension(f, names[i])
            coords[i] = F.gen()

        proper_extension = F.base_field() is not F

        if proper_extension:
            N, from_N, to_N = F.simple_model()
            M, from_M, to_M = N.separable_model()
            coordinate_functions = tuple([to_M(to_N(F(coords[i]))) for i in range(n)])
        else:
            M = F  # is rational function field
            coordinate_functions = tuple([coords[i] for i in range(n)])

        # map to M

        FR = FractionField(I.ring())
        map_to_function_field = hom(FR, M, coordinate_functions)

        # map from M

        def convert(f, i):
            if i == indep:
                i = i - 1
            if i < 0:
                return f._x  # fraction representing rational function field element
            fx = f._x  # polynomial representing function field element
            if not fx:
                fxlist = [fx.base_ring().zero()]
            else:
                fxlist = fx.list()
            coeffs = Sequence(convert(c, i - 1) for c in fxlist)
            B = coeffs.universe()
            S = B[names[i]]
            return S(coeffs)

        z = M.gen()

        if proper_extension:
            Z = FR(convert(from_N(from_M(z)), n - 1))

            def evaluate(f):
                coeffs = f._x.list()
                v = 0
                while coeffs:
                    v = v * Z + coeffs.pop()._x
                return FR(v)
        else:
            def evaluate(f):
                return FR(f._x)

        map_from_function_field = FunctionFieldRingMorphism(Hom(M, FR), evaluate)

        # sanity check
        assert all(map_to_function_field(f).is_zero() for f in I.gens())
        assert map_to_function_field(map_from_function_field(z)) == z

        return M, map_to_function_field, map_from_function_field

    @lazy_attribute
    def _function_field(self):
        """
        Return the abstract function field of the curve.

        TESTS::

            sage: A.<x,y,z> = AffineSpace(GF(11), 3)
            sage: C = Curve([x*z - y^2, y - z^2, x - y*z], A)
            sage: C._function_field
            Function field in z defined by z^3 + 10*x
        """
        return self._nonsingular_model[0]

    @lazy_attribute
    def _map_to_function_field(self):
        """
        Return the map to the function field of the curve.

        TESTS::

            sage: A.<x,y,z> = AffineSpace(GF(11), 3)
            sage: C = Curve([x*z - y^2, y - z^2, x - y*z], A)
            sage: C._map_to_function_field
            Ring morphism:
              From: Fraction Field of Multivariate Polynomial Ring in x, y, z
                    over Finite Field of size 11
              To:   Function field in z defined by z^3 + 10*x
              Defn: x |--> x
                    y |--> z^2
                    z |--> z
        """
        return self._nonsingular_model[1]

    @lazy_attribute
    def _coordinate_functions(self):
        """
        Return the coordinate functions of the curve.

        TESTS::

            sage: A.<x,y,z> = AffineSpace(GF(11), 3)
            sage: C = Curve([x*z - y^2, y - z^2, x - y*z], A)
            sage: C._coordinate_functions
            [x, z^2, z]
        """
        return self._nonsingular_model[1].im_gens()

    @lazy_attribute
    def _map_from_function_field(self):
        """
        Return the map from the function field of the curve.

        TESTS::

            sage: A.<x,y,z> = AffineSpace(GF(11), 3)
            sage: C = Curve([x*z - y^2, y - z^2, x - y*z], A)
            sage: C._map_from_function_field
            Ring morphism:
              From: Function field in z defined by z^3 + 10*x
              To:   Fraction Field of Multivariate Polynomial Ring in x, y, z
                    over Finite Field of size 11
        """
        return self._nonsingular_model[2]

    @lazy_attribute
    def _singularities(self):
        """
        Return a list of the pairs of a singular closed point and the places
        above it.

        TESTS::

            sage: # needs sage.rings.finite_rings
            sage: A.<x,y> = AffineSpace(GF(7^2), 2)
            sage: C = Curve(x^2 - x^4 - y^4)
            sage: C._singularities              # long time
            [(Point (x, y),
              [Place (x, 1/x*y^3 + 1/x*y^2 + 1), Place (x, 1/x*y^3 + 1/x*y^2 + 6)])]
        """
        to_F = self._map_to_function_field
        sing = self.singular_subscheme()

        funcs = []
        for p in sing.defining_polynomials():
            f = to_F(p)
            if not f.is_zero():
                funcs.append(f)

        if funcs:
            f = funcs.pop()
            places = f.zeros()
            for f in funcs:
                places = [p for p in places if f.valuation(p) > 0]
        else:
            places = []

        points = []
        for place in places:
            p = self.place_to_closed_point(place)

            for q, places in points:
                if p == q:
                    places.append(place)
                    break
            else:  # new singularity
                points.append((p, [place]))

        return points

    def singular_closed_points(self):
        """
        Return the singular closed points of the curve.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: A.<x,y> = AffineSpace(GF(7^2), 2)
            sage: C = Curve(x^2 - x^4 - y^4)
            sage: C.singular_closed_points()
            [Point (x, y)]

        ::

            sage: A.<x,y,z> = AffineSpace(GF(11), 3)
            sage: C = Curve([x*z - y^2, y - z^2, x - y*z], A)
            sage: C.singular_closed_points()
            []
        """
        return [p for p, _ in self._singularities]

    @cached_method
    def place_to_closed_point(self, place):
        """
        Return the closed point on the place.

        INPUT:

        - ``place`` -- a place of the function field of the curve

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: A.<x,y> = AffineSpace(GF(4), 2)
            sage: C = Curve(x^5 + y^5 + x*y + 1)
            sage: F = C.function_field()
            sage: pls = F.places(1)
            sage: C.place_to_closed_point(pls[-1])
            Point (x + 1, y + 1)
            sage: C.place_to_closed_point(pls[-2])
            Point (x + 1, y + 1)
        """
        F = self.function_field()

        A = self.ambient_space()
        R = A.coordinate_ring().change_ring(order='degrevlex')

        coords = self._coordinate_functions

        if any(f.valuation(place) < 0 for f in coords):
            raise ValueError("the place is at infinity")

        k, from_k, to_k = place.residue_field()
        V, from_V, to_V = k.vector_space(F.constant_base_field(), map=True)

        # implement an FGLM-like algorithm
        e = [0 for i in range(R.ngens())]
        basis = [R.one()]
        basis_vecs = [to_V(k.one())]  # represent as a vector

        gens = []
        gens_lts = []
        terminate = False
        while True:  # check FGLM termination condition
            # compute next exponent in degree reverse lexicographical order
            j = R.ngens() - 1
            while j > 0 and not e[j]:
                j -= 1

            if not j:  # j is zero
                if terminate:
                    break
                terminate = True
                d = e[0]
                e[0] = 0
                e[-1] = d + 1
            else:
                e[j] -= 1
                e[j - 1] += 1

            m = R.monomial(*e)
            if any(g.divides(m) for g in gens_lts):
                continue

            prod = 1
            for i in range(R.ngens()):
                prod *= coords[i]**e[i]
            vec = to_V(to_k(prod))  # represent as a vector
            mat = matrix(basis_vecs)
            try:
                s = mat.solve_left(vec)
            except ValueError:  # no solution
                basis.append(m)
                basis_vecs.append(vec)
                terminate = False
                continue

            gens.append(m - sum([s[i] * basis[i] for i in range(len(basis))]))
            gens_lts.append(m)

        prime = R.ideal(gens).groebner_basis().ideal()

        return self._closed_point(self, prime, len(basis))

    def places_at_infinity(self):
        """
        Return the places of the curve at infinity.

        EXAMPLES::

            sage: A.<x,y> = AffineSpace(QQ, 2)
            sage: C = Curve(x^3 - y^2 - x^4 - y^4)
            sage: C.places_at_infinity()
            [Place (1/x, 1/x^2*y, 1/x^3*y^2, 1/x^4*y^3)]

        ::

            sage: # needs sage.rings.finite_rings
            sage: F = GF(9)
            sage: A2.<x,y> = AffineSpace(F, 2)
            sage: C = A2.curve(y^3 + y - x^4)
            sage: C.places_at_infinity()
            [Place (1/x, 1/x^3*y^2)]

        ::

            sage: A.<x,y,z> = AffineSpace(GF(11), 3)
            sage: C = Curve([x*z - y^2, y - z^2, x - y*z], A)
            sage: C.places_at_infinity()
            [Place (1/x, 1/x*z^2)]
        """
        return list({p for f in self._coordinate_functions if f for p in f.poles()})

    def places_on(self, point):
        """
        Return the places on the closed point.

        INPUT:

        - ``point`` -- a closed point of the curve

        OUTPUT: list of the places of the function field of the curve

        EXAMPLES::

            sage: A.<x,y> = AffineSpace(QQ, 2)
            sage: C = Curve(x^3 - y^2 - x^4 - y^4)
            sage: C.singular_closed_points()
            [Point (x, y)]
            sage: p, = _
            sage: C.places_on(p)
            [Place (x, y, y^2, 1/x*y^3 + 1/x*y)]

        ::

            sage: # needs sage.rings.finite_rings
            sage: k.<a> = GF(9)
            sage: A.<x,y> = AffineSpace(k, 2)
            sage: C = Curve(y^2 - x^5 - x^4 - 2*x^3 - 2*x - 2)
            sage: pts = C.closed_points()
            sage: pts
            [Point (x, y + (a + 1)),
             Point (x, y + (-a - 1)),
             Point (x + (a + 1), y + (a - 1)),
             Point (x + (a + 1), y + (-a + 1)),
             Point (x - 1, y + (a + 1)),
             Point (x - 1, y + (-a - 1)),
             Point (x + (-a - 1), y + a),
             Point (x + (-a - 1), y + (-a)),
             Point (x + 1, y + 1),
             Point (x + 1, y - 1)]
            sage: p1, p2, p3 = pts[:3]
            sage: C.places_on(p1)
            [Place (x, y + a + 1)]
            sage: C.places_on(p2)
            [Place (x, y + 2*a + 2)]
            sage: C.places_on(p3)
            [Place (x + a + 1, y + a + 2)]

        ::

            sage: # needs sage.rings.finite_rings
            sage: F.<a> = GF(8)
            sage: P.<x,y,z> = ProjectiveSpace(F, 2)
            sage: Cp = Curve(x^3*y + y^3*z + x*z^3)
            sage: C = Cp.affine_patch(0)
        """
        phi = self._map_to_function_field
        gs = [phi(g) for g in point.prime_ideal().gens()]
        fs = [g for g in gs if not g.is_zero()]
        f = fs.pop()
        places = []
        for p in f.zeros():
            if all(f.valuation(p) > 0 for f in fs):
                places.append(p)
        return places

    def parametric_representation(self, place, name=None):
        """
        Return a power series representation of the branch of the
        curve given by ``place``.

        INPUT:

        - ``place`` -- a place on the curve

        EXAMPLES::

            sage: A.<x,y> = AffineSpace(QQ, 2)
            sage: C = Curve(x^2 + y^2 -1)
            sage: p = C(0,1)
            sage: p.closed_point()
            Point (x, y - 1)
            sage: pl = _.place()
            sage: C.parametric_representation(pl)
            (s + ..., 1 - 1/2*s^2 - 1/8*s^4 - 1/16*s^6 + ...)

        ::

            sage: # needs sage.rings.finite_rings
            sage: A.<x,y> = AffineSpace(GF(7^2), 2)
            sage: C = Curve(x^2 - x^4 - y^4)
            sage: p, = C.singular_closed_points()
            sage: b1, b2 = p.places()
            sage: xs, ys = C.parametric_representation(b1)
            sage: f = xs^2 - xs^4 - ys^4
            sage: [f.coefficient(i) for i in range(5)]
            [0, 0, 0, 0, 0]
            sage: xs, ys = C.parametric_representation(b2)
            sage: f = xs^2 - xs^4 - ys^4
            sage: [f.coefficient(i) for i in range(5)]
            [0, 0, 0, 0, 0]
        """
        F = place.function_field()
        F_place = F.completion(place, prec=infinity, name=name)

        return tuple(F_place._expand_lazy(c) for c in self._coordinate_functions)


class IntegralAffineCurve_finite_field(IntegralAffineCurve):
    """
    Integral affine curves.

    INPUT:

    - ``A`` -- an ambient space in which the curve lives

    - ``X`` -- list of polynomials that define the curve

    EXAMPLES::

        sage: A.<x,y,z> = AffineSpace(GF(11), 3)
        sage: C = Curve([x*z - y^2, y - z^2, x - y*z], A); C
        Affine Curve over Finite Field of size 11
         defined by -y^2 + x*z, -z^2 + y, -y*z + x
        sage: C.function_field()
        Function field in z defined by z^3 + 10*x
    """
    _point = IntegralAffineCurvePoint_finite_field

    def places(self, degree=1):
        """
        Return all places on the curve of the ``degree``.

        INPUT:

        - ``degree`` -- positive integer

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: F = GF(9)
            sage: A2.<x,y> = AffineSpace(F, 2)
            sage: C = A2.curve(y^3 + y - x^4)
            sage: C.places()
            [Place (1/x, 1/x^3*y^2),
             Place (x, y),
             Place (x, y + z2 + 1),
             Place (x, y + 2*z2 + 2),
             Place (x + z2, y + 2),
             Place (x + z2, y + z2),
             Place (x + z2, y + 2*z2 + 1),
             Place (x + z2 + 1, y + 1),
             Place (x + z2 + 1, y + z2 + 2),
             Place (x + z2 + 1, y + 2*z2),
             Place (x + 2*z2 + 1, y + 2),
             Place (x + 2*z2 + 1, y + z2),
             Place (x + 2*z2 + 1, y + 2*z2 + 1),
             Place (x + 2, y + 1),
             Place (x + 2, y + z2 + 2),
             Place (x + 2, y + 2*z2),
             Place (x + 2*z2, y + 2),
             Place (x + 2*z2, y + z2),
             Place (x + 2*z2, y + 2*z2 + 1),
             Place (x + 2*z2 + 2, y + 1),
             Place (x + 2*z2 + 2, y + z2 + 2),
             Place (x + 2*z2 + 2, y + 2*z2),
             Place (x + z2 + 2, y + 2),
             Place (x + z2 + 2, y + z2),
             Place (x + z2 + 2, y + 2*z2 + 1),
             Place (x + 1, y + 1),
             Place (x + 1, y + z2 + 2),
             Place (x + 1, y + 2*z2)]
        """
        F = self.function_field()
        return F.places(degree)

    @cached_method(do_pickle=True)
    def closed_points(self, degree=1):
        """
        Return a list of the closed points of ``degree`` of the curve.

        INPUT:

        - ``degree`` -- positive integer

        EXAMPLES::

            sage: A.<x,y> = AffineSpace(GF(7), 2)
            sage: C = Curve(x^2 - x^4 - y^4)
            sage: C.closed_points()
            [Point (x, y),
             Point (x + 1, y),
             Point (x + 2, y + 2),
             Point (x + 2, y - 2),
             Point (x - 2, y + 2),
             Point (x - 2, y - 2),
             Point (x - 1, y)]
        """
        F = self.function_field()
        places_above = F.places(degree)

        points = []

        # consider singular points
        for p in self.singular_closed_points():
            if p.degree() == degree:
                points.append(p)
            for place in p.places():
                if place.degree() == degree:
                    places_above.remove(place)

        for place in places_above:
            try:
                p = self.place_to_closed_point(place)
            except ValueError:  # place is at infinity
                continue
            assert p.degree() == degree  # sanity check
            points.append(p)

        return points


class IntegralAffinePlaneCurve(IntegralAffineCurve, AffinePlaneCurve_field):
    _point = IntegralAffinePlaneCurvePoint


class IntegralAffinePlaneCurve_finite_field(AffinePlaneCurve_finite_field, IntegralAffineCurve_finite_field):
    """
    Integral affine plane curve over a finite field.

    EXAMPLES::

        sage: # needs sage.rings.finite_rings
        sage: A.<x,y> = AffineSpace(GF(8), 2)
        sage: C = Curve(x^5 + y^5 + x*y + 1); C
        Affine Plane Curve over Finite Field in z3 of size 2^3
         defined by x^5 + y^5 + x*y + 1
        sage: C.function_field()
        Function field in y defined by y^5 + x*y + x^5 + 1
    """
    _point = IntegralAffinePlaneCurvePoint_finite_field
