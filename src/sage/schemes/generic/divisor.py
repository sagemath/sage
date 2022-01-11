"""
Divisors on schemes

AUTHORS:

- William Stein

- David Kohel

- David Joyner

- Volker Braun (2010-07-16): Documentation, doctests, coercion fixes, bugfixes.

EXAMPLES::

    sage: x,y,z = ProjectiveSpace(2, GF(5), names='x,y,z').gens()
    sage: C = Curve(y^2*z^7 - x^9 - x*z^8)
    sage: pts = C.rational_points(); pts
    [(0 : 0 : 1), (0 : 1 : 0), (2 : 2 : 1), (2 : 3 : 1), (3 : 1 : 1), (3 : 4 : 1)]
    sage: D1 = C.divisor(pts[0])*3
    sage: D2 = C.divisor(pts[1])
    sage: D3 = 10*C.divisor(pts[5])
    sage: D1.parent() is D2.parent()
    True
    sage: D = D1 - D2 + D3; D
    3*(x, y) - (x, z) + 10*(x + 2*z, y + z)
    sage: D[1][0]
    -1
    sage: D[1][1]
    Ideal (x, z) of Multivariate Polynomial Ring in x, y, z over Finite Field of size 5
    sage: C.divisor([(3, pts[0]), (-1, pts[1]), (10,pts[5])])
    3*(x, y) - (x, z) + 10*(x + 2*z, y + z)
"""
#*******************************************************************************
#  Copyright (C) 2010 Volker Braun <vbraun.name@gmail.com>
#  Copyright (C) 2005 David Kohel <kohel@maths.usyd.edu.au>
#  Copyright (C) 2005 William Stein
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*******************************************************************************

from sage.misc.latex import latex
from sage.misc.repr import repr_lincomb
from sage.misc.search import search
from sage.rings.integer_ring import ZZ
from sage.structure.formal_sum import FormalSum

from .morphism import is_SchemeMorphism
from sage.schemes.affine.affine_space import is_AffineSpace
from sage.schemes.projective.projective_space import is_ProjectiveSpace


def CurvePointToIdeal(C,P):
    r"""
    Return the vanishing ideal of a point on a curve.

    EXAMPLES::

        sage: x,y = AffineSpace(2, QQ, names='xy').gens()
        sage: C = Curve(y^2 - x^9 - x)
        sage: from sage.schemes.generic.divisor import CurvePointToIdeal
        sage: CurvePointToIdeal(C, (0,0))
        Ideal (x, y) of Multivariate Polynomial Ring in x, y over Rational Field
    """
    A = C.ambient_space()
    R = A.coordinate_ring()
    n = A.ngens()
    x = A.gens()
    polys = [ ]
    m = n-1
    while m > 0 and P[m] == 0:
        m += -1
    if is_ProjectiveSpace(A):
        a_m = P[m]
        x_m = x[m]
        for i in range(m):
            ai = P[i]
            if ai == 0:
                polys.append(x[i])
            else:
                polys.append(a_m*x[i]-ai*x_m)
    elif is_AffineSpace(A):
        for i in range(m+1):
            ai = P[i]
            if ai == 0:
                polys.append(x[i])
            else:
                polys.append(x[i]-ai)
    for i in range(m+1,n):
        polys.append(x[i])
    return R.ideal(polys)


def is_Divisor(x):
    r"""
    Test whether ``x`` is an instance of :class:`Divisor_generic`

    INPUT:

    - ``x`` -- anything.

    OUTPUT:

    ``True`` or ``False``.

    EXAMPLES::

        sage: from sage.schemes.generic.divisor import is_Divisor
        sage: x,y = AffineSpace(2, GF(5), names='xy').gens()
        sage: C = Curve(y^2 - x^9 - x)
        sage: is_Divisor( C.divisor([]) )
        True
        sage: is_Divisor("Ceci n'est pas un diviseur")
        False
    """
    return isinstance(x, Divisor_generic)


class Divisor_generic(FormalSum):
    r"""
    A Divisor.
    """

    def __init__(self, v, parent, check=True, reduce=True):
        r"""
        Construct a :class:`Divisor_generic`.

        INPUT:

        INPUT:

        - ``v`` -- object. Usually a list of pairs
          ``(coefficient,divisor)``.

        - ``parent`` -- FormalSums(R) module (default: FormalSums(ZZ))

        - ``check`` -- bool (default: True). Whether to coerce
          coefficients into base ring. Setting it to ``False`` can
          speed up construction.

        - ``reduce`` -- reduce (default: True). Whether to combine
          common terms. Setting it to ``False`` can speed up
          construction.

        .. WARNING::

            The coefficients of the divisor must be in the base ring
            and the terms must be reduced. If you set ``check=False``
            and/or ``reduce=False`` it is your responsibility to pass
            a valid object ``v``.

        EXAMPLES::

            sage: from sage.schemes.generic.divisor import Divisor_generic
            sage: from sage.schemes.generic.divisor_group import DivisorGroup
            sage: Divisor_generic( [(4,5)], DivisorGroup(Spec(ZZ)), False, False)
            4*V(5)
        """
        FormalSum.__init__(self, v, parent, check, reduce)

    def _latex_(self):
        r"""
        Return a LaTeX representation of ``self``.

        OUTPUT:

        - string.

        TESTS::

            sage: R.<x, y> = ZZ[]
            sage: S = Spec(R)
            sage: from sage.schemes.generic.divisor import Divisor_generic
            sage: from sage.schemes.generic.divisor_group import DivisorGroup
            sage: Div = DivisorGroup(S)
            sage: D = Divisor_generic([(4, x), (-5, y), (1, x+2*y)], Div)
            sage: D._latex_()
            '\\mathrm{V}\\left(x + 2 y\\right)
            + 4 \\mathrm{V}\\left(x\\right)
            - 5 \\mathrm{V}\\left(y\\right)'
        """
        # The code is copied from _repr_ with latex adjustments
        terms = list(self)
        # We sort the terms by variety. The order is "reversed" to keep it
        # straight - as the test above demonstrates, it results in the first
        # generator being in front of the second one
        terms.sort(key=lambda x: x[1], reverse=True)
        return repr_lincomb([(r"\mathrm{V}\left(%s\right)" % latex(v), c) for c,v in terms],
                            is_latex=True)

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        OUTPUT:

        - string.

        TESTS::

            sage: R.<x, y> = ZZ[]
            sage: S = Spec(R)
            sage: from sage.schemes.generic.divisor import Divisor_generic
            sage: from sage.schemes.generic.divisor_group import DivisorGroup
            sage: Div = DivisorGroup(S)
            sage: D = Divisor_generic([(4, x), (-5, y), (1, x+2*y)], Div)
            sage: D._repr_()
            'V(x + 2*y) + 4*V(x) - 5*V(y)'
        """
        # The default representation coming from formal sums does not look
        # very nice for divisors
        terms = list(self)
        # We sort the terms by variety. The order is "reversed" to keep it
        # straight - as the test above demonstrates, it results in the first
        # generator being in front of the second one
        terms.sort(key=lambda x: x[1], reverse=True)
        return repr_lincomb([("V(%s)" % v, c) for c,v in terms])

    def scheme(self):
        """
        Return the scheme that this divisor is on.

        EXAMPLES::

            sage: A.<x, y> = AffineSpace(2, GF(5))
            sage: C = Curve(y^2 - x^9 - x)
            sage: pts = C.rational_points(); pts
            [(0, 0), (2, 2), (2, 3), (3, 1), (3, 4)]
            sage: D = C.divisor(pts[0])*3 - C.divisor(pts[1]); D
            3*(x, y) - (x - 2, y - 2)
            sage: D.scheme()
            Affine Plane Curve over Finite Field of size 5 defined by -x^9 + y^2 - x
        """
        return self.parent().scheme()


class Divisor_curve(Divisor_generic):
    r"""
    For any curve `C`, use ``C.divisor(v)`` to
    construct a divisor on `C`. Here `v` can be either


    -  a rational point on `C`

    -  a list of rational points

    -  a list of 2-tuples `(c,P)`, where `c` is an
       integer and `P` is a rational point.


    TODO: Divisors shouldn't be restricted to rational points. The
    problem is that the divisor group is the formal sum of the group of
    points on the curve, and there's no implemented notion of point on
    `E/K` that has coordinates in `L`. This is what
    should be implemented, by adding an appropriate class to
    ``schemes/generic/morphism.py``.

    EXAMPLES::

        sage: E = EllipticCurve([0, 0, 1, -1, 0])
        sage: P = E(0,0)
        sage: 10*P
        (161/16 : -2065/64 : 1)
        sage: D = E.divisor(P)
        sage: D
        (x, y)
        sage: 10*D
        10*(x, y)
        sage: E.divisor([P, P])
        2*(x, y)
        sage: E.divisor([(3,P), (-4,5*P)])
        3*(x, y) - 4*(x - 1/4*z, y + 5/8*z)
    """
    def __init__(self, v, parent=None, check=True, reduce=True):
        """
        Construct a divisor on a curve.

        INPUT:

        - ``v`` -- a list of pairs ``(c, P)``, where ``c`` is an
           integer and ``P`` is a point on a curve. The P's must all
           lie on the same curve.


        - To create the divisor 0 use ``[(0, P)]``, so as to give the curve.

        EXAMPLES::

            sage: E = EllipticCurve([0, 0, 1, -1, 0])
            sage: P = E(0,0)
            sage: from sage.schemes.generic.divisor import Divisor_curve
            sage: from sage.schemes.generic.divisor_group import DivisorGroup
            sage: Divisor_curve([(1,P)], parent=DivisorGroup(E))
            (x, y)
        """
        from sage.schemes.generic.divisor_group import DivisorGroup_curve
        if not isinstance(v, (list, tuple)):
            v = [(1,v)]

        if parent is None:
            if v:
                t = v[0]
                if isinstance(t, tuple) and len(t) == 2:
                    try:
                        C = t[1].scheme()
                    except (TypeError, AttributeError):
                        raise TypeError("Argument v (= %s) must consist of multiplicities and points on a scheme.")
                else:
                    try:
                        C = t.scheme()
                    except TypeError:
                        raise TypeError("Argument v (= %s) must consist of multiplicities and points on a scheme.")
                parent = DivisorGroup_curve(C)
            else:
                raise TypeError("Argument v (= %s) must consist of multiplicities and points on a scheme.")
        else:
            if not isinstance(parent, DivisorGroup_curve):
                raise TypeError("parent (of type %s) must be a DivisorGroup_curve" % type(parent))
            C = parent.scheme()

        if len(v) < 1:
            check = False
        know_points = False
        if check:
            w = []
            points = []
            know_points = True
            for t in v:
                if isinstance(t, tuple) and len(t) == 2:
                    n = ZZ(t[0])
                    I = t[1]
                    points.append((n,I))
                else:
                    n = ZZ(1)
                    I = t
                if is_SchemeMorphism(I):
                    I = CurvePointToIdeal(C,I)
                else:
                    know_points = False
                w.append((n,I))
            v = w
        Divisor_generic.__init__(
            self, v, check=False, reduce=True, parent=parent)

        if know_points:
            self._points = points

    def _repr_(self):
        r"""
        Return a string representation.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: E = EllipticCurve([0, 0, 1, -1, 0])
            sage: E.divisor( E(0,0) )._repr_()
            '(x, y)'
        """
        return repr_lincomb([(tuple(I.gens()), c) for c, I in self])

    def support(self):
        """
        Return the support of this divisor, which is the set of points that
        occur in this divisor with nonzero coefficients.

        EXAMPLES::

            sage: x,y = AffineSpace(2, GF(5), names='xy').gens()
            sage: C = Curve(y^2 - x^9 - x)
            sage: pts = C.rational_points(); pts
            [(0, 0), (2, 2), (2, 3), (3, 1), (3, 4)]
            sage: D = C.divisor_group()([(3,pts[0]), (-1, pts[1])]); D
            3*(x, y) - (x - 2, y - 2)
            sage: D.support()
            [(0, 0), (2, 2)]

        TESTS:

        This checks that :trac:`10732` is fixed::

            sage: R.<x, y, z> = GF(5)[]
            sage: C = Curve(x^7 + y^7 + z^7)
            sage: pts = C.rational_points()
            sage: D = C.divisor([(2, pts[0])])
            sage: D.support()
            [(0 : 4 : 1)]
            sage: (D + D).support()
            [(0 : 4 : 1)]
            sage: E = C.divisor([(-3, pts[1]), (1, pts[2])])
            sage: (D - 2*E).support()
            [(0 : 4 : 1), (1 : 2 : 1), (2 : 1 : 1)]
            sage: (D - D).support()
            []
        """
        try:
            return self._support
        except AttributeError:
            try:
                pts = self._points
            except AttributeError:
                # TODO: in the next line, we should probably replace
                # rational_points() with irreducible_components()
                # once Sage can deal with divisors that are not only
                # rational points (see trac #16225)
                self._points = [(m, self.scheme().ambient_space().subscheme(p).rational_points()[0]) for (m, p) in self]
                pts = self._points
            self._support = [s[1] for s in pts]
            return self._support


    def coefficient(self, P):
        """
        Return the coefficient of a given point P in this divisor.

        EXAMPLES::

            sage: x,y = AffineSpace(2, GF(5), names='xy').gens()
            sage: C = Curve(y^2 - x^9 - x)
            sage: pts = C.rational_points(); pts
            [(0, 0), (2, 2), (2, 3), (3, 1), (3, 4)]
            sage: D = C.divisor(pts[0])
            sage: D.coefficient(pts[0])
            1
            sage: D = C.divisor([(3,pts[0]), (-1,pts[1])]); D
            3*(x, y) - (x - 2, y - 2)
            sage: D.coefficient(pts[0])
            3
            sage: D.coefficient(pts[1])
            -1
        """
        P = self.parent().scheme()(P)
        if not(P in self.support()):
            return self.base_ring().zero()
        t, i = search(self.support(), P)
        assert t
        try:
            return self._points[i][0]
        except AttributeError:
                raise NotImplementedError

