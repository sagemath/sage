r"""
Isomorphisms between Weierstrass models of elliptic curves

AUTHORS:

- Robert Bradshaw (2007): initial version
- John Cremona (Jan 2008): isomorphisms, automorphisms and twists
  in all characteristics
- Lorenz Panny (2021): :class:`EllipticCurveHom` interface
"""
# ****************************************************************************
#   Copyright (C) 2007 Robert Bradshaw <robertwb@math.washington.edu>
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

from sage.structure.element import get_coercion_model

from .constructor import EllipticCurve
from sage.schemes.elliptic_curves.hom import EllipticCurveHom
from sage.structure.richcmp import (richcmp_method, richcmp, richcmp_not_equal,
                                    op_NE)
from sage.structure.sequence import Sequence
from sage.rings.all import PolynomialRing


@richcmp_method
class baseWI(object):
    r"""
    This class implements the basic arithmetic of isomorphisms between
    Weierstrass models of elliptic curves.

    These are specified by lists of the form `[u,r,s,t]` (with
    `u\not=0`) which specifies a transformation `(x,y) \mapsto (x',y')`
    where

            `(x,y) = (u^2x'+r , u^3y' + su^2x' + t).`

    INPUT:

    - ``u,r,s,t`` (default (1,0,0,0)) -- standard parameters of an
      isomorphism between Weierstrass models.

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.weierstrass_morphism import *
        sage: baseWI()
        (1, 0, 0, 0)
        sage: baseWI(2,3,4,5)
        (2, 3, 4, 5)
        sage: R.<u,r,s,t> = QQ[]
        sage: baseWI(u,r,s,t)
        (u, r, s, t)
    """
    def __init__(self, u=1, r=0, s=0, t=0):
        r"""
        Constructor: check for valid parameters (defaults to identity)

        INPUT:

        - ``u,r,s,t`` (default (1,0,0,0)) -- standard parameters of an
          isomorphism between Weierstrass models.

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.weierstrass_morphism import *
            sage: baseWI()
            (1, 0, 0, 0)
            sage: baseWI(2,3,4,5)
            (2, 3, 4, 5)
            sage: R.<u,r,s,t> = QQ[]
            sage: baseWI(u,r,s,t)
            (u, r, s, t)
        """
        if u == 0:
            raise ValueError("u!=0 required for baseWI")
        self.u = u
        self.r = r
        self.s = s
        self.t = t

    def __richcmp__(self, other, op):
        """
        Standard comparison function.

        The ordering is just lexicographic on the tuple `(u,r,s,t)`.

        .. NOTE::

            In a list of automorphisms, there is no guarantee that the
            identity will be first!

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.weierstrass_morphism import baseWI
            sage: baseWI(1,2,3,4) == baseWI(1,2,3,4)
            True
            sage: baseWI(1,2,3,4) != baseWI(1,2,3,4)
            False
            sage: baseWI(1,2,3,4) < baseWI(1,2,3,5)
            True
            sage: baseWI(1,2,3,4) > baseWI(1,2,3,4)
            False

        It will never return equality if ``other`` is of another type::

            sage: baseWI() == 1
            False
        """
        if not isinstance(other, baseWI):
            return (op == op_NE)
        return richcmp(self.tuple(), other.tuple(), op)

    def tuple(self):
        r"""
        Return the parameters `u,r,s,t` as a tuple.

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.weierstrass_morphism import *
            sage: w = baseWI(2,3,4,5)
            sage: w.tuple()
            (2, 3, 4, 5)
        """
        return (self.u, self.r, self.s, self.t)

    def __mul__(self, other):
        r"""
        Return the composition of this isomorphism and another.

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.weierstrass_morphism import *
            sage: baseWI(1,2,3,4)*baseWI(5,6,7,8)
            (5, 56, 22, 858)
            sage: baseWI()*baseWI(1,2,3,4)*baseWI()
            (1, 2, 3, 4)
        """
        u1, r1, s1, t1 = other.tuple()
        u2, r2, s2, t2 = self.tuple()
        return baseWI(u1 * u2,
                      (u1**2) * r2 + r1,
                      u1 * s2 + s1,
                      (u1**3) * t2 + s1 * (u1**2) * r2 + t1)

    def __invert__(self):
        r"""
        Return the inverse of this isomorphism.

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.weierstrass_morphism import *
            sage: w = baseWI(2,3,4,5)
            sage: ~w
            (1/2, -3/4, -2, 7/8)
            sage: w*~w
            (1, 0, 0, 0)
            sage: ~w*w
            (1, 0, 0, 0)
            sage: R.<u,r,s,t> = QQ[]
            sage: w = baseWI(u,r,s,t)
            sage: ~w
            (1/u, (-r)/u^2, (-s)/u, (r*s - t)/u^3)
            sage: ~w*w
            (1, 0, 0, 0)
        """
        u, r, s, t = self.tuple()
        return baseWI(1/u, -r/(u**2), -s/u, (r*s-t)/(u**3))

    def __repr__(self):
        r"""
        Return the string representation  of this isomorphism.

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.weierstrass_morphism import *
            sage: baseWI(2,3,4,5)
            (2, 3, 4, 5)
        """
        return repr(self.tuple())

    def is_identity(self):
        r"""
        Return ``True`` if this is the identity isomorphism.

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.weierstrass_morphism import *
            sage: w = baseWI(); w.is_identity()
            True
            sage: w = baseWI(2,3,4,5); w.is_identity()
            False
        """
        return self.tuple() == (1, 0, 0, 0)

    def __call__(self, EorP):
        r"""
        Base application of isomorphisms to curves and points.

        A baseWI `w` may be applied to a list `[a1,a2,a3,a4,a6]`
        representing the `a`-invariants of an elliptic curve `E`,
        returning the `a`-invariants of `w(E)`; or to `P=[x,y]` or
        `P=[x,y,z]` representing a point in `\mathbb{A}^2` or
        `\mathbb{P}^2`, returning the transformed point.

        INPUT:

        - ``EorP`` -- either an elliptic curve, or a point on an elliptic curve.

        OUTPUT:

        The transformed curve or point.

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.weierstrass_morphism import *
            sage: E = EllipticCurve([0,0,1,-7,6])
            sage: w = baseWI(2,3,4,5)
            sage: w(E.ainvs())
            [4, -7/4, 11/8, -3/2, -9/32]
            sage: P = E(-2,3)
            sage: w(P.xy())
            [-5/4, 9/4]
            sage: EllipticCurve(w(E.ainvs()))(w(P.xy()))
            (-5/4 : 9/4 : 1)
        """
        u, r, s, t = self.tuple()
        if len(EorP) == 5:
            a1, a2, a3, a4, a6 = EorP
            a6 += r*(a4 + r*(a2 + r)) - t*(a3 + r*a1 + t)
            a4 += -s*a3 + 2*r*a2 - (t + r*s)*a1 + 3*r*r - 2*s*t
            a3 += r*a1 + t + t
            a2 += -s*a1 + 3*r - s*s
            a1 += 2*s
            return [a1/u, a2/u**2, a3/u**3, a4/u**4, a6/u**6]
        if len(EorP) == 2:
            x, y = EorP
            x -= r
            y -= (s*x+t)
            return [x/u**2, y/u**3]
        if len(EorP) == 3:
            x, y, z = EorP
            x -= r*z
            y -= (s*x+t*z)
            return [x/u**2, y/u**3, z]
        raise ValueError("baseWI(a) only for a=(x,y), (x:y:z) or (a1,a2,a3,a4,a6)")


def isomorphisms(E, F, JustOne=False):
    r"""
    Return one or all isomorphisms between two elliptic curves.

    INPUT:

    - ``E``, ``F`` (EllipticCurve) -- Two elliptic curves.

    - ``JustOne`` (bool) If ``True``, returns one isomorphism, or ``None`` if
      the curves are not isomorphic.  If ``False``, returns a (possibly
      empty) list of isomorphisms.

    OUTPUT:

    Either ``None``, or a 4-tuple `(u,r,s,t)` representing an isomorphism,
    or a list of these.

    .. note::

        This function is not intended for users, who should use the
        interface provided by ``ell_generic``.

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.weierstrass_morphism import *
        sage: isomorphisms(EllipticCurve_from_j(0),EllipticCurve('27a3'))
        [(-1, 0, 0, -1), (1, 0, 0, 0)]
        sage: isomorphisms(EllipticCurve_from_j(0),EllipticCurve('27a3'),JustOne=True)
        (1, 0, 0, 0)
        sage: isomorphisms(EllipticCurve_from_j(0),EllipticCurve('27a1'))
        []
        sage: isomorphisms(EllipticCurve_from_j(0),EllipticCurve('27a1'),JustOne=True)

    TESTS:

    Check that :trac:`32632` is fixed::

        sage: z8 = GF(2^8).gen()
        sage: E1 = EllipticCurve([z8, z8, z8, z8, z8])
        sage: isomorphisms(E1, E1)
        [(1, 0, 0, 0), (1, 0, z8, z8)]
        sage: E2 = EllipticCurve([z8^2, 0, 0, 0, z8^7 + z8^4])
        sage: isomorphisms(E1, E2)
        [(z8^7 + z8^3 + z8^2 + z8, 1, 1, z8^7 + z8^3 + z8^2 + z8 + 1),
         (z8^7 + z8^3 + z8^2 + z8, 1, z8 + 1, z8^7 + z8^3 + z8^2 + z8 + 1)]
    """
    from .ell_generic import is_EllipticCurve
    if not is_EllipticCurve(E) or not is_EllipticCurve(F):
        raise ValueError("arguments are not elliptic curves")
    K = E.base_ring()

    j = E.j_invariant()
    if j != F.j_invariant():
        if JustOne:
            return None
        return []

    from sage.rings.polynomial.polynomial_ring import polygen
    x = polygen(K, 'x')

    a1E, a2E, a3E, a4E, a6E = E.ainvs()
    a1F, a2F, a3F, a4F, a6F = F.ainvs()

    char = K.characteristic()

    if char == 2:
        if j == 0:
            ulist = (x**3-(a3E/a3F)).roots(multiplicities=False)
            ans = []
            for u in ulist:
                slist = (x**4+a3E*x+(a2F**2+a4F)*u**4+a2E**2+a4E).roots(multiplicities=False)
                for s in slist:
                    r = s**2+a2E+a2F*u**2
                    tlist = (x**2 + a3E*x + r**3 + a2E*r**2 + a4E*r + a6E + a6F*u**6).roots(multiplicities=False)
                    for t in tlist:
                        if JustOne:
                            return (u, r, s, t)
                        ans.append((u, r, s, t))
            if JustOne:
                return None
            ans.sort()
            return ans
        else:
            ans = []
            u = a1E/a1F
            r = (a3E+a3F*u**3)/a1E
            slist = [s[0] for s in (x**2+a1E*x+(r+a2E+a2F*u**2)).roots()]
            for s in slist:
                t = (a4E+a4F*u**4 + s*a3E + r*s*a1E + r**2) / a1E
                if JustOne:
                    return (u, r, s, t)
                ans.append((u, r, s, t))
            if JustOne:
                return None
            ans.sort()
            return ans

    b2E, b4E, b6E, b8E = E.b_invariants()
    b2F, b4F, b6F, b8F = F.b_invariants()

    if char == 3:
        if j == 0:
            ulist = (x**4-(b4E/b4F)).roots(multiplicities=False)
            ans = []
            for u in ulist:
                s = a1E-a1F*u
                t = a3E-a3F*u**3
                rlist = (x**3-b4E*x+(b6E-b6F*u**6)).roots(multiplicities=False)
                for r in rlist:
                    if JustOne:
                        return (u, r, s, t+r*a1E)
                    ans.append((u, r, s, t+r*a1E))
            if JustOne:
                return None
            ans.sort()
            return ans
        else:
            ulist = (x**2 - b2E / b2F).roots(multiplicities=False)
            ans = []
            for u in ulist:
                r = (b4F * u**4 - b4E) / b2E
                s = (a1E - a1F * u)
                t = (a3E - a3F * u**3 + a1E * r)
                if JustOne:
                    return (u, r, s, t)
                ans.append((u, r, s, t))
            if JustOne:
                return None
            ans.sort()
            return ans

# now char!=2,3:
    c4E, c6E = E.c_invariants()
    c4F, c6F = F.c_invariants()

    if j == 0:
        m, um = 6, c6E/c6F
    elif j == 1728:
        m, um = 4, c4E/c4F
    else:
        m, um = 2, (c6E*c4F)/(c6F*c4E)
    ulist = (x**m-um).roots(multiplicities=False)
    ans = []
    for u in ulist:
        s = (a1F*u - a1E)/2
        r = (a2F*u**2 + a1E*s + s**2 - a2E)/3
        t = (a3F*u**3 - a1E*r - a3E)/2
        if JustOne:
            return (u, r, s, t)
        ans.append((u, r, s, t))
    if JustOne:
        return None
    ans.sort()
    return ans


class WeierstrassIsomorphism(EllipticCurveHom, baseWI):
    r"""
    Class representing a Weierstrass isomorphism between two elliptic curves.
    """
    def __init__(self, E=None, urst=None, F=None):
        r"""
        Constructor for WeierstrassIsomorphism class,

        INPUT:

        - ``E`` -- an EllipticCurve, or None (see below).

        - ``urst`` -- a 4-tuple `(u,r,s,t)`, or None (see below).

        - ``F`` -- an EllipticCurve, or None (see below).

        Given two Elliptic Curves ``E`` and ``F`` (represented by
        Weierstrass models as usual), and a transformation ``urst``
        from ``E`` to ``F``, construct an isomorphism from ``E`` to
        ``F``.  An exception is raised if ``urst(E)!=F``.  At most one
        of ``E``, ``F``, ``urst`` can be None.  If ``F==None`` then
        ``F`` is constructed as ``urst(E)``.  If ``E==None`` then
        ``E`` is constructed as ``urst^-1(F)``.  If ``urst==None``
        then an isomorphism from ``E`` to ``F`` is constructed if
        possible, and an exception is raised if they are not
        isomorphic.  Otherwise ``urst`` can be a tuple of length 4 or
        a object of type ``baseWI``.

        Users will not usually need to use this class directly, but instead use
        methods such as ``isomorphism`` of elliptic curves.

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.weierstrass_morphism import *
            sage: WeierstrassIsomorphism(EllipticCurve([0,1,2,3,4]),(-1,2,3,4))
            Elliptic-curve morphism:
            From: Elliptic Curve defined by y^2 + 2*y = x^3 + x^2 + 3*x + 4 over Rational Field
            To:   Elliptic Curve defined by y^2 - 6*x*y - 10*y = x^3 - 2*x^2 - 11*x - 2 over Rational Field
            Via:  (u,r,s,t) = (-1, 2, 3, 4)
            sage: E = EllipticCurve([0,1,2,3,4])
            sage: F = EllipticCurve(E.cremona_label())
            sage: WeierstrassIsomorphism(E,None,F)
            Elliptic-curve morphism:
            From: Elliptic Curve defined by y^2 + 2*y = x^3 + x^2 + 3*x + 4 over Rational Field
            To:   Elliptic Curve defined by y^2  = x^3 + x^2 + 3*x + 5 over Rational Field
            Via:  (u,r,s,t) = (1, 0, 0, -1)
            sage: w = WeierstrassIsomorphism(None,(1,0,0,-1),F)
            sage: w._domain==E
            True
        """
        from .ell_generic import is_EllipticCurve

        if E is not None:
            if not is_EllipticCurve(E):
                raise ValueError("First argument must be an elliptic curve or None")
        if F is not None:
            if not is_EllipticCurve(F):
                raise ValueError("Third argument must be an elliptic curve or None")
        if urst is not None:
            if len(urst) != 4:
                raise ValueError("Second argument must be [u,r,s,t] or None")
        if len([par for par in [E, urst, F] if par is not None]) < 2:
            raise ValueError("At most 1 argument can be None")

        if F is None:  # easy case
            baseWI.__init__(self, *urst)
            F = EllipticCurve(baseWI.__call__(self, list(E.a_invariants())))

        elif E is None:  # easy case in reverse
            baseWI.__init__(self, *urst)
            inv_urst = baseWI.__invert__(self)
            E = EllipticCurve(baseWI.__call__(inv_urst, list(F.a_invariants())))

        elif urst is None:  # try to construct the morphism
            urst = isomorphisms(E, F, True)
            if urst is None:
                raise ValueError("Elliptic curves not isomorphic.")
            baseWI.__init__(self, *urst)

        else:  # none of the parameters is None:
            baseWI.__init__(self, *urst)
            if F != EllipticCurve(baseWI.__call__(self, list(E.a_invariants()))):
                raise ValueError("second argument is not an isomorphism from first argument to third argument")

        base_ring = get_coercion_model().common_parent(E.base_ring(), F.base_ring(), *urst)
        self._mpoly_ring = PolynomialRing(base_ring, ['x','y'])
        self._poly_ring = PolynomialRing(base_ring, ['x'])

        self._domain = E
        self._codomain = F
        EllipticCurveHom.__init__(self, self._domain, self._codomain)

    def _richcmp_(self, other, op):
        r"""
        Standard comparison function for the WeierstrassIsomorphism class.

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.weierstrass_morphism import *
            sage: E = EllipticCurve('389a1')
            sage: F = E.change_weierstrass_model(1,2,3,4)
            sage: w1 = E.isomorphism_to(F)
            sage: w1 == w1
            True
            sage: w2 = F.automorphisms()[0] *w1
            sage: w1 == w2
            False

            sage: E = EllipticCurve_from_j(GF(7)(0))
            sage: F = E.change_weierstrass_model(2,3,4,5)
            sage: a = E.isomorphisms(F)
            sage: b = [w*a[0] for w in F.automorphisms()]
            sage: b.sort()
            sage: a == b
            True
            sage: c = [a[0]*w for w in E.automorphisms()]
            sage: c.sort()
            sage: a == c
            True
        """
        if isinstance(other, WeierstrassIsomorphism):
            lx = self._domain
            rx = other._domain
            if lx != rx:
                return richcmp_not_equal(lx, rx, op)

            lx = self._codomain
            rx = other._codomain
            if lx != rx:
                return richcmp_not_equal(lx, rx, op)

            return baseWI.__richcmp__(self, other, op)

        return EllipticCurveHom._richcmp_(self, other, op)

    def _eval(self, P):
        r"""
        Less strict evaluation method for internal use.

        In particular, this can be used to evaluate ``self`` at a
        point defined over an extension field.

        INPUT: a sequence of 3 coordinates defining a point on ``self``

        OUTPUT: the result of evaluating ``self'' at the given point

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.weierstrass_morphism import WeierstrassIsomorphism
            sage: E = EllipticCurve([i,0]); E
            Elliptic Curve defined by y^2 = x^3 + I*x over Number Field in I with defining polynomial x^2 + 1 with I = 1*I
            sage: iso = WeierstrassIsomorphism(E, (i,1,2,3))
            sage: P = E.change_ring(QQbar).lift_x(QQbar.random_element())
            sage: Q = iso._eval(P)
            sage: Q.curve()
            Elliptic Curve defined by y^2 + (-4*I)*x*y + 6*I*y = x^3 + x^2 + (I-9)*x + (-I+8) over Algebraic Field
            sage: y = next(filter(bool, iter(QQbar.random_element, None)))  # sample until nonzero
            sage: iso._eval((0, y, 0)) == 0
            True
        """
        if self._domain.defining_polynomial()(*P):
            raise ValueError(f'{P} not on {self._domain}')

        Q = baseWI.__call__(self, P)
        k = Sequence(tuple(P) + tuple(Q)).universe()
        return self._codomain.base_extend(k).point(Q)

    def __call__(self, P):
        r"""
        Call function for WeierstrassIsomorphism class.

        INPUT:

        - ``P`` (Point) -- a point on the domain curve.

        OUTPUT:

        (Point) the transformed point on the codomain curve.

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.weierstrass_morphism import *
            sage: E = EllipticCurve('37a1')
            sage: w = WeierstrassIsomorphism(E,(2,3,4,5))
            sage: P = E(0,-1)
            sage: w(P)
            (-3/4 : 3/4 : 1)
            sage: w(P).curve() == E.change_weierstrass_model((2,3,4,5))
            True
        """
        if P[2] == 0:
            return self._codomain(0)
        return self._codomain.point(baseWI.__call__(self,
                                                          tuple(P._coords)),
                                          check=False)

    def __invert__(self):
        r"""
        Return the inverse of this WeierstrassIsomorphism.

        EXAMPLES::

            sage: E = EllipticCurve('5077')
            sage: F = E.change_weierstrass_model([2,3,4,5]); F
            Elliptic Curve defined by y^2 + 4*x*y + 11/8*y = x^3 - 7/4*x^2 - 3/2*x - 9/32 over Rational Field
            sage: w = E.isomorphism_to(F)
            sage: P = E(-2,3,1)
            sage: w(P)
            (-5/4 : 9/4 : 1)
            sage: ~w
            Elliptic-curve morphism:
              From: Elliptic Curve defined by y^2 + 4*x*y + 11/8*y = x^3 - 7/4*x^2 - 3/2*x - 9/32 over Rational Field
              To:   Elliptic Curve defined by y^2 + y = x^3 - 7*x + 6 over Rational Field
              Via:  (u,r,s,t) = (1/2, -3/4, -2, 7/8)
            sage: Q = w(P); Q
            (-5/4 : 9/4 : 1)
            sage: (~w)(Q)
            (-2 : 3 : 1)
        """
        winv = baseWI.__invert__(self).tuple()
        return WeierstrassIsomorphism(self._codomain, winv, self._domain)

    @staticmethod
    def _composition_impl(left, right):
        r"""
        Return the composition of a ``WeierstrassIsomorphism``
        with another elliptic-curve morphism.

        Called by :meth:`EllipticCurveHom._composition_`.

        EXAMPLES::

            sage: E1 = EllipticCurve('5077')
            sage: E2 = E1.change_weierstrass_model([2,3,4,5])
            sage: w1 = E1.isomorphism_to(E2)
            sage: E3 = E2.change_weierstrass_model([6,7,8,9])
            sage: w2 = E2.isomorphism_to(E3)
            sage: P = E1(-2,3,1)
            sage: (w2*w1)(P) == w2(w1(P))
            True

        TESTS:

        We should return ``NotImplemented`` when passed a combination of
        elliptic-curve morphism types that we don't handle here::

            sage: E = EllipticCurve([1,0])
            sage: phi = E.isogeny(E(0,0))
            sage: w1._composition_impl(phi.dual(), phi)
            NotImplemented
        """
        if isinstance(left, WeierstrassIsomorphism) and isinstance(right, WeierstrassIsomorphism):
            if left._domain != right._codomain:
                raise ValueError("Domain of first argument must equal codomain of second")
            w = baseWI.__mul__(left, right)
            return WeierstrassIsomorphism(right._domain, w.tuple(), left._codomain)

        return NotImplemented

    def __repr__(self):
        r"""
        Return the string representation of this WeierstrassIsomorphism.

        OUTPUT:

        (string) The underlying morphism, together with an extra line
        showing the `(u,r,s,t)` parameters.

        EXAMPLES::

            sage: E1 = EllipticCurve('5077')
            sage: E2 = E1.change_weierstrass_model([2,3,4,5])
            sage: E1.isomorphism_to(E2)
            Elliptic-curve morphism:
            From: Elliptic Curve defined by y^2 + y = x^3 - 7*x + 6 over Rational Field
            To:   Elliptic Curve defined by y^2 + 4*x*y + 11/8*y = x^3 - 7/4*x^2 - 3/2*x - 9/32 over Rational Field
            Via:  (u,r,s,t) = (2, 3, 4, 5)
        """
        return EllipticCurveHom.__repr__(self) + "\n  Via:  (u,r,s,t) = " + baseWI.__repr__(self)

    # EllipticCurveHom methods

    def degree(self):
        """
        Return the degree as a rational map of this isomorphism.

        Isomorphisms always have degree `1` by definition.

        EXAMPLES::

            sage: E1 = EllipticCurve([1,2,3,4,5])
            sage: E2 = EllipticCurve_from_j(E1.j_invariant())
            sage: E1.isomorphism_to(E2).degree()
            1
        """
        return 1

    def rational_maps(self):
        """
        Return the pair of rational maps defining this isomorphism.

        EXAMPLES::

            sage: E1 = EllipticCurve([11,22,33,44,55])
            sage: E2 = EllipticCurve_from_j(E1.j_invariant())
            sage: iso = E1.isomorphism_to(E2); iso
            Elliptic-curve morphism:
              From: Elliptic Curve defined by y^2 + 11*x*y + 33*y = x^3 + 22*x^2 + 44*x + 55 over Rational Field
              To:   Elliptic Curve defined by y^2 + x*y = x^3 + x^2 - 684*x + 6681 over Rational Field
              Via:  (u,r,s,t) = (1, -17, -5, 77)
            sage: iso.rational_maps()
            (x + 17, 5*x + y + 8)
            sage: f = E2.defining_polynomial()(*iso.rational_maps(), 1)
            sage: I = E1.defining_ideal()
            sage: x,y,z = I.ring().gens()
            sage: f in I + Ideal(z-1)
            True

        ::

            sage: E = EllipticCurve(GF(65537), [1,1,1,1,1])
            sage: w = E.isomorphism_to(E.short_weierstrass_model())
            sage: f,g = w.rational_maps()
            sage: P = E.random_point()
            sage: w(P).xy() == (f(P.xy()), g(P.xy()))
            True

        TESTS::

            sage: iso.rational_maps()[0].parent()
            Multivariate Polynomial Ring in x, y over Rational Field
            sage: iso.rational_maps()[1].parent()
            Multivariate Polynomial Ring in x, y over Rational Field
        """
        return tuple(baseWI.__call__(self, self._mpoly_ring.gens()))

    def x_rational_map(self):
        """
        Return the `x`-coordinate rational map of this isomorphism.

        EXAMPLES::

            sage: E1 = EllipticCurve([11,22,33,44,55])
            sage: E2 = EllipticCurve_from_j(E1.j_invariant())
            sage: iso = E1.isomorphism_to(E2); iso
            Elliptic-curve morphism:
              From: Elliptic Curve defined by y^2 + 11*x*y + 33*y = x^3 + 22*x^2 + 44*x + 55 over Rational Field
              To:   Elliptic Curve defined by y^2 + x*y = x^3 + x^2 - 684*x + 6681 over Rational Field
              Via:  (u,r,s,t) = (1, -17, -5, 77)
            sage: iso.x_rational_map()
            x + 17
            sage: iso.x_rational_map() == iso.rational_maps()[0]
            True

        TESTS::

            sage: iso.x_rational_map().parent()
            Univariate Polynomial Ring in x over Rational Field
        """
        x, = self._poly_ring.gens()
        return (x - self.r) / self.u**2

    def kernel_polynomial(self):
        """
        Return the kernel polynomial of this isomorphism.

        Isomorphisms have trivial kernel by definition, hence this
        method always returns `1`.

        EXAMPLES::

            sage: E1 = EllipticCurve([11,22,33,44,55])
            sage: E2 = EllipticCurve_from_j(E1.j_invariant())
            sage: iso = E1.isomorphism_to(E2)
            sage: iso.kernel_polynomial()
            1
            sage: psi = E1.isogeny(iso.kernel_polynomial(), codomain=E2); psi
            Isogeny of degree 1 from Elliptic Curve defined by y^2 + 11*x*y + 33*y = x^3 + 22*x^2 + 44*x + 55 over Rational Field to Elliptic Curve defined by y^2 + x*y = x^3 + x^2 - 684*x + 6681 over Rational Field
            sage: psi in {iso, -iso}
            True

        TESTS::

            sage: iso.kernel_polynomial().parent()
            Univariate Polynomial Ring in x over Rational Field
        """
        return self._poly_ring(1)

    def dual(self):
        """
        Return the dual isogeny of this isomorphism.

        For isomorphisms, the dual is just the inverse.

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.weierstrass_morphism import WeierstrassIsomorphism
            sage: E = EllipticCurve(QuadraticField(-3), [0,1])
            sage: w = WeierstrassIsomorphism(E, (CyclotomicField(3).gen(),0,0,0))
            sage: (w.dual() * w).rational_maps()
            (x, y)

        ::

            sage: E1 = EllipticCurve([11,22,33,44,55])
            sage: E2 = E1.short_weierstrass_model()
            sage: iso = E1.isomorphism_to(E2)
            sage: iso.dual() == ~iso
            True
        """
        return ~self

    def __neg__(self):
        """
        Return the negative of this isomorphism, i.e., its composition
        with the negation map `[-1]`.

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.weierstrass_morphism import WeierstrassIsomorphism
            sage: E = EllipticCurve([11,22,33,44,55])
            sage: w = WeierstrassIsomorphism(E, (66,77,88,99))
            sage: -w
            Elliptic-curve morphism:
              From: Elliptic Curve defined by y^2 + 11*x*y + 33*y = x^3 + 22*x^2 + 44*x + 55 over Rational Field
              To:   Elliptic Curve defined by y^2 + 17/6*x*y + 49/13068*y = x^3 - 769/396*x^2 - 3397/862488*x + 44863/7513995456 over Rational Field
              Via:  (u,r,s,t) = (-66, 77, -99, -979)
            sage: -(-w) == w
            True

        ::

            sage: from sage.schemes.elliptic_curves.weierstrass_morphism import WeierstrassIsomorphism
            sage: E = EllipticCurve(QuadraticField(-3), [0,1])
            sage: w = WeierstrassIsomorphism(E, (CyclotomicField(3).gen(),0,0,0))
            sage: w.tuple()
            (zeta3, 0, 0, 0)
            sage: (-w).tuple()
            (-zeta3, 0, 0, 0)
            sage: (-w)^3 == -(w^3)
            True

        ::

            sage: from sage.schemes.elliptic_curves.weierstrass_morphism import WeierstrassIsomorphism
            sage: E = EllipticCurve(QuadraticField(-1), [1,0])
            sage: t = WeierstrassIsomorphism(E, (i,0,0,0))
            sage: -t^2 == WeierstrassIsomorphism(E, (1,0,0,0))
            True
        """
        a1,_,a3,_,_ = self._domain.a_invariants()
        w = baseWI(-1, 0, -a1, -a3)
        urst = baseWI.__mul__(self, w).tuple()
        return WeierstrassIsomorphism(self._domain, urst, self._codomain)

