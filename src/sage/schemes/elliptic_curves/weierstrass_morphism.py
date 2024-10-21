# sage_setup: distribution = sagemath-schemes
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
from sage.structure.richcmp import (richcmp, richcmp_not_equal, op_EQ, op_NE)
from sage.structure.sequence import Sequence
from sage.rings.integer import Integer
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing


class baseWI:
    r"""
    This class implements the basic arithmetic of isomorphisms between
    Weierstrass models of elliptic curves.

    These are specified by lists of the form `[u,r,s,t]` (with `u \neq 0`)
    which specifies a transformation `(x,y) \mapsto (x',y')` where

            `(x,y) = (u^2x'+r , u^3y' + su^2x' + t).`

    INPUT:

    - ``u``, ``r``, ``s``, ``t`` -- (default: `1`, `0`, `0`, `0`); standard
      parameters of an isomorphism between Weierstrass models

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
        Constructor: check for valid parameters (defaults to identity).

        INPUT:

        - ``u``, ``r``, ``s``, ``t`` -- (default: `1`, `0`, `0`, `0`); standard
          parameters of an isomorphism between Weierstrass models

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
        if not u:
            raise ValueError("u!=0 required for baseWI")
        self.u = u
        self.r = r
        self.s = s
        self.t = t

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
        return baseWI(1/u, -r/u**2, -s/u, (r*s-t)/u**3)

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

        - ``EorP`` -- either an elliptic curve, or a point on an elliptic curve

        OUTPUT: the transformed curve or point

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


def _isomorphisms(E, F):
    r"""
    Enumerate all isomorphisms between two elliptic curves,
    as a generator object.

    INPUT:

    - ``E``, ``F`` -- two elliptic curves

    OUTPUT:

    A generator object producing 4-tuples `(u,r,s,t)` representing an isomorphism.

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.weierstrass_morphism import _isomorphisms
        sage: list(_isomorphisms(EllipticCurve_from_j(0), EllipticCurve('27a3')))
        [(1, 0, 0, 0), (-1, 0, 0, -1)]
        sage: list(_isomorphisms(EllipticCurve_from_j(0), EllipticCurve('27a1')))
        []

    TESTS:

    Check that :issue:`32632` is fixed::

        sage: z8 = GF(2^8).gen()
        sage: E1 = EllipticCurve([z8, z8, z8, z8, z8])
        sage: list(_isomorphisms(E1, E1))
        [(1, 0, 0, 0), (1, 0, z8, z8)]
        sage: E2 = EllipticCurve([z8^2, 0, 0, 0, z8^7 + z8^4])
        sage: list(_isomorphisms(E1, E2))
        [(z8^7 + z8^3 + z8^2 + z8, 1, 1, z8^7 + z8^3 + z8^2 + z8 + 1),
         (z8^7 + z8^3 + z8^2 + z8, 1, z8 + 1, z8^7 + z8^3 + z8^2 + z8 + 1)]

    Random testing::

        sage: p = random_prime(100)
        sage: F = GF(p).algebraic_closure()
        sage: j = F.random_element()
        sage: while j in (0, 1728):     # skip the hard case
        ....:     j = F.random_element()
        sage: j = F(choice((0, 1728)))  # long time -- do the hard case
        sage: E = EllipticCurve_from_j(j)
        sage: u,r,s,t = (F^4).random_element()
        sage: u = u or 1
        sage: E = E.change_weierstrass_model(u,r,s,t)
        sage: Aut = E.automorphisms()
        sage: len(set(Aut)) == len(Aut)
        True
        sage: all(-a in Aut for a in Aut)
        True
        sage: len(Aut) in (2, 4, 6, 12, 24)
        True
        sage: j = E.j_invariant()
        sage: {
        ....:      2: j not in (0, 1728),
        ....:      4: p >= 5 and j == 1728,
        ....:      6: p >= 5 and j == 0,
        ....:     12: p == 3 and j == 0,  # note 1728 == 0
        ....:     24: p == 2 and j == 0,  # note 1728 == 0
        ....: }[len(Aut)]
        True
        sage: u,r,s,t = (F^4).random_element()
        sage: u = u or 1
        sage: F = E.change_weierstrass_model(u,r,s,t)
        sage: Iso = E.isomorphisms(F)
        sage: len(set(Iso)) == len(Iso)
        True
        sage: all(-f in Iso for f in Iso)
        True
        sage: len(Iso) == len(Aut)
        True
        sage: all({iso2*iso1 for iso1 in Iso} == set(Aut) for iso2 in F.isomorphisms(E))
        True
    """
    from .ell_generic import EllipticCurve_generic
    if not isinstance(E, EllipticCurve_generic) or not isinstance(F, EllipticCurve_generic):
        raise ValueError("arguments are not elliptic curves")

    j = E.j_invariant()
    if j != F.j_invariant():
        return

    K = E.base_ring()

    from sage.rings.polynomial.polynomial_ring import polygen
    x = polygen(K, 'x')

    a1E, a2E, a3E, a4E, a6E = E.ainvs()
    a1F, a2F, a3F, a4F, a6F = F.ainvs()

    char = K.characteristic()

    if char == 2:
        if j == 0:
            ulist = (x**3 - a3E/a3F).roots(multiplicities=False)
            for u in ulist:
                slist = (x**4 + a3E*x + (a2F**2 + a4F)*u**4 + a2E**2 + a4E).roots(multiplicities=False)
                for s in slist:
                    r = s**2 + a2E + a2F*u**2
                    tlist = (x**2 + a3E*x + r**3 + a2E*r**2 + a4E*r + a6E + a6F*u**6).roots(multiplicities=False)
                    for t in tlist:
                        yield (u, r, s, t)
        else:
            u = a1E/a1F
            r = (a3E + a3F*u**3)/a1E
            slist = (x**2 + a1E*x + r + a2E + a2F*u**2).roots(multiplicities=False)
            for s in slist:
                t = (a4E + a4F*u**4 + s*a3E + r*s*a1E + r**2) / a1E
                yield (u, r, s, t)
        return

    b2E, b4E, b6E, b8E = E.b_invariants()
    b2F, b4F, b6F, b8F = F.b_invariants()

    if char == 3:
        if j == 0:
            ulist = (x**4 - b4E/b4F).roots(multiplicities=False)
            for u in ulist:
                s = a1E - a1F*u
                t = a3E - a3F*u**3
                rlist = (x**3 - b4E*x + b6E - b6F*u**6).roots(multiplicities=False)
                for r in rlist:
                    yield (u, r, s, t + r*a1E)
        else:
            ulist = (x**2 - b2E/b2F).roots(multiplicities=False)
            for u in ulist:
                r = (b4F * u**4 - b4E) / b2E
                s = a1E - a1F * u
                t = a3E - a3F * u**3 + a1E * r
                yield (u, r, s, t)
        return

    # now char != 2,3:

    c4E, c6E = E.c_invariants()
    c4F, c6F = F.c_invariants()

    if j == 0:
        m, um = 6, c6E/c6F
    elif j == 1728:
        m, um = 4, c4E/c4F
    else:
        m, um = 2, (c6E*c4F)/(c6F*c4E)
    ulist = (x**m - um).roots(multiplicities=False)
    for u in ulist:
        s = (a1F*u - a1E)/2
        r = (a2F*u**2 + a1E*s + s**2 - a2E)/3
        t = (a3F*u**3 - a1E*r - a3E)/2
        yield (u, r, s, t)


class WeierstrassIsomorphism(EllipticCurveHom, baseWI):
    r"""
    Class representing a Weierstrass isomorphism between two elliptic curves.

    INPUT:

    - ``E`` -- an ``EllipticCurve``, or ``None`` (see below)

    - ``urst`` -- a 4-tuple `(u,r,s,t)`, a :class:`baseWI` object,
      or ``None`` (see below)

    - ``F`` -- an ``EllipticCurve``, or ``None`` (see below)

    Given two Elliptic Curves ``E`` and ``F`` (represented by Weierstrass
    models as usual), and a transformation ``urst`` from ``E`` to ``F``,
    construct an isomorphism from ``E`` to ``F``.
    An exception is raised if ``urst(E) != F``.  At most one of ``E``,
    ``F``, ``urst`` can be ``None``.  In this case, the missing input is
    constructed from the others in such a way that ``urst(E) == F`` holds,
    and an exception is raised if this is impossible (typically because
    ``E`` and ``F`` are not isomorphic).

    Users will not usually need to use this class directly, but instead use
    methods such as
    :meth:`~sage.schemes.elliptic_curves.ell_generic.EllipticCurve_generic.isomorphism_to`
    or
    :meth:`~sage.schemes.elliptic_curves.ell_generic.EllipticCurve_generic.isomorphisms`.

    Explicitly, the isomorphism defined by `(u,r,s,t)` maps a point `(x,y)`
    to the point

    .. MATH::

        ((x-r) / u^2, \; (y - s(x-r) - t) / u^3) .

    If the domain `E` has Weierstrass coefficients `[a_1,a_2,a_3,a_4,a_6]`,
    the codomain `F` is given by

    .. MATH::

        a_1' &= (a_1 + 2s) / u \\
        a_2' &= (a_2 - a_1s + 3r - s^2) / u^2 \\
        a_3' &= (a_3 + a_1r + 2t) / u^3 \\
        a_4' &= (a_4 + 2a_2r - a_1(rs+t) - a_3s + 3r^2 - 2st) / u^4 \\
        a_6' &= (a_6 - a_1rt + a_2r^2 - a_3t + a_4r + r^3 - t^2) / u^6 .

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.weierstrass_morphism import *
        sage: WeierstrassIsomorphism(EllipticCurve([0,1,2,3,4]), (-1,2,3,4))
        Elliptic-curve morphism:
          From: Elliptic Curve defined by y^2 + 2*y = x^3 + x^2 + 3*x + 4 over Rational Field
          To:   Elliptic Curve defined by y^2 - 6*x*y - 10*y = x^3 - 2*x^2 - 11*x - 2 over Rational Field
          Via:  (u,r,s,t) = (-1, 2, 3, 4)
        sage: E = EllipticCurve([0,1,2,3,4])
        sage: F = EllipticCurve(E.cremona_label())
        sage: WeierstrassIsomorphism(E, None, F)
        Elliptic-curve morphism:
          From: Elliptic Curve defined by y^2 + 2*y = x^3 + x^2 + 3*x + 4 over Rational Field
          To:   Elliptic Curve defined by y^2  = x^3 + x^2 + 3*x + 5 over Rational Field
          Via:  (u,r,s,t) = (1, 0, 0, -1)
        sage: w = WeierstrassIsomorphism(None, (1,0,0,-1), F)
        sage: w._domain == E
        True
    """
    def __init__(self, E=None, urst=None, F=None):
        r"""
        Constructor for the ``WeierstrassIsomorphism`` class.

        TESTS:

        Check for :issue:`33215`::

            sage: from sage.schemes.elliptic_curves.weierstrass_morphism import WeierstrassIsomorphism
            sage: E = EllipticCurve(GF(71^2), [5,5])                                    # needs sage.rings.finite_rings
            sage: iso = WeierstrassIsomorphism(E, (1,2,3,4))                            # needs sage.rings.finite_rings
            sage: ~iso  # indirect doctest                                              # needs sage.rings.finite_rings
            Elliptic-curve morphism:
              From: Elliptic Curve defined by y^2 + 6*x*y + 8*y = x^3 + 68*x^2 + 64*x + 7 over Finite Field in z2 of size 71^2
              To:   Elliptic Curve defined by y^2 = x^3 + 5*x + 5 over Finite Field in z2 of size 71^2
              Via:  (u,r,s,t) = (1, 69, 68, 2)

        Test for :issue:`33312`::

            sage: type(iso.degree())                                                    # needs sage.rings.finite_rings
            <class 'sage.rings.integer.Integer'>
        """
        from .ell_generic import EllipticCurve_generic

        if E is not None:
            if not isinstance(E, EllipticCurve_generic):
                raise ValueError("first argument must be an elliptic curve or None")
        if F is not None:
            if not isinstance(F, EllipticCurve_generic):
                raise ValueError("third argument must be an elliptic curve or None")
        if urst is not None:
            if len(urst) != 4:
                raise ValueError("second argument must be [u,r,s,t] or None")
        if len([par for par in [E, urst, F] if par is not None]) < 2:
            raise ValueError("at most 1 argument can be None")

        inps = []
        if E is not None:
            inps.append(E.base_ring())
        if F is not None:
            inps.append(F.base_ring())
        if urst is not None:
            inps += list(urst)
        base_ring = get_coercion_model().common_parent(*inps)

        if urst is not None:
            urst = Sequence(urst, base_ring)

        if F is None:  # easy case
            baseWI.__init__(self, *urst)
            F = EllipticCurve(baseWI.__call__(self, list(E.a_invariants())))

        elif E is None:  # easy case in reverse
            baseWI.__init__(self, *urst)
            inv_urst = baseWI.__invert__(self)
            E = EllipticCurve(baseWI.__call__(inv_urst, list(F.a_invariants())))

        elif urst is None:  # try to construct the morphism
            try:
                urst = next(_isomorphisms(E, F))
            except StopIteration:
                raise ValueError("elliptic curves not isomorphic")
            baseWI.__init__(self, *urst)

        else:  # none of the parameters is None:
            baseWI.__init__(self, *urst)
            if F != EllipticCurve(baseWI.__call__(self, list(E.a_invariants()))):
                raise ValueError("second argument is not an isomorphism from first argument to third argument")

        self._mpoly_ring = PolynomialRing(base_ring, ['x','y'])
        self._poly_ring = PolynomialRing(base_ring, ['x'])
        self._xyfield = self._mpoly_ring.fraction_field()
        self._xfield = self._poly_ring.fraction_field()

        self._domain = E
        self._codomain = F
        self._degree = Integer(1)
        EllipticCurveHom.__init__(self, self._domain, self._codomain)

    @staticmethod
    def _comparison_impl(left, right, op):
        r"""
        Compare an isomorphism to another elliptic-curve morphism.

        Called by :meth:`EllipticCurveHom._richcmp_`.

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.weierstrass_morphism import *
            sage: E = EllipticCurve('389a1')
            sage: F = E.change_weierstrass_model(1,2,3,4)
            sage: w1 = E.isomorphism_to(F)
            sage: w1 == w1
            True
            sage: w2 = F.automorphisms()[1] * w1
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
        if not isinstance(left, WeierstrassIsomorphism) or not isinstance(right, WeierstrassIsomorphism):
            return NotImplemented

        lx = left._domain
        rx = right._domain
        if lx != rx:
            return richcmp_not_equal(lx, rx, op)

        lx = left._codomain
        rx = right._codomain
        if lx != rx:
            return richcmp_not_equal(lx, rx, op)

        if op in (op_EQ, op_NE):
            return richcmp(left.tuple(), right.tuple(), op)

        # This makes sure that the identity and negation morphisms
        # come first in a sorted list of WeierstrassIsomorphisms.
        # More generally, we're making sure that a morphism and its
        # negative appear next to each other, and that those pairs
        # of isomorphisms satisfying u=+-1 come first.
        def _sorting_key(iso):
            v, w = iso.tuple(), (-iso).tuple()
            i = 0 if (1,0,0,0) in (v,w) else 1
            j = 0 if v[0] == 1 else 1 if w[0] == 1 else 2
            return (i,) + min(v,w) + (j,) + v

        return richcmp(_sorting_key(left), _sorting_key(right), op)

    def _eval(self, P):
        r"""
        Less strict evaluation method for internal use.

        In particular, this can be used to evaluate ``self`` at a
        point defined over an extension field.

        INPUT: a sequence of 3 coordinates defining a point on ``self``

        OUTPUT: the result of evaluating ``self`` at the given point

        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: from sage.schemes.elliptic_curves.weierstrass_morphism import WeierstrassIsomorphism
            sage: E = EllipticCurve([i, 0]); E
            Elliptic Curve defined by y^2 = x^3 + I*x
             over Number Field in I with defining polynomial x^2 + 1 with I = 1*I
            sage: iso = WeierstrassIsomorphism(E, (i,1,2,3))
            sage: P = E.change_ring(QQbar).lift_x(QQbar.random_element())
            sage: Q = iso._eval(P)
            sage: Q.curve()
            Elliptic Curve defined by y^2 + (-4*I)*x*y + 6*I*y = x^3 + x^2 + (I-9)*x + (-I+8)
             over Algebraic Field
            sage: y = next(filter(bool, iter(QQbar.random_element, None)))  # sample until nonzero
            sage: iso._eval((0, y, 0)) == 0
            True
        """
        if self._domain.defining_polynomial()(*P):
            raise ValueError(f'{P} not on {self._domain}')
        k = Sequence(P).universe()

        Q = baseWI.__call__(self, P)
        return self._codomain.base_extend(k).point(Q)

    def _call_(self, P):
        r"""
        Call function for WeierstrassIsomorphism class.

        INPUT:

        - ``P`` -- Point; a point on the domain curve

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

        TESTS:

        Check that copying the order over works::

            sage: # needs sage.rings.finite_rings
            sage: E = EllipticCurve(GF(431^2), [1,0])
            sage: i = next(a for a in E.automorphisms() if a^2 == -a^24)
            sage: P,_ = E.gens()
            sage: P._order
            432
            sage: i(P)._order
            432
            sage: E(i(P))._order
            432

        Check that the isomorphism cannot be evaluated on points outside
        its domain (see :issue:`35799`)::

            sage: # needs sage.rings.finite_rings
            sage: E = EllipticCurve(GF(101), [1,1])
            sage: f = E.automorphisms()[0]
            sage: EE = EllipticCurve(GF(101), [5,5])
            sage: P = EE.lift_x(2)
            sage: P in f.domain()
            False
            sage: f(P)
            Traceback (most recent call last):
            ...
            TypeError: (2 : 15 : 1) fails to convert into the map's
            domain Elliptic Curve defined by y^2 = x^3 + x + 1 over
            Finite Field of size 101, but a `pushforward` method is
            not properly implemented
        """
        if P[2] == 0:
            return self._codomain(0)
        res = baseWI.__call__(self, tuple(P._coords))
        Q = self._codomain.point(res, check=False)
        if hasattr(P, '_order'):
            Q._order = P._order
        return Q

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

            sage: E1 = EllipticCurve([1,0])
            sage: phi = E1.isogeny(E1(0,0))
            sage: E2 = phi.codomain()
            sage: psi = E2.isogeny(E2(0,0))
            sage: w1._composition_impl(psi, phi)
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

    def rational_maps(self):
        """
        Return the pair of rational maps defining this isomorphism.

        EXAMPLES::

            sage: E1 = EllipticCurve([11,22,33,44,55])
            sage: E2 = EllipticCurve_from_j(E1.j_invariant())
            sage: iso = E1.isomorphism_to(E2); iso
            Elliptic-curve morphism:
              From: Elliptic Curve defined by y^2 + 11*x*y + 33*y = x^3 + 22*x^2 + 44*x + 55
                    over Rational Field
              To:   Elliptic Curve defined by y^2 + x*y = x^3 + x^2 - 684*x + 6681
                    over Rational Field
              Via:  (u,r,s,t) = (1, -17, -5, 77)
            sage: iso.rational_maps()
            (x + 17, 5*x + y + 8)
            sage: f = E2.defining_polynomial()(*iso.rational_maps(), 1)
            sage: I = E1.defining_ideal()
            sage: x,y,z = I.ring().gens()
            sage: f in I + Ideal(z-1)
            True

        ::

            sage: # needs sage.rings.finite_rings
            sage: E = EllipticCurve(GF(65537), [1,1,1,1,1])
            sage: w = E.isomorphism_to(E.short_weierstrass_model())
            sage: f,g = w.rational_maps()
            sage: P = E.random_point()
            sage: w(P).xy() == (f(P.xy()), g(P.xy()))
            True

        TESTS:

        Check for :issue:`34811`::

            sage: iso.rational_maps()[0].parent()
            Fraction Field of Multivariate Polynomial Ring in x, y over Rational Field
            sage: iso.rational_maps()[1].parent()
            Fraction Field of Multivariate Polynomial Ring in x, y over Rational Field
        """
        return tuple(baseWI.__call__(self, self._xyfield.gens()))

    def x_rational_map(self):
        """
        Return the `x`-coordinate rational map of this isomorphism.

        EXAMPLES::

            sage: E1 = EllipticCurve([11,22,33,44,55])
            sage: E2 = EllipticCurve_from_j(E1.j_invariant())
            sage: iso = E1.isomorphism_to(E2); iso
            Elliptic-curve morphism:
              From: Elliptic Curve defined by y^2 + 11*x*y + 33*y = x^3 + 22*x^2 + 44*x + 55
                    over Rational Field
              To:   Elliptic Curve defined by y^2 + x*y = x^3 + x^2 - 684*x + 6681
                    over Rational Field
              Via:  (u,r,s,t) = (1, -17, -5, 77)
            sage: iso.x_rational_map()
            x + 17
            sage: iso.x_rational_map() == iso.rational_maps()[0]
            True

        TESTS:

        Check for :issue:`34811`::

            sage: iso.x_rational_map().parent()
            Fraction Field of Univariate Polynomial Ring in x over Rational Field
        """
        x, = self._xfield.gens()
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
            Isogeny of degree 1
             from Elliptic Curve defined by y^2 + 11*x*y + 33*y = x^3 + 22*x^2 + 44*x + 55
                  over Rational Field
               to Elliptic Curve defined by y^2 + x*y = x^3 + x^2 - 684*x + 6681
                  over Rational Field
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
            sage: E = EllipticCurve(QuadraticField(-3), [0,1])                          # needs sage.rings.number_field
            sage: w = WeierstrassIsomorphism(E, (CyclotomicField(3).gen(),0,0,0))       # needs sage.rings.number_field
            sage: (w.dual() * w).rational_maps()                                        # needs sage.rings.number_field
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
              To:   Elliptic Curve defined by y^2 + 17/6*x*y + 49/13068*y = x^3 - 769/396*x^2 - 3397/862488*x + 44863/7513995456
                    over Rational Field
              Via:  (u,r,s,t) = (-66, 77, -99, -979)
            sage: -(-w) == w
            True

        ::

            sage: # needs sage.rings.number_field
            sage: from sage.schemes.elliptic_curves.weierstrass_morphism import WeierstrassIsomorphism
            sage: K.<a> = QuadraticField(-3)
            sage: E = EllipticCurve(K, [0,1])
            sage: w = WeierstrassIsomorphism(E, (CyclotomicField(3).gen(),0,0,0))
            sage: w.tuple()
            (1/2*a - 1/2, 0, 0, 0)
            sage: (-w).tuple()
            (-1/2*a + 1/2, 0, 0, 0)
            sage: (-w)^3 == -(w^3)
            True

        ::

            sage: from sage.schemes.elliptic_curves.weierstrass_morphism import WeierstrassIsomorphism, identity_morphism
            sage: E = EllipticCurve(QuadraticField(-1), [1,0])                          # needs sage.rings.number_field
            sage: t = WeierstrassIsomorphism(E, (i,0,0,0))                              # needs sage.rings.number_field
            sage: -t^2 == identity_morphism(E)                                          # needs sage.rings.number_field
            True
        """
        a1,_,a3,_,_ = self._domain.a_invariants()
        w = baseWI(-1, 0, -a1, -a3)
        urst = baseWI.__mul__(self, w).tuple()
        return WeierstrassIsomorphism(self._domain, urst, self._codomain)

    def scaling_factor(self):
        r"""
        Return the Weierstrass scaling factor associated to this
        Weierstrass isomorphism.

        The scaling factor is the constant `u` (in the base field)
        such that `\varphi^* \omega_2 = u \omega_1`, where
        `\varphi: E_1\to E_2` is this isomorphism and `\omega_i` are
        the standard Weierstrass differentials on `E_i` defined by
        `\mathrm dx/(2y+a_1x+a_3)`.

        EXAMPLES::

            sage: E = EllipticCurve(QQbar, [0,1])                                       # needs sage.rings.number_field
            sage: all(f.scaling_factor() == f.formal()[1] for f in E.automorphisms())   # needs sage.rings.number_field
            True

        ALGORITHM: The scaling factor equals the `u` component of
        the tuple `(u,r,s,t)` defining the isomorphism.
        """
        return self.u

    def inseparable_degree(self):
        r"""
        Return the inseparable degree of this Weierstrass isomorphism.

        For isomorphisms, this method always returns one.

        TESTS::

            sage: from sage.schemes.elliptic_curves.weierstrass_morphism import WeierstrassIsomorphism
            sage: WeierstrassIsomorphism.inseparable_degree(None)
            1
        """
        return Integer(1)

    def is_identity(self):
        r"""
        Check if this Weierstrass isomorphism is the identity.

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.weierstrass_morphism import WeierstrassIsomorphism
            sage: p = 97
            sage: Fp = GF(p)
            sage: E = EllipticCurve(Fp, [1, 28])
            sage: ws = WeierstrassIsomorphism(E, None, E)
            sage: ws.is_identity()
            False

        ::

            sage: from sage.schemes.elliptic_curves.weierstrass_morphism import WeierstrassIsomorphism
            sage: p = 97
            sage: Fp = GF(p)
            sage: E = EllipticCurve(Fp, [1, 28])
            sage: ws = WeierstrassIsomorphism(E, (1, 0, 0, 0), None)
            sage: ws.is_identity()
            True
        """
        return self.tuple() == (1, 0, 0, 0)

    def order(self):
        r"""
        Compute the order of this Weierstrass isomorphism if it is an automorphism.

        A :exc:`ValueError` is raised if the domain is not equal to the codomain.

        A :exc:`NotImplementedError` is raised if the order of the automorphism is not 1, 2, 3, 4 or 6.

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.weierstrass_morphism import *
            sage: p = 97
            sage: Fp = GF(p)
            sage: E = EllipticCurve(Fp, [1, 28])
            sage: ws = WeierstrassIsomorphism(E, None, E)
            sage: ws.order()
            2

        TESTS::

            sage: from sage.schemes.elliptic_curves.weierstrass_morphism import *
            sage: p = 97
            sage: Fp = GF(p)
            sage: E = EllipticCurve(Fp, [1, 28])
            sage: ws = WeierstrassIsomorphism(E, None, E)
            sage: ws.order()
            2
            sage: E1 = EllipticCurve(Fp, [1, 69])
            sage: ws = E.isomorphism_to(E1)
            sage: ws.order()
            Traceback (most recent call last):
            ...
            ValueError: the domain is different from the codomain

        ::

            sage: E = EllipticCurve_from_j(Fp(0))
            sage: ws = WeierstrassIsomorphism(E, (Fp(36), 0, 0, 0), None)
            sage: ws.order()
            6
            sage: ws2 = ws*ws
            sage: ws2.order()
            3
            sage: F2_bar = GF(2).algebraic_closure()
            sage: E = EllipticCurve_from_j(F2_bar(0))
            sage: ws = WeierstrassIsomorphism(E, None, E)
            sage: ws.order()
            3
        """
        # Check if it is an actual endomorphism
        if self._domain != self._codomain:
            raise ValueError("the domain is different from the codomain")

        if self.is_identity():
            return Integer(1)

        ws2 = WeierstrassIsomorphism._composition_impl(self, self)
        if ws2.is_identity():
            return Integer(2)

        ws3 = WeierstrassIsomorphism._composition_impl(self, ws2)
        if ws3.is_identity():
            return Integer(3)

        ws4 = WeierstrassIsomorphism._composition_impl(ws2, ws2)
        if ws4.is_identity():
            return Integer(4)

        ws6 = WeierstrassIsomorphism._composition_impl(ws2, ws4)
        if ws6.is_identity():
            return Integer(6)

        raise NotImplementedError("the order of the endomorphism is not 1, 2, 3, 4 or 6")


def identity_morphism(E):
    r"""
    Given an elliptic curve `E`, return the identity morphism
    on `E` as a :class:`WeierstrassIsomorphism`.

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.weierstrass_morphism import identity_morphism
        sage: E = EllipticCurve([5,6,7,8,9])
        sage: id_ = identity_morphism(E)
        sage: id_.rational_maps()
        (x, y)
    """
    R = E.base_ring()
    zero = R.zero()
    return WeierstrassIsomorphism(E, (R.one(), zero, zero, zero))


def negation_morphism(E):
    r"""
    Given an elliptic curve `E`, return the negation endomorphism
    `[-1]` of `E` as a :class:`WeierstrassIsomorphism`.

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.weierstrass_morphism import negation_morphism
        sage: E = EllipticCurve([5,6,7,8,9])
        sage: neg = negation_morphism(E)
        sage: neg.rational_maps()
        (x, -5*x - y - 7)
    """
    R = E.base_ring()
    return WeierstrassIsomorphism(E, (-R.one(), R.zero(), -E.a1(), -E.a3()))
