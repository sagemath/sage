r"""
Homomorphisms of rings

We give a large number of examples of ring homomorphisms.

EXAMPLES:

Natural inclusion `\ZZ \hookrightarrow \QQ`::

    sage: H = Hom(ZZ, QQ)
    sage: phi = H([1])
    sage: phi(10)
    10
    sage: phi(3/1)
    3
    sage: phi(2/3)
    Traceback (most recent call last):
    ...
    TypeError: 2/3 fails to convert into the map's domain Integer Ring, but a `pushforward` method is not properly implemented

There is no homomorphism in the other direction::

    sage: H = Hom(QQ, ZZ)
    sage: H([1])
    Traceback (most recent call last):
    ...
    ValueError: relations do not all (canonically) map to 0 under map determined by images of generators

EXAMPLES:

Reduction to finite field::

    sage: H = Hom(ZZ, GF(9, 'a'))
    sage: phi = H([1])
    sage: phi(5)
    2
    sage: psi = H([4])
    sage: psi(5)
    2

Map from single variable polynomial ring::

    sage: R.<x> = ZZ[]
    sage: phi = R.hom([2], GF(5))
    sage: phi
    Ring morphism:
      From: Univariate Polynomial Ring in x over Integer Ring
      To:   Finite Field of size 5
      Defn: x |--> 2
    sage: phi(x + 12)
    4

Identity map on the real numbers::

    sage: f = RR.hom([RR(1)]); f
    Ring endomorphism of Real Field with 53 bits of precision
      Defn: 1.00000000000000 |--> 1.00000000000000
    sage: f(2.5)
    2.50000000000000
    sage: f = RR.hom( [2.0] )
    Traceback (most recent call last):
    ...
    ValueError: relations do not all (canonically) map to 0 under map determined by images of generators

Homomorphism from one precision of field to another.

From smaller to bigger doesn't make sense::

    sage: R200 = RealField(200)
    sage: f = RR.hom( R200 )
    Traceback (most recent call last):
    ...
    TypeError: natural coercion morphism from Real Field with 53 bits of precision to Real Field with 200 bits of precision not defined

From bigger to small does::

    sage: f = RR.hom( RealField(15) )
    sage: f(2.5)
    2.500
    sage: f(RR.pi())
    3.142

Inclusion map from the reals to the complexes::

    sage: i = RR.hom([CC(1)]); i
    Ring morphism:
      From: Real Field with 53 bits of precision
      To:   Complex Field with 53 bits of precision
      Defn: 1.00000000000000 |--> 1.00000000000000
    sage: i(RR('3.1'))
    3.10000000000000

A map from a multivariate polynomial ring to itself::

    sage: R.<x,y,z> = PolynomialRing(QQ,3)
    sage: phi = R.hom([y,z,x^2]); phi
    Ring endomorphism of Multivariate Polynomial Ring in x, y, z over Rational Field
      Defn: x |--> y
            y |--> z
            z |--> x^2
    sage: phi(x+y+z)
    x^2 + y + z

An endomorphism of a quotient of a multi-variate polynomial ring::

    sage: R.<x,y> = PolynomialRing(QQ)
    sage: S.<a,b> = quo(R, ideal(1 + y^2))
    sage: phi = S.hom([a^2, -b])
    sage: phi
    Ring endomorphism of Quotient of Multivariate Polynomial Ring in x, y over Rational Field by the ideal (y^2 + 1)
      Defn: a |--> a^2
            b |--> -b
    sage: phi(b)
    -b
    sage: phi(a^2 + b^2)
    a^4 - 1

The reduction map from the integers to the integers modulo 8, viewed as
a quotient ring::

    sage: R = ZZ.quo(8*ZZ)
    sage: pi = R.cover()
    sage: pi
    Ring morphism:
      From: Integer Ring
      To:   Ring of integers modulo 8
      Defn: Natural quotient map
    sage: pi.domain()
    Integer Ring
    sage: pi.codomain()
    Ring of integers modulo 8
    sage: pi(10)
    2
    sage: pi.lift()
    Set-theoretic ring morphism:
      From: Ring of integers modulo 8
      To:   Integer Ring
      Defn: Choice of lifting map
    sage: pi.lift(13)
    5

Inclusion of ``GF(2)`` into ``GF(4,'a')``::

    sage: k = GF(2)
    sage: i = k.hom(GF(4, 'a'))
    sage: i
    Ring morphism:
      From: Finite Field of size 2
      To:   Finite Field in a of size 2^2
      Defn: 1 |--> 1
    sage: i(0)
    0
    sage: a = i(1); a.parent()
    Finite Field in a of size 2^2

We next compose the inclusion with reduction from the integers to
``GF(2)``::

    sage: pi = ZZ.hom(k)
    sage: pi
    Natural morphism:
      From: Integer Ring
      To:   Finite Field of size 2
    sage: f = i * pi
    sage: f
    Composite map:
      From: Integer Ring
      To:   Finite Field in a of size 2^2
      Defn:   Natural morphism:
              From: Integer Ring
              To:   Finite Field of size 2
            then
              Ring morphism:
              From: Finite Field of size 2
              To:   Finite Field in a of size 2^2
              Defn: 1 |--> 1
    sage: a = f(5); a
    1
    sage: a.parent()
    Finite Field in a of size 2^2

Inclusion from `\QQ` to the 3-adic field::

    sage: phi = QQ.hom(Qp(3, print_mode = 'series'))
    sage: phi
    Ring morphism:
      From: Rational Field
      To:   3-adic Field with capped relative precision 20
    sage: phi.codomain()
    3-adic Field with capped relative precision 20
    sage: phi(394)
    1 + 2*3 + 3^2 + 2*3^3 + 3^4 + 3^5 + O(3^20)

An automorphism of a quotient of a univariate polynomial ring::

    sage: R.<x> = PolynomialRing(QQ)
    sage: S.<sqrt2> = R.quo(x^2-2)
    sage: sqrt2^2
    2
    sage: (3+sqrt2)^10
    993054*sqrt2 + 1404491
    sage: c = S.hom([-sqrt2])
    sage: c(1+sqrt2)
    -sqrt2 + 1

Note that Sage verifies that the morphism is valid::

    sage: (1 - sqrt2)^2
    -2*sqrt2 + 3
    sage: c = S.hom([1-sqrt2])    # this is not valid
    Traceback (most recent call last):
    ...
    ValueError: relations do not all (canonically) map to 0 under map determined by images of generators

Endomorphism of power series ring::

    sage: R.<t> = PowerSeriesRing(QQ); R
    Power Series Ring in t over Rational Field
    sage: f = R.hom([t^2]); f
    Ring endomorphism of Power Series Ring in t over Rational Field
      Defn: t |--> t^2
    sage: R.set_default_prec(10)
    sage: s = 1/(1 + t); s
    1 - t + t^2 - t^3 + t^4 - t^5 + t^6 - t^7 + t^8 - t^9 + O(t^10)
    sage: f(s)
    1 - t^2 + t^4 - t^6 + t^8 - t^10 + t^12 - t^14 + t^16 - t^18 + O(t^20)

Frobenius on a power series ring over a finite field::

    sage: R.<t> = PowerSeriesRing(GF(5))
    sage: f = R.hom([t^5]); f
    Ring endomorphism of Power Series Ring in t over Finite Field of size 5
      Defn: t |--> t^5
    sage: a = 2 + t + 3*t^2 + 4*t^3 + O(t^4)
    sage: b = 1 + t + 2*t^2 + t^3 + O(t^5)
    sage: f(a)
    2 + t^5 + 3*t^10 + 4*t^15 + O(t^20)
    sage: f(b)
    1 + t^5 + 2*t^10 + t^15 + O(t^25)
    sage: f(a*b)
    2 + 3*t^5 + 3*t^10 + t^15 + O(t^20)
    sage: f(a)*f(b)
    2 + 3*t^5 + 3*t^10 + t^15 + O(t^20)

Homomorphism of Laurent series ring::

    sage: R.<t> = LaurentSeriesRing(QQ, 10)
    sage: f = R.hom([t^3 + t]); f
    Ring endomorphism of Laurent Series Ring in t over Rational Field
      Defn: t |--> t + t^3
    sage: s = 2/t^2 + 1/(1 + t); s
    2*t^-2 + 1 - t + t^2 - t^3 + t^4 - t^5 + t^6 - t^7 + t^8 - t^9 + O(t^10)
    sage: f(s)
    2*t^-2 - 3 - t + 7*t^2 - 2*t^3 - 5*t^4 - 4*t^5 + 16*t^6 - 9*t^7 + O(t^8)
    sage: f = R.hom([t^3]); f
    Ring endomorphism of Laurent Series Ring in t over Rational Field
      Defn: t |--> t^3
    sage: f(s)
    2*t^-6 + 1 - t^3 + t^6 - t^9 + t^12 - t^15 + t^18 - t^21 + t^24 - t^27 + O(t^30)

Note that the homomorphism must result in a converging Laurent
series, so the valuation of the image of the generator must be
positive::

    sage: R.hom([1/t])
    Traceback (most recent call last):
    ...
    ValueError: relations do not all (canonically) map to 0 under map determined by images of generators
    sage: R.hom([1])
    Traceback (most recent call last):
    ...
    ValueError: relations do not all (canonically) map to 0 under map determined by images of generators

Complex conjugation on cyclotomic fields::

    sage: K.<zeta7> = CyclotomicField(7)
    sage: c = K.hom([1/zeta7]); c
    Ring endomorphism of Cyclotomic Field of order 7 and degree 6
      Defn: zeta7 |--> -zeta7^5 - zeta7^4 - zeta7^3 - zeta7^2 - zeta7 - 1
    sage: a = (1+zeta7)^5; a
    zeta7^5 + 5*zeta7^4 + 10*zeta7^3 + 10*zeta7^2 + 5*zeta7 + 1
    sage: c(a)
    5*zeta7^5 + 5*zeta7^4 - 4*zeta7^2 - 5*zeta7 - 4
    sage: c(zeta7 + 1/zeta7)       # this element is obviously fixed by inversion
    -zeta7^5 - zeta7^4 - zeta7^3 - zeta7^2 - 1
    sage: zeta7 + 1/zeta7
    -zeta7^5 - zeta7^4 - zeta7^3 - zeta7^2 - 1

Embedding a number field into the reals::

    sage: R.<x> = PolynomialRing(QQ)
    sage: K.<beta> = NumberField(x^3 - 2)
    sage: alpha = RR(2)^(1/3); alpha
    1.25992104989487
    sage: i = K.hom([alpha],check=False); i
    Ring morphism:
      From: Number Field in beta with defining polynomial x^3 - 2
      To:   Real Field with 53 bits of precision
      Defn: beta |--> 1.25992104989487
    sage: i(beta)
    1.25992104989487
    sage: i(beta^3)
    2.00000000000000
    sage: i(beta^2 + 1)
    2.58740105196820

An example from Jim Carlson::

    sage: K = QQ # by the way :-)
    sage: R.<a,b,c,d> = K[]; R
    Multivariate Polynomial Ring in a, b, c, d over Rational Field
    sage: S.<u> = K[]; S
    Univariate Polynomial Ring in u over Rational Field
    sage: f = R.hom([0,0,0,u], S); f
    Ring morphism:
      From: Multivariate Polynomial Ring in a, b, c, d over Rational Field
      To:   Univariate Polynomial Ring in u over Rational Field
      Defn: a |--> 0
            b |--> 0
            c |--> 0
            d |--> u
    sage: f(a+b+c+d)
    u
    sage: f( (a+b+c+d)^2 )
    u^2

TESTS::

    sage: H = Hom(ZZ, QQ)
    sage: H == loads(dumps(H))
    True

::

    sage: K.<zeta7> = CyclotomicField(7)
    sage: c = K.hom([1/zeta7])
    sage: c == loads(dumps(c))
    True

::

    sage: R.<t> = PowerSeriesRing(GF(5))
    sage: f = R.hom([t^5])
    sage: f == loads(dumps(f))
    True

We define the identity map in many possible ways. These should all
compare equal::

    sage: k = GF(2)
    sage: R.<x> = k[]
    sage: F4.<a> = R.quo(x^2+x+1)
    sage: H = End(F4)

    sage: from sage.rings.morphism import *
    sage: phi1 = H.identity(); phi1
    Identity endomorphism of Univariate Quotient Polynomial Ring in a over Finite Field of size 2 with modulus x^2 + x + 1
    sage: phi2 = H([a]); phi2
    Ring endomorphism of Univariate Quotient Polynomial Ring in a over Finite Field of size 2 with modulus x^2 + x + 1
      Defn: a |--> a
    sage: phi3 = RingHomomorphism_from_base(H, R.hom([x])); phi3
    Ring endomorphism of Univariate Quotient Polynomial Ring in a over Finite Field of size 2 with modulus x^2 + x + 1
      Defn: Induced from base ring by
            Ring endomorphism of Univariate Polynomial Ring in x over Finite Field of size 2 (using GF2X)
              Defn: x |--> x
    sage: phi4 = RingHomomorphism_cover(H); phi4
    Ring endomorphism of Univariate Quotient Polynomial Ring in a over Finite Field of size 2 with modulus x^2 + x + 1
      Defn: Natural quotient map
    sage: phi5 = F4.frobenius_endomorphism() ^ 2; phi5
    Frobenius endomorphism x |--> x^(2^2) of Univariate Quotient Polynomial Ring in a over Finite Field of size 2 with modulus x^2 + x + 1
    sage: maps = [phi1, phi2, phi3, phi4, phi5]
    sage: for f in maps:
    ....:     for g in maps:
    ....:         if f != g:
    ....:             print("{} != {}".format(f, g))
"""

#*****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import print_function, absolute_import

from cpython.object cimport Py_EQ, Py_NE

from . import ideal
import sage.structure.all
from sage.structure.richcmp cimport (richcmp, rich_to_bool,
        richcmp_not_equal)


def is_RingHomomorphism(phi):
    """
    Return ``True`` if ``phi`` is of type :class:`RingHomomorphism`.

    EXAMPLES::

        sage: f = Zmod(8).cover()
        sage: sage.rings.morphism.is_RingHomomorphism(f)
        doctest:warning
        ...
        DeprecationWarning: is_RingHomomorphism() should not be used anymore. Check whether the category_for() your morphism is a subcategory of Rings() instead.
        See http://trac.sagemath.org/23204 for details.
        True
        sage: sage.rings.morphism.is_RingHomomorphism(2/3)
        False
    """
    sage.misc.superseded.deprecation(23204, "is_RingHomomorphism() should not be used anymore. Check whether the category_for() your morphism is a subcategory of Rings() instead.")
    # We use the category framework to determine whether something is a ring homomorphism.
    from sage.categories.map import Map
    from sage.categories.all import Rings
    return isinstance(phi, Map) and phi.category_for().is_subcategory(Rings())


cdef class RingMap(Morphism):
    """
    Set-theoretic map between rings.

    TESTS:

    This is an abstract base class that is not directly instantiated,
    but we will do so anyway as a test::

        sage: f = sage.rings.morphism.RingMap(ZZ.Hom(ZZ))
        sage: parent(f)
        Set of Homomorphisms from Integer Ring to Integer Ring
        sage: type(f)
        <type 'sage.rings.morphism.RingMap'>
    """
    def _repr_type(self):
        """
        TESTS::

            sage: f = sage.rings.morphism.RingMap(ZZ.Hom(ZZ))
            sage: type(f)
            <type 'sage.rings.morphism.RingMap'>
            sage: f._repr_type()
            'Set-theoretic ring'
            sage: f
            Set-theoretic ring endomorphism of Integer Ring
        """
        return "Set-theoretic ring"


cdef class RingMap_lift(RingMap):
    r"""
    Given rings `R` and `S` such that for any
    `x \in R` the function ``x.lift()`` is an
    element that naturally coerces to `S`, this returns the
    set-theoretic ring map `R \to S` sending `x` to
    ``x.lift()``.

    EXAMPLES::

        sage: R.<x,y> = QQ[]
        sage: S.<xbar,ybar> = R.quo( (x^2 + y^2, y) )
        sage: S.lift()
        Set-theoretic ring morphism:
          From: Quotient of Multivariate Polynomial Ring in x, y over Rational Field by the ideal (x^2 + y^2, y)
          To:   Multivariate Polynomial Ring in x, y over Rational Field
          Defn: Choice of lifting map
        sage: S.lift() == 0
        False

    Since :trac:`11068`, it is possible to create
    quotient rings of non-commutative rings by two-sided
    ideals. It was needed to modify :class:`RingMap_lift`
    so that rings can be accepted that are no instances
    of :class:`sage.rings.ring.Ring`, as in the following
    example::

        sage: MS = MatrixSpace(GF(5),2,2)
        sage: I = MS*[MS.0*MS.1,MS.2+MS.3]*MS
        sage: Q = MS.quo(I)
        sage: Q.0*Q.1   # indirect doctest
        [0 1]
        [0 0]
    """
    def __init__(self, R, S):
        """
        Create a lifting ring map.

        EXAMPLES::

            sage: f = Zmod(8).lift()          # indirect doctest
            sage: f(3)
            3
            sage: type(f(3))
            <type 'sage.rings.integer.Integer'>
            sage: type(f)
            <type 'sage.rings.morphism.RingMap_lift'>

        An invalid example::

            sage: GF9.<one, a> = GaussianIntegers().quotient(3)
            sage: from sage.rings.morphism import RingMap_lift
            sage: RingMap_lift(GF9, ZZ)
            Traceback (most recent call last):
            ...
            TypeError: no canonical coercion from Number Field in I with defining polynomial x^2 + 1 with I = 1*I to Integer Ring
        """
        self.S = <Parent?>S
        x = <Element?>R(0).lift()
        f = self.S.coerce_map_from(x._parent)
        if f is None:
            raise TypeError(f"no canonical coercion from {x._parent} to {S}")
        self.to_S = f

        from sage.categories.sets_cat import Sets
        H = R.Hom(S, Sets())
        RingMap.__init__(self, H)

    cdef _update_slots(self, dict _slots):
        """
        Helper for copying and pickling.

        EXAMPLES::

            sage: f = Zmod(8).lift()
            sage: g = copy(f)    # indirect doctest
            sage: g(3) == f(3)
            True
            sage: f == g
            True
            sage: f is g
            False
        """
        self.S = _slots['S']
        self.to_S = _slots['to_S']
        Morphism._update_slots(self, _slots)

    cdef dict _extra_slots(self):
        """
        Helper for copying and pickling.

        EXAMPLES::

            sage: f = Zmod(8).lift()
            sage: g = copy(f)    # indirect doctest
            sage: g(3) == f(3)
            True
        """
        slots = Morphism._extra_slots(self)
        slots['S'] = self.S
        slots['to_S'] = self.to_S
        return slots

    cpdef _richcmp_(self, other, int op):
        """
        Compare a ring lifting maps ``self`` to ``other``.

        Ring lifting maps never compare equal to any other data type.
        If ``other`` is a ring lifting maps, the parents of ``self`` and
        ``other`` are compared.

        EXAMPLES::

            sage: f = Zmod(8).lift()
            sage: g = Zmod(10).lift()
            sage: f == f
            True
            sage: f == g
            False

        Verify that :trac:`5758` has been fixed::

            sage: Zmod(8).lift() == 1
            False
        """
        if not isinstance(other, RingMap_lift):
            # Generic comparison
            return RingMap._richcmp_(self, other, op)
        # Two lifting maps with the same parent must be equal
        return rich_to_bool(op, 0)

    def __hash__(self):
        """
        Return the hash of this morphism.

        TESTS::

            sage: f = Zmod(8).lift()
            sage: type(f)
            <type 'sage.rings.morphism.RingMap_lift'>
            sage: hash(f) == hash(f)
            True
            sage: {f: 1}[f]
            1
            sage: g = Zmod(10).lift()
            sage: hash(f) == hash(g)
            False
        """
        return hash((self.domain(), self.codomain()))

    def _repr_defn(self):
        """
        Used in printing out lifting maps.

        EXAMPLES::

            sage: f = Zmod(8).lift()
            sage: f._repr_defn()
            'Choice of lifting map'
            sage: f
            Set-theoretic ring morphism:
              From: Ring of integers modulo 8
              To:   Integer Ring
              Defn: Choice of lifting map
        """
        return "Choice of lifting map"

    cpdef Element _call_(self, x):
        """
        Evaluate this function at ``x``.

        EXAMPLES::

            sage: f = Zmod(8).lift()
            sage: type(f)
            <type 'sage.rings.morphism.RingMap_lift'>
            sage: f(-1)                       # indirect doctest
            7
            sage: type(f(-1))
            <type 'sage.rings.integer.Integer'>
        """
        return self.to_S(x.lift())


cdef class RingHomomorphism(RingMap):
    """
    Homomorphism of rings.
    """
    def __init__(self, parent):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: f = ZZ.hom(Zp(3)); f
            Ring morphism:
              From: Integer Ring
              To:   3-adic Ring with capped relative precision 20

        TESTS::

            sage: isinstance(f, sage.rings.morphism.RingHomomorphism)
            True

        """
        from .homset import RingHomset_generic
        if not isinstance(parent, RingHomset_generic):
            raise TypeError("parent must be a ring homset")
        RingMap.__init__(self, parent)

    def _repr_type(self):
        """
        Used internally in printing this morphism.

        TESTS::

            sage: ZZ.hom(Zp(3))._repr_type()
            'Ring'

        """
        return "Ring"

    def _set_lift(self, lift):
        r"""
        Used internally to define a lifting homomorphism associated to
        this homomorphism, which goes in the other direction.  I.e.,
        if ``self`` is from `R` to `S`, then the lift must be a set-theoretic
        map from `S` to `R` such that ``self(lift(x)) == x``.

        INPUT:

        - ``lift`` -- a ring map

        OUTPUT:

        Changes the state of ``self``.

        EXAMPLES::

            sage: R = ZZ.quo(3*ZZ)
            sage: pi = R.cover() # indirect doctest
            sage: pi.lift()
            Set-theoretic ring morphism:
              From: Ring of integers modulo 3
              To:   Integer Ring
              Defn: Choice of lifting map

        """
        if lift.domain() != self.codomain():
            raise TypeError("lift must have correct domain")
        if lift.codomain() != self.domain():
            raise TypeError("lift must have correct codomain")
        self._lift = lift

    cdef _update_slots(self, dict _slots):
        """
        Helper for copying and pickling.

        EXAMPLES::

            sage: f = ZZ.hom(Zmod(6))
            sage: g = copy(f)    # indirect doctest
            sage: g == f
            True
            sage: g is f
            False
            sage: g(7)
            1
        """
        if '_lift' in _slots:
            self._lift = _slots['_lift']
        Morphism._update_slots(self, _slots)

    cdef dict _extra_slots(self):
        """
        Helper for copying and pickling.

        EXAMPLES::

            sage: f = ZZ.hom(Zmod(6))
            sage: g = copy(f)    # indirect doctest
            sage: g == f
            True
            sage: g is f
            False
            sage: g(7)
            1
        """
        slots = Morphism._extra_slots(self)
        try:
            slots['_lift'] = self._lift
        except AttributeError:
            pass
        return slots

    def _composition_(self, right, homset):
        """
        If ``homset`` is a homset of rings and ``right`` is a
        ring homomorphism given by the images of generators,
        (indirectly in the case of homomorphisms from relative
        number fields), the composition with ``self`` will be
        of the appropriate type.

        Otherwise, a formal composite map is returned.

        EXAMPLES::

            sage: R.<x,y> = QQ[]
            sage: S.<a,b> = QQ[]
            sage: f = R.hom([a+b,a-b])
            sage: g = S.hom(Frac(S))
            sage: g*f # indirect doctest
            Composite map:
              From: Multivariate Polynomial Ring in x, y over Rational Field
              To:   Fraction Field of Multivariate Polynomial Ring in a, b over Rational Field
              Defn:   Ring morphism:
                      From: Multivariate Polynomial Ring in x, y over Rational Field
                      To:   Multivariate Polynomial Ring in a, b over Rational Field
                      Defn: x |--> a + b
                            y |--> a - b
                    then
                      Coercion map:
                      From: Multivariate Polynomial Ring in a, b over Rational Field
                      To:   Fraction Field of Multivariate Polynomial Ring in a, b over Rational Field

        When ``right`` is defined by the images of generators, the
        result has the type of a homomorphism between its domain and
        codomain::

            sage: C = CyclotomicField(24)
            sage: f = End(C)[1]
            sage: type(f*f) == type(f)
            True

        An example where the domain of ``right`` is a relative number field::

            sage: PQ.<X> = QQ[]
            sage: K.<a, b> = NumberField([X^2 - 2, X^2 - 3])
            sage: e, u, v, w = End(K)
            sage: u*v
            Relative number field endomorphism of Number Field in a with defining polynomial X^2 - 2 over its base field
              Defn: a |--> -a
                    b |--> b

        An example where ``right`` is not a ring homomorphism::

            sage: from sage.categories.morphism import SetMorphism
            sage: h = SetMorphism(Hom(R,S,Rings()), lambda p: p[0])
            sage: g*h
            Composite map:
              From: Multivariate Polynomial Ring in x, y over Rational Field
              To:   Fraction Field of Multivariate Polynomial Ring in a, b over Rational Field
              Defn:   Generic morphism:
                      From: Multivariate Polynomial Ring in x, y over Rational Field
                      To:   Multivariate Polynomial Ring in a, b over Rational Field
                    then
                      Coercion map:
                      From: Multivariate Polynomial Ring in a, b over Rational Field
                      To:   Fraction Field of Multivariate Polynomial Ring in a, b over Rational Field

        AUTHORS:

        - Simon King (2010-05)
        - Francis Clarke (2011-02)
        """
        from sage.categories.morphism import IdentityMorphism
        from sage.categories.rings import Rings
        if isinstance(right, IdentityMorphism):
            return self
        if homset.homset_category().is_subcategory(Rings()):
            if isinstance(right, RingHomomorphism_im_gens):
                try:
                    return homset([self(g) for g in right.im_gens()], False)
                except ValueError:
                    pass
            from sage.rings.number_field.morphism import RelativeNumberFieldHomomorphism_from_abs
            if isinstance(right, RelativeNumberFieldHomomorphism_from_abs):
                try:
                    return homset(self*right.abs_hom())
                except ValueError:
                    pass
        return sage.categories.map.Map._composition_(self, right, homset)

    def pushforward(self, I):
        """
        Returns the pushforward of the ideal `I` under this ring
        homomorphism.

        EXAMPLES::

            sage: R.<x,y> = QQ[]; S.<xx,yy> = R.quo([x^2,y^2]);  f = S.cover()
            sage: f.pushforward(R.ideal([x,3*x+x*y+y^2]))
            Ideal (xx, xx*yy + 3*xx) of Quotient of Multivariate Polynomial Ring in x, y over Rational Field by the ideal (x^2, y^2)
        """
        if not ideal.is_Ideal(I):
            raise TypeError("I must be an ideal")
        R = self.codomain()
        return R.ideal([self(y) for y in I.gens()])

    def inverse_image(self, I):
        """
        Return the inverse image of the ideal `I` under this ring
        homomorphism.

        EXAMPLES:

        This is not implemented in any generality yet::

            sage: f = ZZ.hom(Zp(2))
            sage: f.inverse_image(ZZ.ideal(2))
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def lift(self, x=None):
        """
        Return a lifting homomorphism associated to this homomorphism, if
        it has been defined.

        If ``x`` is not ``None``, return the value of the lift morphism on
        ``x``.

        EXAMPLES::

            sage: R.<x,y> = QQ[]
            sage: f = R.hom([x,x])
            sage: f(x+y)
            2*x
            sage: f.lift()
            Traceback (most recent call last):
            ...
            ValueError: no lift map defined
            sage: g = R.hom(R)
            sage: f._set_lift(g)
            sage: f.lift() == g
            True
            sage: f.lift(x)
            x
        """
        if self._lift is None:
            raise ValueError("no lift map defined")
        if x is None:
            return self._lift
        return self._lift(x)


cdef class RingHomomorphism_coercion(RingHomomorphism):
    r"""
    A ring homomorphism that is a coercion.

    .. WARNING:;

        This class is obsolete. Set the category of your morphism to a
        subcategory of ``Rings`` instead.

    TESTS:

        sage: from sage.rings.morphism import RingHomomorphism_coercion
        sage: parent = Hom(ZZ,ZZ)
        sage: f = parent.__make_element_class__(RingHomomorphism_coercion)(parent)
        doctest:warning
        ...
        DeprecationWarning: Set the category of your morphism to a subcategory of Rings instead.
        See http://trac.sagemath.org/23204 for details.
        sage: TestSuite(f).run()

    """
    def __init__(self, parent, check = True):
        r"""
        TESTS:

            sage: from sage.rings.morphism import RingHomomorphism_coercion
            sage: parent = Hom(ZZ,ZZ)
            sage: f = parent.__make_element_class__(RingHomomorphism_coercion)(parent) # py2
            sage: f = parent.__make_element_class__(RingHomomorphism_coercion)(parent) # py3
            doctest:warning
            ...
            DeprecationWarning: Set the category of your morphism to a subcategory of Rings instead.
            See http://trac.sagemath.org/23204 for details.
            sage: isinstance(f, RingHomomorphism_coercion)
            True

        """
        sage.misc.superseded.deprecation(23204, "Set the category of your morphism to a subcategory of Rings instead.")

        RingHomomorphism.__init__(self, parent)
        # putting in check allows us to define subclasses of RingHomomorphism_coercion that implement _coerce_map_from
        if check and not self.codomain().has_coerce_map_from(self.domain()):
            raise TypeError("Natural coercion morphism from %s to %s not defined."%(self.domain(), self.codomain()))

    def _repr_type(self):
        """
        Used internally when printing this.

        EXAMPLES::

            sage: from sage.rings.morphism import RingHomomorphism_coercion
            sage: parent = Hom(ZZ,ZZ)
            sage: f = parent.__make_element_class__(RingHomomorphism_coercion)(parent)
            sage: f._repr_type()
            'Ring Coercion'

        """
        return "Ring Coercion"

    cpdef _richcmp_(self, other, int op):
        """
        Compare a ring coercion morphism ``self`` to ``other``.

        Ring coercion morphisms never compare equal to any other data type. If
        other is a ring coercion morphism, the parents of ``self`` and
        ``other`` are compared.

        EXAMPLES::

            sage: from sage.rings.morphism import RingHomomorphism_coercion
            sage: parent = Hom(ZZ,ZZ)
            sage: f = parent.__make_element_class__(RingHomomorphism_coercion)(parent)
            sage: f == f
            True
            sage: f != f
            False
        """
        if not isinstance(other, RingHomomorphism_coercion):
            # Generic comparison
            return RingMap._richcmp_(self, other, op)
        # Two coercion maps with the same parent must be equal
        return rich_to_bool(op, 0)

    def __hash__(self):
        """
        Return the hash of this morphism.

        TESTS::

            sage: from sage.rings.morphism import RingHomomorphism_coercion
            sage: parent = Hom(ZZ,ZZ)
            sage: f = parent.__make_element_class__(RingHomomorphism_coercion)(parent)
            sage: g = parent.__make_element_class__(RingHomomorphism_coercion)(parent)
            sage: hash(f) == hash(g)
            True

        """
        return hash((self.domain(), self.codomain()))

    cpdef Element _call_(self, x):
        """
        Evaluate this coercion morphism at ``x``.

        EXAMPLES::

            sage: from sage.rings.morphism import RingHomomorphism_coercion
            sage: parent = Hom(ZZ,ZZ)
            sage: f = parent.__make_element_class__(RingHomomorphism_coercion)(parent)
            sage: f(0)
            0

        """
        return self.codomain().coerce(x)


cdef class RingHomomorphism_im_gens(RingHomomorphism):
    """
    A ring homomorphism determined by the images of generators.
    """
    def __init__(self, parent, im_gens, check=True):
        """
        EXAMPLES::

            sage: R.<x,y> = QQ[]
            sage: phi = R.hom([x,x+y]); phi
            Ring endomorphism of Multivariate Polynomial Ring in x, y over Rational Field
              Defn: x |--> x
                    y |--> x + y
            sage: type(phi)
            <type 'sage.rings.morphism.RingHomomorphism_im_gens'>

        Here's another example where the domain isn't free::

            sage: S.<xx,yy> = R.quotient(x - y)
            sage: phi = S.hom([xx+1,xx+1])

        Note that one has to specify valid images::

            sage: phi = S.hom([xx+1,xx-1])
            Traceback (most recent call last):
            ...
            TypeError: images do not define a valid homomorphism

        There is a check option, but it may be ignored in some cases
        -- it's purpose isn't so you can lie to Sage, but to sometimes
        speed up creation of a homomorphism::

            sage: phi = S.hom([xx+1,xx-1],check=False)
            Traceback (most recent call last):
            ...
            TypeError: images do not define a valid homomorphism
        """
        RingHomomorphism.__init__(self, parent)
        if not isinstance(im_gens, sage.structure.sequence.Sequence_generic):
            if not isinstance(im_gens, (tuple, list)):
                im_gens = [im_gens]
            im_gens = sage.structure.all.Sequence(im_gens, parent.codomain(),
                                                  check=check, immutable=True)
        if check:
            if len(im_gens) != parent.domain().ngens():
                raise ValueError("number of images must equal number of generators")
            t = parent.domain()._is_valid_homomorphism_(parent.codomain(), im_gens)
            if not t:
                raise ValueError("relations do not all (canonically) map to 0 under map determined by images of generators")
        if not im_gens.is_immutable():
            import copy
            im_gens = copy.copy(im_gens)
            im_gens.set_immutable()
        self.__im_gens = im_gens

    def im_gens(self):
        """
        Return the images of the generators of the domain.

        OUTPUT:

        - ``list`` -- a copy of the list of gens (it is safe to change this)

        EXAMPLES::

            sage: R.<x,y> = QQ[]
            sage: f = R.hom([x,x+y])
            sage: f.im_gens()
            [x, x + y]

        We verify that the returned list of images of gens is a copy,
        so changing it doesn't change ``f``::

            sage: f.im_gens()[0] = 5
            sage: f.im_gens()
            [x, x + y]
        """
        return list(self.__im_gens)

    cdef _update_slots(self, dict _slots):
        """
        Helper for copying and pickling.

        EXAMPLES::

            sage: R.<x,y> = QQ[]
            sage: f = R.hom([x,x+y])
            sage: g = copy(f)   # indirect doctest
            sage: g == f
            True
            sage: g is f
            False
            sage: g(y)
            x + y
        """
        self.__im_gens = _slots['__im_gens']
        RingHomomorphism._update_slots(self, _slots)

    cdef dict _extra_slots(self):
        """
        Helper for copying and pickling.

        EXAMPLES::

            sage: R.<x,y> = QQ[]
            sage: f = R.hom([x,x+y])
            sage: g = copy(f)   # indirect doctest
            sage: g == f
            True
            sage: g is f
            False
            sage: g(y)
            x + y
        """
        slots = RingHomomorphism._extra_slots(self)
        slots['__im_gens'] = self.__im_gens
        return slots

    cpdef _richcmp_(self, other, int op):
        r"""
        EXAMPLES:

        A single variate quotient over `\QQ`::

            sage: R.<x> = QQ[]
            sage: Q.<a> = R.quotient(x^2 + x + 1)
            sage: f1 = R.hom([a])
            sage: f2 = R.hom([a + a^2 + a + 1])
            sage: f1 == f2
            True
            sage: f1 == R.hom([a^2])
            False
            sage: f1(x^3 + x)
            a + 1
            sage: f2(x^3 + x)
            a + 1

        TESTS::

            sage: loads(dumps(f2)) == f2
            True

        ::

            sage: R.<x,y> = QQ[]; f = R.hom([x,x+y]); g = R.hom([y,x])
            sage: f == g             # indirect doctest
            False

        EXAMPLES:

        A multivariate quotient over a finite field::

            sage: R.<x,y> = GF(7)[]
            sage: Q.<a,b> = R.quotient([x^2 + x + 1, y^2 + y + 1])
            sage: f1 = R.hom([a, b])
            sage: f2 = R.hom([a + a^2 + a + 1, b + b^2 + b + 1])
            sage: f1 == f2
            True
            sage: f1 == R.hom([b,a])
            False
            sage: x^3 + x + y^2
            x^3 + y^2 + x
            sage: f1(x^3 + x + y^2)
            a - b
            sage: f2(x^3 + x + y^2)
            a - b

        TESTS::

            sage: loads(dumps(f2)) == f2
            True

        This was fixed in :trac:`24277`::

            sage: H = End(QQ)
            sage: H(1) == H.identity()
            True
        """
        if not isinstance(other, RingHomomorphism_im_gens):
            # Generic comparison
            return RingMap._richcmp_(self, other, op)
        # Check equality using the images of the generators.
        self_im = self.__im_gens
        other_im = (<RingHomomorphism_im_gens>other).__im_gens
        return richcmp(self_im, other_im, op)

    def __hash__(self):
        """
        Return the hash of this morphism.

        EXAMPLES::

            sage: R.<x> = ZZ[]
            sage: s = R.hom([x+1])
            sage: type(s)
            <type 'sage.rings.morphism.RingHomomorphism_im_gens'>
            sage: hash(s) == hash(s)
            True
            sage: {s: 1}[s]
            1
        """
        return hash(self.__im_gens)

    def _repr_defn(self):
        """
        Used in constructing string representation of ``self``.

        EXAMPLES::

            sage: R.<x,y> = QQ[]; f = R.hom([x^2,x+y])
            sage: print(f._repr_defn())
            x |--> x^2
            y |--> x + y
        """
        D = self.domain()
        ig = self.__im_gens
        return '\n'.join(['%s |--> %s'%(D.gen(i), ig[i]) for\
                       i in range(D.ngens())])

    cpdef Element _call_(self, x):
        """
        Evaluate this homomorphism at ``x``.

        EXAMPLES::

            sage: R.<x,y,z> = ZZ[]; f = R.hom([2*x,z,y])
            sage: f(x+2*y+3*z)             # indirect doctest
            2*x + 3*y + 2*z
        """
        return x._im_gens_(self.codomain(), self.im_gens())


cdef class RingHomomorphism_from_base(RingHomomorphism):
    """
    A ring homomorphism determined by a ring homomorphism of the base ring.

    AUTHOR:

    - Simon King (initial version, 2010-04-30)

    EXAMPLES:

    We define two polynomial rings and a ring homomorphism::

        sage: R.<x,y> = QQ[]
        sage: S.<z> = QQ[]
        sage: f = R.hom([2*z,3*z],S)

    Now we construct polynomial rings based on ``R`` and ``S``, and let
    ``f`` act on the coefficients::

        sage: PR.<t> = R[]
        sage: PS = S['t']
        sage: Pf = PR.hom(f,PS)
        sage: Pf
        Ring morphism:
          From: Univariate Polynomial Ring in t over Multivariate Polynomial Ring in x, y over Rational Field
          To:   Univariate Polynomial Ring in t over Univariate Polynomial Ring in z over Rational Field
          Defn: Induced from base ring by
                Ring morphism:
                  From: Multivariate Polynomial Ring in x, y over Rational Field
                  To:   Univariate Polynomial Ring in z over Rational Field
                  Defn: x |--> 2*z
                        y |--> 3*z
        sage: p = (x - 4*y + 1/13)*t^2 + (1/2*x^2 - 1/3*y^2)*t + 2*y^2 + x
        sage: Pf(p)
        (-10*z + 1/13)*t^2 - z^2*t + 18*z^2 + 2*z

    Similarly, we can construct the induced homomorphism on a matrix ring over
    our polynomial rings::

        sage: MR = MatrixSpace(R,2,2)
        sage: MS = MatrixSpace(S,2,2)
        sage: M = MR([x^2 + 1/7*x*y - y^2, - 1/2*y^2 + 2*y + 1/6, 4*x^2 - 14*x, 1/2*y^2 + 13/4*x - 2/11*y])
        sage: Mf = MR.hom(f,MS)
        sage: Mf
        Ring morphism:
          From: Full MatrixSpace of 2 by 2 dense matrices over Multivariate Polynomial Ring in x, y over Rational Field
          To:   Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in z over Rational Field
          Defn: Induced from base ring by
                Ring morphism:
                  From: Multivariate Polynomial Ring in x, y over Rational Field
                  To:   Univariate Polynomial Ring in z over Rational Field
                  Defn: x |--> 2*z
                        y |--> 3*z
        sage: Mf(M)
        [           -29/7*z^2 -9/2*z^2 + 6*z + 1/6]
        [       16*z^2 - 28*z   9/2*z^2 + 131/22*z]

    The construction of induced homomorphisms is recursive, and so we have::

        sage: MPR = MatrixSpace(PR, 2)
        sage: MPS = MatrixSpace(PS, 2)
        sage: M = MPR([(- x + y)*t^2 + 58*t - 3*x^2 + x*y, (- 1/7*x*y - 1/40*x)*t^2 + (5*x^2 + y^2)*t + 2*y, (- 1/3*y + 1)*t^2 + 1/3*x*y + y^2 + 5/2*y + 1/4, (x + 6*y + 1)*t^2])
        sage: MPf = MPR.hom(f,MPS); MPf
        Ring morphism:
          From: Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in t over Multivariate Polynomial Ring in x, y over Rational Field
          To:   Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in t over Univariate Polynomial Ring in z over Rational Field
          Defn: Induced from base ring by
                Ring morphism:
                  From: Univariate Polynomial Ring in t over Multivariate Polynomial Ring in x, y over Rational Field
                  To:   Univariate Polynomial Ring in t over Univariate Polynomial Ring in z over Rational Field
                  Defn: Induced from base ring by
                        Ring morphism:
                          From: Multivariate Polynomial Ring in x, y over Rational Field
                          To:   Univariate Polynomial Ring in z over Rational Field
                          Defn: x |--> 2*z
                                y |--> 3*z
        sage: MPf(M)
        [                    z*t^2 + 58*t - 6*z^2 (-6/7*z^2 - 1/20*z)*t^2 + 29*z^2*t + 6*z]
        [    (-z + 1)*t^2 + 11*z^2 + 15/2*z + 1/4                           (20*z + 1)*t^2]
    """
    def __init__(self, parent, underlying):
        """
        Initialize ``self``.

        TESTS::

            sage: from sage.rings.morphism import RingHomomorphism_from_base
            sage: R.<x> = ZZ[]
            sage: f = R.hom([2*x],R)
            sage: P = MatrixSpace(R,2).Hom(MatrixSpace(R,2))
            sage: g = RingHomomorphism_from_base(P,f)
            sage: g
            Ring endomorphism of Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in x over Integer Ring
              Defn: Induced from base ring by
                    Ring endomorphism of Univariate Polynomial Ring in x over Integer Ring
                      Defn: x |--> 2*x

        Note that an induced homomorphism only makes sense if domain and
        codomain are constructed in a compatible way. So, the following
        results in an error::

            sage: P = MatrixSpace(R,2).Hom(R['t'])
            sage: g = RingHomomorphism_from_base(P,f)
            Traceback (most recent call last):
            ...
            ValueError: domain (Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in x over Integer Ring) and codomain (Univariate Polynomial Ring in t over Univariate Polynomial Ring in x over Integer Ring) must have the same functorial construction over their base rings
        """
        RingHomomorphism.__init__(self, parent)
        if underlying.domain() != parent.domain().base():
            raise ValueError("The given homomorphism has to have the domain %s"%parent.domain().base())
        if underlying.codomain() != parent.codomain().base():
            raise ValueError("The given homomorphism has to have the codomain %s"%parent.codomain().base())
        if parent.domain().construction()[0] != parent.codomain().construction()[0]:
            raise ValueError(f"domain ({parent.domain()}) and codomain ({parent.codomain()}) must have the same functorial construction over their base rings")
        self.__underlying = underlying

    def underlying_map(self):
        """
        Return the underlying homomorphism of the base ring.

        EXAMPLES::

            sage: R.<x,y> = QQ[]
            sage: S.<z> = QQ[]
            sage: f = R.hom([2*z,3*z],S)
            sage: MR = MatrixSpace(R,2)
            sage: MS = MatrixSpace(S,2)
            sage: g = MR.hom(f,MS)
            sage: g.underlying_map() == f
            True
        """
        return self.__underlying

    cdef _update_slots(self, dict _slots):
        """
        Helper for copying and pickling.

        EXAMPLES::

            sage: R.<x,y> = QQ[]
            sage: S.<z> = QQ[]
            sage: f = R.hom([2*z,3*z],S)
            sage: PR.<t> = R[]
            sage: PS = S['t']
            sage: phi = PR.hom(f,PS)
            sage: type(phi)
            <type 'sage.rings.morphism.RingHomomorphism_from_base'>
            sage: psi = copy(phi); psi    # indirect doctest
            Ring morphism:
              From: Univariate Polynomial Ring in t over Multivariate Polynomial Ring in x, y over Rational Field
              To:   Univariate Polynomial Ring in t over Univariate Polynomial Ring in z over Rational Field
              Defn: Induced from base ring by
                    Ring morphism:
                      From: Multivariate Polynomial Ring in x, y over Rational Field
                      To:   Univariate Polynomial Ring in z over Rational Field
                      Defn: x |--> 2*z
                            y |--> 3*z
            sage: psi(x*t)
            2*z*t
        """
        self.__underlying = _slots['__underlying']
        RingHomomorphism._update_slots(self, _slots)

    cdef dict _extra_slots(self):
        """
        Helper for copying and pickling.

        EXAMPLES::

            sage: R.<x,y> = QQ[]
            sage: S.<z> = QQ[]
            sage: f = R.hom([2*z,3*z],S)
            sage: PR.<t> = R[]
            sage: PS = S['t']
            sage: phi = PR.hom(f,PS)
            sage: type(phi)
            <type 'sage.rings.morphism.RingHomomorphism_from_base'>
            sage: psi = copy(phi); psi    # indirect doctest
            Ring morphism:
              From: Univariate Polynomial Ring in t over Multivariate Polynomial Ring in x, y over Rational Field
              To:   Univariate Polynomial Ring in t over Univariate Polynomial Ring in z over Rational Field
              Defn: Induced from base ring by
                    Ring morphism:
                      From: Multivariate Polynomial Ring in x, y over Rational Field
                      To:   Univariate Polynomial Ring in z over Rational Field
                      Defn: x |--> 2*z
                            y |--> 3*z
            sage: psi(x*t)
            2*z*t
        """
        slots = RingHomomorphism._extra_slots(self)
        slots['__underlying'] = self.__underlying
        return slots

    cpdef _richcmp_(self, other, int op):
        r"""
        EXAMPLES:

        A multivariate polynomial ring over a single variate quotient over
        `\QQ`::

            sage: R.<x> = QQ[]
            sage: Q.<a> = R.quotient(x^2 + x + 1)
            sage: f1 = R.hom([a])
            sage: f2 = R.hom([a + a^2 + a + 1])
            sage: PR.<s,t> = R[]
            sage: PQ = Q['s','t']
            sage: f1P = PR.hom(f1,PQ)
            sage: f2P = PR.hom(f2,PQ)
            sage: f1P == f2P
            True

        TESTS::

            sage: f1P == loads(dumps(f1P))
            True

            sage: R.<x,y> = QQ[]; f = R.hom([x,x+y]); g = R.hom([y,x])
            sage: S.<z> = R[]
            sage: fS = S.hom(f,S); gS = S.hom(g,S)
            sage: fS != gS   # indirect doctest
            True

        EXAMPLES:

        A matrix ring over a multivariate quotient over a finite field::

            sage: R.<x,y> = GF(7)[]
            sage: Q.<a,b> = R.quotient([x^2 + x + 1, y^2 + y + 1])
            sage: f1 = R.hom([a, b])
            sage: f2 = R.hom([a + a^2 + a + 1, b + b^2 + b + 1])
            sage: MR = MatrixSpace(R,2)
            sage: MQ = MatrixSpace(Q,2)
            sage: f1M = MR.hom(f1,MQ)
            sage: f2M = MR.hom(f2,MQ)
            sage: f1M == f2M
            True

        TESTS::

            sage: f1M == loads(dumps(f1M))
            True
        """
        if not isinstance(other, RingHomomorphism_from_base):
            # Generic comparison
            return RingMap._richcmp_(self, other, op)
        self_underlying = self.__underlying
        other_underlying = (<RingHomomorphism_from_base>other).__underlying
        return richcmp(self_underlying, other_underlying, op)

    def _repr_defn(self):
        """
        Used in constructing string representation of ``self``.

        EXAMPLES:

        We use a matrix ring over univariate polynomial ring over the fraction field
        over a multivariate polynomial ring::

            sage: R1.<x,y> = ZZ[]
            sage: f = R1.hom([x+y,x-y])
            sage: R2 = MatrixSpace(FractionField(R1)['t'],2)
            sage: g = R2.hom(f,R2)
            sage: g         #indirect doctest
            Ring endomorphism of Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in t over Fraction Field of Multivariate Polynomial Ring in x, y over Integer Ring
              Defn: Induced from base ring by
                    Ring endomorphism of Univariate Polynomial Ring in t over Fraction Field of Multivariate Polynomial Ring in x, y over Integer Ring
                      Defn: Induced from base ring by
                            Ring endomorphism of Fraction Field of Multivariate Polynomial Ring in x, y over Integer Ring
                              Defn: x |--> x + y
                                    y |--> x - y
        """
        U = repr(self.__underlying).split('\n')
        return 'Induced from base ring by\n'+'\n'.join(U)

    cpdef Element _call_(self, x):
        """
        Evaluate this homomorphism at ``x``.

        EXAMPLES::

            sage: R1.<x,y> = ZZ[]
            sage: f = R1.hom([x+y,x-y])
            sage: f(2*x + y + 2) # indirect doctest
            3*x + y + 2
        """
        P = self.codomain()
        try:
            return P(dict([(a, self.__underlying(b)) for a,b in x.dict().items()]))
        except Exception:
            pass
        try:
            return P([self.__underlying(b) for b in x])
        except Exception:
            pass
        try:
            return P(self.__underlying(x.numerator()))/P(self.__underlying(x.denominator()))
        except Exception:
            raise TypeError("invalid argument %s" % repr(x))


cdef class RingHomomorphism_cover(RingHomomorphism):
    r"""
    A homomorphism induced by quotienting a ring out by an ideal.

    EXAMPLES::

        sage: R.<x,y> = PolynomialRing(QQ, 2)
        sage: S.<a,b> = R.quo(x^2 + y^2)
        sage: phi = S.cover(); phi
        Ring morphism:
          From: Multivariate Polynomial Ring in x, y over Rational Field
          To:   Quotient of Multivariate Polynomial Ring in x, y over Rational Field by the ideal (x^2 + y^2)
          Defn: Natural quotient map
        sage: phi(x+y)
        a + b
    """
    def __init__(self, parent):
        """
        Create a covering ring homomorphism, induced by quotienting out by an
        ideal.

        EXAMPLES::

            sage: f = Zmod(6).cover(); f    # implicit test
            Ring morphism:
              From: Integer Ring
              To:   Ring of integers modulo 6
              Defn: Natural quotient map
            sage: type(f)
            <type 'sage.rings.morphism.RingHomomorphism_cover'>
        """
        RingHomomorphism.__init__(self, parent)

    cpdef Element _call_(self, x):
        """
        Evaluate this covering homomorphism at ``x``, which just involves
        coercing ``x`` into the domain, then codomain.

        EXAMPLES::

            sage: f = Zmod(6).cover()
            sage: type(f)
            <type 'sage.rings.morphism.RingHomomorphism_cover'>
            sage: f(-5)                 # indirect doctest
            1

        TESTS:

        We verify that calling directly raises the expected error
        (just coercing into the codomain), but calling with __call__
        (the second call below) gives a TypeError since 1/2 can't be
        coerced into the domain. ::

            sage: f._call_(1/2)
            Traceback (most recent call last):
            ...
            ZeroDivisionError: inverse of Mod(2, 6) does not exist
            sage: f(1/2)
            Traceback (most recent call last):
            ...
            TypeError: 1/2 fails to convert into the map's domain Integer Ring, but a `pushforward` method is not properly implemented
        """
        return self.codomain()(x)

    def _repr_defn(self):
        """
        Used internally for printing covering morphisms.

        EXAMPLES::

            sage: f = Zmod(6).cover()
            sage: f._repr_defn()
            'Natural quotient map'
            sage: type(f)
            <type 'sage.rings.morphism.RingHomomorphism_cover'>
        """
        return "Natural quotient map"

    def kernel(self):
        """
        Return the kernel of this covering morphism, which is the ideal that
        was quotiented out by.

        EXAMPLES::

            sage: f = Zmod(6).cover()
            sage: f.kernel()
            Principal ideal (6) of Integer Ring
        """
        return self.codomain().defining_ideal()

    cpdef _richcmp_(self, other, int op):
        """
        Compare ``self`` to ``other``.

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(QQ, 2)
            sage: S.<a,b> = R.quo(x^2 + y^2)
            sage: phi = S.cover()
            sage: phi == loads(dumps(phi))
            True
            sage: phi == R.quo(x^2 + y^3).cover()
            False
        """
        if not isinstance(other, RingHomomorphism_cover):
            # Generic comparison
            return RingMap._richcmp_(self, other, op)
        # Two cover maps with the same parent must be equal
        return rich_to_bool(op, 0)

    def __hash__(self):
        """
        Return the hash of this morphism.

        TESTS::

            sage: R.<x,y> = PolynomialRing(QQ, 2)
            sage: S.<a,b> = R.quo(x^2 + y^2)
            sage: phi = S.cover()
            sage: type(phi)
            <type 'sage.rings.morphism.RingHomomorphism_cover'>
            sage: hash(phi) == hash(phi)
            True
            sage: {phi: 1}[phi]
            1
        """
        return hash((self.domain(), self.codomain()))


cdef class RingHomomorphism_from_quotient(RingHomomorphism):
    r"""
    A ring homomorphism with domain a generic quotient ring.

    INPUT:

    -  ``parent`` -- a ring homset ``Hom(R,S)``

    -  ``phi`` -- a ring homomorphism ``C --> S``, where ``C`` is the
       domain of ``R.cover()``

    OUTPUT: a ring homomorphism

    The domain `R` is a quotient object `C \to R`, and
    ``R.cover()`` is the ring homomorphism
    `\varphi: C \to R`. The condition on the elements
    ``im_gens`` of `S` is that they define a
    homomorphism `C \to S` such that each generator of the
    kernel of `\varphi` maps to `0`.

    EXAMPLES::

        sage: R.<x, y, z> = PolynomialRing(QQ, 3)
        sage: S.<a, b, c> = R.quo(x^3 + y^3 + z^3)
        sage: phi = S.hom([b, c, a]); phi
        Ring endomorphism of Quotient of Multivariate Polynomial Ring in x, y, z over Rational Field by the ideal (x^3 + y^3 + z^3)
          Defn: a |--> b
                b |--> c
                c |--> a
        sage: phi(a+b+c)
        a + b + c
        sage: loads(dumps(phi)) == phi
        True

    Validity of the homomorphism is determined, when possible, and a
    ``TypeError`` is raised if there is no homomorphism sending the
    generators to the given images::

        sage: S.hom([b^2, c^2, a^2])
        Traceback (most recent call last):
        ...
        TypeError: images do not define a valid homomorphism
    """
    def __init__(self, parent, phi):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: R.<x,y> = QQ[]; S.<xx,yy> = R.quo([x^2,y^2]); S.hom([yy,xx])
            Ring endomorphism of Quotient of Multivariate Polynomial Ring in x, y over Rational Field by the ideal (x^2, y^2)
              Defn: xx |--> yy
                    yy |--> xx
        """
        RingHomomorphism.__init__(self, parent)
        R = parent.domain()
        pi = R.cover()  # the covering map, which should be a RingHomomorphism
        if not isinstance(pi, RingHomomorphism):
            raise TypeError("pi should be a ring homomorphism")
        if not isinstance(phi, RingHomomorphism):
            raise TypeError("phi should be a ring homomorphism")
        if pi.domain() != phi.domain():
            raise ValueError("Domain of phi must equal domain of covering (%s != %s)." % (pi.domain(), phi.domain()))
        for x in pi.kernel().gens():
            if phi(x) != 0:
                raise ValueError("relations do not all (canonically) map to 0 under map determined by images of generators")
        self._lift = pi.lift()
        self.phi = phi

    cdef _update_slots(self, dict _slots):
        """
        Helper for copying and pickling.

        EXAMPLES::

            sage: R.<x, y, z> = PolynomialRing(QQ, 3)
            sage: S.<a, b, c> = R.quo(x^3 + y^3 + z^3)
            sage: phi = S.hom([b, c, a]); phi
            Ring endomorphism of Quotient of Multivariate Polynomial Ring in x, y, z over Rational Field by the ideal (x^3 + y^3 + z^3)
              Defn: a |--> b
                    b |--> c
                    c |--> a
            sage: phi(a+b+c)
            a + b + c
            sage: psi = copy(phi)    # indirect doctest
            sage: psi == phi
            True
            sage: psi is phi
            False
            sage: psi(a) == phi(a)
            True

        """
        self.phi = _slots['phi']
        RingHomomorphism._update_slots(self, _slots)

    cdef dict _extra_slots(self):
        """
        Helper for copying and pickling.

        EXAMPLES::

            sage: R.<x, y, z> = PolynomialRing(QQ, 3)
            sage: S.<a, b, c> = R.quo(x^3 + y^3 + z^3)
            sage: phi = S.hom([b, c, a]); phi
            Ring endomorphism of Quotient of Multivariate Polynomial Ring in x, y, z over Rational Field by the ideal (x^3 + y^3 + z^3)
              Defn: a |--> b
                    b |--> c
                    c |--> a
            sage: phi(a+b+c)
            a + b + c
            sage: psi = copy(phi)    # indirect doctest
            sage: psi == phi
            True
            sage: psi is phi
            False
            sage: psi(a) == phi(a)
            True
        """
        slots = RingHomomorphism._extra_slots(self)
        slots['phi'] = self.phi
        return slots

    def _phi(self):
        """
        Underlying morphism used to define this quotient map, i.e.,
        morphism from the cover of the domain.

        EXAMPLES::

            sage: R.<x,y> = QQ[]; S.<xx,yy> = R.quo([x^2,y^2]); f = S.hom([yy,xx])
            sage: f._phi()
            Ring morphism:
              From: Multivariate Polynomial Ring in x, y over Rational Field
              To:   Quotient of Multivariate Polynomial Ring in x, y over Rational Field by the ideal (x^2, y^2)
              Defn: x |--> yy
                    y |--> xx
        """
        return self.phi

    def morphism_from_cover(self):
        """
        Underlying morphism used to define this quotient map, i.e.,
        the morphism from the cover of the domain.

        EXAMPLES::

            sage: R.<x,y> = QQ[]; S.<xx,yy> = R.quo([x^2,y^2])
            sage: S.hom([yy,xx]).morphism_from_cover()
            Ring morphism:
              From: Multivariate Polynomial Ring in x, y over Rational Field
              To:   Quotient of Multivariate Polynomial Ring in x, y over Rational Field by the ideal (x^2, y^2)
              Defn: x |--> yy
                    y |--> xx
        """
        return self.phi

    cpdef _richcmp_(self, other, int op):
        """
        Compare ``self`` to ``other``.

        EXAMPLES::

            sage: R.<x, y, z> = PolynomialRing(GF(19), 3)
            sage: S.<a, b, c> = R.quo(x^3 + y^3 + z^3)
            sage: phi = S.hom([b, c, a])
            sage: psi = S.hom([c, b, a])
            sage: f = S.hom([b, c, a + a^3 + b^3 + c^3])
            sage: phi == psi
            False
            sage: phi == f
            True
        """
        if not isinstance(other, RingHomomorphism_from_quotient):
            # Generic comparison
            return RingMap._richcmp_(self, other, op)
        # Generic comparison
        self_phi = self.phi
        other_phi = (<RingHomomorphism_from_quotient>other).phi
        return richcmp(self_phi, other_phi, op)

    def __hash__(self):
        """
        Return the hash of this morphism.

        EXAMPLES::

            sage: R.<x, y, z> = PolynomialRing(GF(19), 3)
            sage: S.<a, b, c> = R.quo(x^3 + y^3 + z^3)
            sage: phi = S.hom([b, c, a])
            sage: type(phi)
            <type 'sage.rings.morphism.RingHomomorphism_from_quotient'>
            sage: hash(phi) == hash(phi)
            True
            sage: {phi: 1}[phi]
            1
        """
        return hash(self.phi)

    def _repr_defn(self):
        """
        Used internally for printing this function.

        EXAMPLES::

            sage: R.<x,y> = QQ[]; S.<xx,yy> = R.quo([x^2,y^2]); f = S.hom([yy,xx])
            sage: print(f._repr_defn())
            xx |--> yy
            yy |--> xx
        """
        D = self.domain()
        ig = self.phi.im_gens()
        return '\n'.join(['%s |--> %s'%(D.gen(i), ig[i]) for\
                          i in range(D.ngens())])

    cpdef Element _call_(self, x):
        """
        Evaluate this function at ``x``.

        EXAMPLES::

            sage: R.<x,y> = QQ[]; S.<xx,yy> = R.quo([x^2,y^2]); f = S.hom([yy,xx])
            sage: f(3*x + (1/2)*y)   # indirect doctest
            1/2*xx + 3*yy
        """
        return self.phi(self.lift(x))


cdef class FrobeniusEndomorphism_generic(RingHomomorphism):
    """
    A class implementing Frobenius endomorphisms on rings of prime
    characteristic.
    """
    def __init__(self, domain, n=1):
        """
        INPUT:

        -  ``domain`` -- a ring

        -  ``n`` -- a nonnegative integer (default: 1)

        OUTPUT:

        The `n`-th power of the absolute (arithmetic) Frobenius
        endomorphism on ``domain``

        TESTS::

            sage: from sage.rings.morphism import FrobeniusEndomorphism_generic
            sage: K.<u> = PowerSeriesRing(GF(5))
            sage: FrobeniusEndomorphism_generic(K)
            Frobenius endomorphism x |--> x^5 of Power Series Ring in u over Finite Field of size 5
            sage: FrobeniusEndomorphism_generic(K, 2)
            Frobenius endomorphism x |--> x^(5^2) of Power Series Ring in u over Finite Field of size 5
        """
        from .ring import CommutativeRing
        from sage.categories.homset import Hom
        if not isinstance(domain, CommutativeRing):
            raise TypeError("The base ring must be a commutative ring")
        self._p = domain.characteristic()
        if not self._p.is_prime():
            raise TypeError("the characteristic of the base ring must be prime")
        try:
            n = Integer(n)
        except TypeError:
            raise TypeError("n (=%s) is not a nonnegative integer" % n)
        if n < 0:
            raise TypeError("n (=%s) is not a nonnegative integer" % n)
        self._power = n
        self._q = self._p ** self._power
        RingHomomorphism.__init__(self, Hom(domain, domain))

    def _repr_(self):
        """
        Return a string representation of this endomorphism.

        EXAMPLES::

            sage: K.<u> = PowerSeriesRing(GF(5))
            sage: Frob = K.frobenius_endomorphism(); Frob
            Frobenius endomorphism x |--> x^5 of Power Series Ring in u over Finite Field of size 5

            sage: Frob._repr_()
            'Frobenius endomorphism x |--> x^5 of Power Series Ring in u over Finite Field of size 5'
        """
        if self._power == 0:
            s = "Identity endomorphism"
        elif self._power == 1:
            s = "Frobenius endomorphism x |--> x^%s" % self._p
        else:
            s = "Frobenius endomorphism x |--> x^(%s^%s)" % (self._p, self._power)
        s += " of %s" % self.domain()
        return s

    def _repr_short(self):
        """
        Return a short string representation of this endomorphism.

        EXAMPLES::

            sage: K.<u> = PowerSeriesRing(GF(5))
            sage: Frob = K.frobenius_endomorphism()
            sage: Frob._repr_short()
            'Frob'
            sage: (Frob^2)._repr_short()
            'Frob^2'
        """
        if self._power == 0:
            s = "Identity"
        elif self._power == 1:
            s = "Frob"
        else:
            s = "Frob^%s" % self._power
        return s

    def _latex_(self):
        r"""
        Return a latex representation of this endomorphism.

        EXAMPLES::

            sage: K.<u> = PowerSeriesRing(GF(5))
            sage: Frob = K.frobenius_endomorphism(2)
            sage: Frob._latex_()
            '\\verb"Frob"^{2}'
        """
        if self._power == 0:
            s = '\\verb"id"'
        elif self._power == 1:
            s = '\\verb"Frob"'
        else:
            s = '\\verb"Frob"^{%s}' % self._power
        return s

    cpdef Element _call_ (self, x):
        """
        TESTS::

            sage: K.<u> = PowerSeriesRing(GF(5))
            sage: Frob = K.frobenius_endomorphism()
            sage: Frob(u)
            u^5
            sage: (Frob^2)(1+u)
            1 + u^25
        """
        return x ** self._q

    def power(self):
        """
        Return an integer `n` such that this endomorphism
        is the `n`-th power of the absolute (arithmetic)
        Frobenius.

        EXAMPLES::

            sage: K.<u> = PowerSeriesRing(GF(5))
            sage: Frob = K.frobenius_endomorphism()
            sage: Frob.power()
            1
            sage: (Frob^9).power()
            9
        """
        return self._power

    def __pow__(self, n, ignored):
        """
        Return the `n`-th iterate of this endomorphism.

        EXAMPLES::

            sage: K.<u> = PowerSeriesRing(GF(5))
            sage: Frob = K.frobenius_endomorphism(); Frob
            Frobenius endomorphism x |--> x^5 of Power Series Ring in u over Finite Field of size 5
            sage: Frob^2
            Frobenius endomorphism x |--> x^(5^2) of Power Series Ring in u over Finite Field of size 5
        """
        return self.__class__(self.domain(), self.power()*n)

    def _composition(self, right):
        """
        Return self o right.

        EXAMPLES::

            sage: K.<u> = PowerSeriesRing(GF(5))
            sage: f = K.frobenius_endomorphism(); f
            Frobenius endomorphism x |--> x^5 of Power Series Ring in u over Finite Field of size 5
            sage: g = K.frobenius_endomorphism(2); g
            Frobenius endomorphism x |--> x^(5^2) of Power Series Ring in u over Finite Field of size 5
            sage: f * g
            Frobenius endomorphism x |--> x^(5^3) of Power Series Ring in u over Finite Field of size 5
        """
        if isinstance(right, FrobeniusEndomorphism_generic):
            return self.__class__(self.domain(), self._power + right.power())
        else:
            return RingHomomorphism._composition(self, right)

    def __hash__(self):
        """
        Return a hash of this morphism.

        It is the hash of the triple (domain, codomain, definition)
        where ``definition`` is:

        - a tuple consisting of the images of the generators
          of the domain if domain has generators

        - the string representation of this morphism otherwise

        AUTHOR:

        - Xavier Caruso (2012-07-09)
        """
        domain = self.domain()
        codomain = self.codomain()
        return hash((domain, codomain, ('Frob', self._power)))
