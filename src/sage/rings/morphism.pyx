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
    TypeError: 2/3 fails to convert into the map's domain Integer Ring,
    but a `pushforward` method is not properly implemented

There is no homomorphism in the other direction::

    sage: H = Hom(QQ, ZZ)
    sage: H([1])
    Traceback (most recent call last):
    ...
    ValueError: relations do not all (canonically) map to 0
    under map determined by images of generators

EXAMPLES:

Reduction to finite field::

    sage: # needs sage.rings.finite_rings
    sage: H = Hom(ZZ, GF(9, 'a'))
    sage: phi = H([1])
    sage: phi(5)
    2
    sage: psi = H([4])
    sage: psi(5)
    2

Map from single variable polynomial ring::

    sage: R.<x> = ZZ[]
    sage: phi = R.hom([2], GF(5)); phi
    Ring morphism:
      From: Univariate Polynomial Ring in x over Integer Ring
      To:   Finite Field of size 5
      Defn: x |--> 2
    sage: phi(x + 12)
    4

Identity map on the real numbers::

    sage: # needs sage.rings.real_mpfr
    sage: f = RR.hom([RR(1)]); f
    Ring endomorphism of Real Field with 53 bits of precision
      Defn: 1.00000000000000 |--> 1.00000000000000
    sage: f(2.5)
    2.50000000000000
    sage: f = RR.hom([2.0])
    Traceback (most recent call last):
    ...
    ValueError: relations do not all (canonically) map to 0
    under map determined by images of generators

Homomorphism from one precision of field to another.

From smaller to bigger doesn't make sense::

    sage: R200 = RealField(200)                                                         # needs sage.rings.real_mpfr
    sage: f = RR.hom( R200 )                                                            # needs sage.rings.real_mpfr
    Traceback (most recent call last):
    ...
    TypeError: natural coercion morphism from Real Field with 53 bits of precision
    to Real Field with 200 bits of precision not defined

From bigger to small does::

    sage: f = RR.hom(RealField(15))                                                     # needs sage.rings.real_mpfr
    sage: f(2.5)                                                                        # needs sage.rings.real_mpfr
    2.500
    sage: f(RR.pi())                                                                    # needs sage.rings.real_mpfr
    3.142

Inclusion map from the reals to the complexes::

    sage: # needs sage.rings.real_mpfr
    sage: i = RR.hom([CC(1)]); i
    Ring morphism:
      From: Real Field with 53 bits of precision
      To:   Complex Field with 53 bits of precision
      Defn: 1.00000000000000 |--> 1.00000000000000
    sage: i(RR('3.1'))
    3.10000000000000

A map from a multivariate polynomial ring to itself::

    sage: R.<x,y,z> = PolynomialRing(QQ, 3)
    sage: phi = R.hom([y, z, x^2]); phi
    Ring endomorphism of Multivariate Polynomial Ring in x, y, z over Rational Field
      Defn: x |--> y
            y |--> z
            z |--> x^2
    sage: phi(x + y + z)
    x^2 + y + z

An endomorphism of a quotient of a multi-variate polynomial ring::

    sage: # needs sage.libs.singular
    sage: R.<x,y> = PolynomialRing(QQ)
    sage: S.<a,b> = quo(R, ideal(1 + y^2))
    sage: phi = S.hom([a^2, -b]); phi
    Ring endomorphism of Quotient of Multivariate Polynomial Ring in x, y
     over Rational Field by the ideal (y^2 + 1)
      Defn: a |--> a^2
            b |--> -b
    sage: phi(b)
    -b
    sage: phi(a^2 + b^2)
    a^4 - 1

The reduction map from the integers to the integers modulo 8, viewed as
a quotient ring::

    sage: R = ZZ.quo(8*ZZ)
    sage: pi = R.cover(); pi
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

    sage: # needs sage.rings.finite_rings
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

    sage: # needs sage.rings.finite_rings
    sage: pi = ZZ.hom(k); pi
    Natural morphism:
      From: Integer Ring
      To:   Finite Field of size 2
    sage: f = i * pi; f
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

    sage: # needs sage.rings.padics
    sage: phi = QQ.hom(Qp(3, print_mode='series'))
    sage: phi
    Ring morphism:
      From: Rational Field
      To:   3-adic Field with capped relative precision 20
    sage: phi.codomain()
    3-adic Field with capped relative precision 20
    sage: phi(394)
    1 + 2*3 + 3^2 + 2*3^3 + 3^4 + 3^5 + O(3^20)

An automorphism of a quotient of a univariate polynomial ring::

    sage: # needs sage.libs.pari
    sage: R.<x> = PolynomialRing(QQ)
    sage: S.<sqrt2> = R.quo(x^2 - 2)
    sage: sqrt2^2
    2
    sage: (3+sqrt2)^10
    993054*sqrt2 + 1404491
    sage: c = S.hom([-sqrt2])
    sage: c(1+sqrt2)
    -sqrt2 + 1

Note that Sage verifies that the morphism is valid::

    sage: (1 - sqrt2)^2                                                                 # needs sage.libs.pari
    -2*sqrt2 + 3
    sage: c = S.hom([1 - sqrt2])    # this is not valid                                 # needs sage.libs.pari
    Traceback (most recent call last):
    ...
    ValueError: relations do not all (canonically) map to 0
    under map determined by images of generators

Endomorphism of power series ring::

    sage: R.<t> = PowerSeriesRing(QQ, default_prec=10); R
    Power Series Ring in t over Rational Field
    sage: f = R.hom([t^2]); f
    Ring endomorphism of Power Series Ring in t over Rational Field
      Defn: t |--> t^2
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
    ValueError: relations do not all (canonically) map to 0
    under map determined by images of generators
    sage: R.hom([1])
    Traceback (most recent call last):
    ...
    ValueError: relations do not all (canonically) map to 0
    under map determined by images of generators

Complex conjugation on cyclotomic fields::

    sage: # needs sage.rings.number_field
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

    sage: # needs sage.rings.number_field
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

    sage: K = QQ  # by the way :-)
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
    sage: f(a + b + c + d)
    u
    sage: f((a+b+c+d)^2)
    u^2

TESTS::

    sage: H = Hom(ZZ, QQ)
    sage: H == loads(dumps(H))
    True

::

    sage: # needs sage.rings.number_field
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

    sage: # needs sage.rings.finite_rings
    sage: k = GF(2)
    sage: R.<x> = k[]
    sage: F4.<a> = R.quo(x^2 + x + 1)
    sage: H = End(F4)
    sage: from sage.rings.morphism import *
    sage: phi1 = H.identity(); phi1
    Identity endomorphism of Univariate Quotient Polynomial Ring in a
     over Finite Field of size 2 with modulus x^2 + x + 1
    sage: phi2 = H([a]); phi2
    Ring endomorphism of Univariate Quotient Polynomial Ring in a
     over Finite Field of size 2 with modulus x^2 + x + 1
      Defn: a |--> a
    sage: phi3 = RingHomomorphism_from_base(H, R.hom([x])); phi3                        # needs sage.libs.ntl
    Ring endomorphism of Univariate Quotient Polynomial Ring in a
     over Finite Field of size 2 with modulus x^2 + x + 1
      Defn: Induced from base ring by
            Ring endomorphism of Univariate Polynomial Ring in x
             over Finite Field of size 2 (using GF2X)
              Defn: x |--> x
    sage: phi4 = RingHomomorphism_cover(H); phi4
    Ring endomorphism of Univariate Quotient Polynomial Ring in a
     over Finite Field of size 2 with modulus x^2 + x + 1
      Defn: Natural quotient map
    sage: phi5 = F4.frobenius_endomorphism() ^ 2; phi5
    Frobenius endomorphism x |--> x^(2^2) of
     Univariate Quotient Polynomial Ring in a
     over Finite Field of size 2 with modulus x^2 + x + 1
    sage: maps = [phi1, phi2, phi3, phi4, phi5]                                         # needs sage.libs.ntl
    sage: for f in maps:                                                                # needs sage.libs.ntl
    ....:     for g in maps:
    ....:         if f != g:
    ....:             print("{} != {}".format(f, g))
"""

# ****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.rings import ideal
import sage.structure.all
from sage.structure.richcmp cimport (richcmp, rich_to_bool)
from sage.misc.cachefunc import cached_method
from sage.categories.rings import Rings
from sage.categories.facade_sets import FacadeSets


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
        <class 'sage.rings.morphism.RingMap'>
    """
    def _repr_type(self):
        """
        TESTS::

            sage: f = sage.rings.morphism.RingMap(ZZ.Hom(ZZ))
            sage: type(f)
            <class 'sage.rings.morphism.RingMap'>
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
        sage: S.<xbar,ybar> = R.quo( (x^2 + y^2, y) )                                   # needs sage.libs.singular
        sage: S.lift()                                                                  # needs sage.libs.singular
        Set-theoretic ring morphism:
          From: Quotient of Multivariate Polynomial Ring in x, y
                 over Rational Field by the ideal (x^2 + y^2, y)
          To:   Multivariate Polynomial Ring in x, y over Rational Field
          Defn: Choice of lifting map
        sage: S.lift() == 0                                                             # needs sage.libs.singular
        False

    Since :issue:`11068`, it is possible to create
    quotient rings of non-commutative rings by two-sided
    ideals. It was needed to modify :class:`RingMap_lift`
    so that rings can be accepted that are no instances
    of :class:`sage.rings.ring.Ring`, as in the following
    example::

        sage: # needs sage.modules sage.rings.finite_rings
        sage: MS = MatrixSpace(GF(5), 2, 2)
        sage: I = MS * [MS.0*MS.1, MS.2+MS.3] * MS
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
            <class 'sage.rings.integer.Integer'>
            sage: type(f)
            <class 'sage.rings.morphism.RingMap_lift'>

        An invalid example::

            sage: GF9.<one, a> = GaussianIntegers().quotient(3)                         # needs sage.rings.number_field
            sage: from sage.rings.morphism import RingMap_lift
            sage: RingMap_lift(GF9, ZZ)                                                 # needs sage.rings.number_field
            Traceback (most recent call last):
            ...
            TypeError: no canonical coercion from Number Field in I
            with defining polynomial x^2 + 1 with I = 1*I to Integer Ring
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

        Verify that :issue:`5758` has been fixed::

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
            <class 'sage.rings.morphism.RingMap_lift'>
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
            <class 'sage.rings.morphism.RingMap_lift'>
            sage: f(-1)                       # indirect doctest
            7
            sage: type(f(-1))
            <class 'sage.rings.integer.Integer'>
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

            sage: f = ZZ.hom(Zp(3)); f                                                  # needs sage.rings.padics
            Ring morphism:
              From: Integer Ring
              To:   3-adic Ring with capped relative precision 20

        TESTS::

            sage: isinstance(f, sage.rings.morphism.RingHomomorphism)                   # needs sage.rings.padics
            True
        """
        from sage.rings.homset import RingHomset_generic
        if not isinstance(parent, RingHomset_generic):
            raise TypeError("parent must be a ring homset")
        RingMap.__init__(self, parent)

    def _repr_type(self):
        """
        Used internally in printing this morphism.

        TESTS::

            sage: ZZ.hom(Zp(3))._repr_type()                                            # needs sage.rings.padics
            'Ring'
        """
        return "Ring"

    def _set_lift(self, lift):
        r"""
        Used internally to define a lifting map associated to
        this homomorphism, which goes in the other direction.  I.e.,
        if ``self`` is from `R` to `S`, then the lift must be a set-theoretic
        map from `S` to `R` such that ``self(lift(x)) == x``.

        INPUT:

        - ``lift`` -- a ring map

        OUTPUT: changes the state of ``self``

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

            sage: C = CyclotomicField(24)                                               # needs sage.rings.number_field
            sage: f = End(C)[1]                                                         # needs sage.rings.number_field
            sage: type(f*f) == type(f)                                                  # needs sage.rings.number_field
            True

        An example where the domain of ``right`` is a relative number field::

            sage: PQ.<X> = QQ[]
            sage: K.<a, b> = NumberField([X^2 - 2, X^2 - 3])                            # needs sage.rings.number_field
            sage: e, u, v, w = End(K)                                                   # needs sage.rings.number_field
            sage: u*v                                                                   # needs sage.rings.number_field
            Relative number field endomorphism of
             Number Field in a with defining polynomial X^2 - 2 over its base field
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
                      To:   Fraction Field of
                            Multivariate Polynomial Ring in a, b over Rational Field

        We check that composition works when there is a base map::

            sage: # needs sage.rings.finite_rings
            sage: R.<x> = ZZ[]
            sage: K.<a> = GF(7^2)
            sage: L.<u> = K.extension(x^3 - 3)
            sage: phi = L.hom([u^7], base_map=K.frobenius_endomorphism())
            sage: phi
            Ring endomorphism of Univariate Quotient Polynomial Ring in u
             over Finite Field in a of size 7^2 with modulus u^3 + 4
              Defn: u |--> 2*u
                    with map of base ring
            sage: psi = phi^3; psi
            Ring endomorphism of Univariate Quotient Polynomial Ring in u
             over Finite Field in a of size 7^2 with modulus u^3 + 4
              Defn: u |--> u
                    with map of base ring
            sage: psi(a) == phi(phi(phi(a)))
            True

        It also works when the image of the base map is not contained within the base ring of the codomain::

            sage: S.<x> = QQ[]
            sage: T.<y> = S[]
            sage: cc = S.hom([x + y])
            sage: f = T.hom([x - y], base_map=cc)
            sage: f*f
            Ring endomorphism of Univariate Polynomial Ring in y
             over Univariate Polynomial Ring in x over Rational Field
              Defn: y |--> 2*y
                    with map of base ring
            sage: (f*f).base_map()
            Ring morphism:
              From: Univariate Polynomial Ring in x over Rational Field
              To:   Univariate Polynomial Ring in y over
                    Univariate Polynomial Ring in x over Rational Field
              Defn: x |--> 2*x
                    with map of base ring

            sage: S.<x> = QQ[]
            sage: T.<y> = S[]
            sage: cc = S.hom([x+y])
            sage: f = T.hom([x-y], base_map=cc)
            sage: g = T.hom([x-y])
            sage: (f*g)(x)
            y + x
            sage: f(g(x))
            y + x

        AUTHORS:

        - Simon King (2010-05)
        - Francis Clarke (2011-02)
        - David Roe (2019-10)
        """
        from sage.categories.morphism import IdentityMorphism
        from sage.categories.rings import Rings
        if isinstance(right, IdentityMorphism):
            return self
        if homset.homset_category().is_subcategory(Rings()):
            if isinstance(right, RingHomomorphism_im_gens):
                rbm = right.base_map()
                kwds = {'check': False}
                if rbm is None and isinstance(self, RingHomomorphism_im_gens) and self.base_map() is not None:
                    rbm = right.codomain().coerce_map_from(right.domain().base_ring())
                if rbm is not None:
                    kwds['base_map'] = self * rbm
                try:
                    return homset([self(g) for g in right.im_gens()], **kwds)
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
        Return the pushforward of the ideal `I` under this ring
        homomorphism.

        EXAMPLES::

            sage: R.<x,y> = QQ[]; S.<xx,yy> = R.quo([x^2, y^2]); f = S.cover()          # needs sage.libs.singular
            sage: f.pushforward(R.ideal([x, 3*x + x*y + y^2]))                          # needs sage.libs.singular
            Ideal (xx, xx*yy + 3*xx) of Quotient of Multivariate Polynomial Ring
             in x, y over Rational Field by the ideal (x^2, y^2)
        """
        if not isinstance(I, ideal.Ideal_generic):
            raise TypeError("I must be an ideal")
        R = self.codomain()
        return R.ideal([self(y) for y in I.gens()])

    def inverse_image(self, I):
        """
        Return the inverse image of an ideal or an element in the codomain
        of this ring homomorphism.

        INPUT:

        - ``I`` -- an ideal or element in the codomain

        OUTPUT:

        For an ideal `I` in the codomain, this returns the largest ideal in the
        domain whose image is contained in `I`.

        Given an element `b` in the codomain, this returns an arbitrary element
        `a` in the domain such that ``self(a) = b`` if one such exists.
        The element `a` is unique if this ring homomorphism is injective.

        EXAMPLES::

            sage: R.<x,y,z> = QQ[]
            sage: S.<u,v> = QQ[]
            sage: f = R.hom([u^2, u*v, v^2], S)
            sage: I = S.ideal([u^6, u^5*v, u^4*v^2, u^3*v^3])
            sage: J = f.inverse_image(I); J                                             # needs sage.libs.singular
            Ideal (y^2 - x*z, x*y*z, x^2*z, x^2*y, x^3)
            of Multivariate Polynomial Ring in x, y, z over Rational Field
            sage: f(J) == I                                                             # needs sage.libs.singular
            True

        Under the above homomorphism, there exists an inverse image for
        every element that only involves monomials of even degree::

            sage: [f.inverse_image(p) for p in [u^2, u^4, u*v + u^3*v^3]]               # needs sage.libs.singular
            [x, x^2, x*y*z + y]
            sage: f.inverse_image(u*v^2)                                                # needs sage.libs.singular
            Traceback (most recent call last):
            ...
            ValueError: element u*v^2 does not have preimage

        The image of the inverse image ideal can be strictly smaller than the
        original ideal::

            sage: # needs sage.libs.singular sage.rings.number_field
            sage: S.<u,v> = QQ['u,v'].quotient('v^2 - 2')
            sage: f = QuadraticField(2).hom([v], S)
            sage: I = S.ideal(u + v)
            sage: J = f.inverse_image(I)
            sage: J.is_zero()
            True
            sage: f(J) < I
            True

        Fractional ideals are not yet fully supported::

            sage: # needs sage.rings.number_field
            sage: K.<a> = NumberField(QQ['x']('x^2+2'))
            sage: f = K.hom([-a], K)
            sage: I = K.ideal([a + 1])
            sage: f.inverse_image(I)                                                    # needs sage.libs.singular
            Traceback (most recent call last):
            ...
            NotImplementedError: inverse image not implemented...
            sage: f.inverse_image(K.ideal(0)).is_zero()                                 # needs sage.libs.singular
            True
            sage: f.inverse()(I)                                                        # needs sage.rings.padics
            Fractional ideal (-a + 1)

        ALGORITHM:

        By default, this computes a Gröbner basis of an ideal related to the
        graph of the ring homomorphism.

        REFERENCES:

        - Proposition 2.5.12 [DS2009]_

        TESTS::

            sage: ZZ.hom(Zp(2)).inverse_image(ZZ.ideal(2))                              # needs sage.rings.padics
            Traceback (most recent call last):
            ...
            ValueError: not an ideal or element in codomain 2-adic Ring
            with capped relative precision 20

        ::

            sage: ZZ.hom(Zp(2)).inverse_image(Zp(2).ideal(2))                           # needs sage.rings.padics
            Traceback (most recent call last):
            ...
            NotImplementedError: base rings must be equal
        """
        from sage.categories.ring_ideals import RingIdeals
        B = self.codomain()
        if I in RingIdeals(B):
            return self._inverse_image_ideal(I)
        elif I in B:
            return self._inverse_image_element(I)
        else:
            raise ValueError("not an ideal or element in codomain %s" % B)

    def _inverse_image_ideal(self, I):
        """
        Return the inverse image of an ideal under this ring homomorphism.

        EXAMPLES::

            sage: # needs sage.rings.number_field sage.symbolic
            sage: R.<x,y> = QQbar[]
            sage: f = R.hom([x, QQbar(i) * x + y^2], R)
            sage: I = R.ideal(y^3)
            sage: J = f._inverse_image_ideal(I); J
            Ideal (x^2 + 2*I*x*y - y^2)
             of Multivariate Polynomial Ring in x, y over Algebraic Field
            sage: f(J) <= I
            True

        TESTS:

        Check that :issue:`31367` is fixed::

            sage: A.<t> = QQ[]
            sage: B.<x,y> = QQ['x,y'].quotient('y')                                     # needs sage.libs.singular
            sage: f = A.hom([x], B)                                                     # needs sage.libs.singular
            sage: f.kernel()                                                            # needs sage.libs.singular
            Principal ideal (0) of Univariate Polynomial Ring in t over Rational Field

        ::

            sage: A.<t,u> = QQ[]
            sage: B.<x,y,z> = QQ['x,y,z'].quotient('z')                                 # needs sage.libs.singular
            sage: f = A.hom([x, y], B)                                                  # needs sage.libs.singular
            sage: f.kernel()                                                            # needs sage.libs.singular
            Ideal (0) of Multivariate Polynomial Ring in t, u over Rational Field
        """
        from sage.rings.polynomial.polynomial_quotient_ring import PolynomialQuotientRing_generic
        from sage.rings.quotient_ring import QuotientRing_nc
        from sage.rings.polynomial.multi_polynomial_ring import MPolynomialRing_base
        from sage.rings.polynomial.polynomial_ring import PolynomialRing_generic
        B = self.codomain()
        graph, from_B, to_A = self._graph_ideal()
        Q = graph.ring()
        gens_B = Q.gens()[:B.ngens()]
        if I.is_zero():
            # avoid adding the 0-ideal to the graph ideal in order to benefit
            # from a cached Gröbner basis
            graph_I = graph
        elif isinstance(B, (MPolynomialRing_base,
                            PolynomialRing_generic,
                            QuotientRing_nc,
                            PolynomialQuotientRing_generic)):
            graph_I = graph + from_B(I)
        else:
            # nonzero fractional ideals of number fields not yet supported
            raise NotImplementedError("inverse image not implemented "
                                      "for ideals in %s" % B)
        if isinstance(Q, QuotientRing_nc):
            # elimination_ideal does not work with quotient rings, so
            # switch to the cover ring
            gens_B_lifted = Q.cover_ring().gens()[:B.ngens()]
            graph_I_lifted = Q.cover()._inverse_image_ideal(graph_I)
            preimage = graph_I_lifted.elimination_ideal(gens_B_lifted)
            _, ambient_to_A = to_A
            return ambient_to_A(preimage)
        else:
            preimage = graph_I.elimination_ideal(gens_B)
            return to_A(preimage)

    def _inverse_image_element(self, b):
        """
        Return an element `a` such that ``self(a) = b`` if one such exists.

        TESTS:

        A degenerate case::

            sage: R.<x,y> = QQ['x,y'].quotient(1)                                       # needs sage.libs.singular
            sage: f = R.hom([y, x], R)                                                  # needs sage.libs.singular
            sage: f.inverse_image(x), f.inverse_image(y)  # indirect doctest            # needs sage.libs.singular
            (0, 0)

        Check cases involving quotient rings in which a generator is constant
        (:issue:`31178`)::

            sage: # needs sage.libs.singular
            sage: R.<x,y> = QQ[]
            sage: B.<c,d> = R.quotient(R.ideal(x))
            sage: g = R.hom([d^2, d^3], B)
            sage: g.inverse_image(d)
            Traceback (most recent call last):
            ...
            ValueError: element d does not have preimage
            sage: g.inverse_image(d^2)
            x
            sage: g.inverse_image(d^3)
            y
            sage: A.<a,b> = R.quotient(R.ideal(y^2 - x^3))
            sage: h = A.hom([d^2, d^3], B)
            sage: h.inverse_image(d^2)
            a

        Check that quotient rings are handled correctly (:issue:`33217`)::

            sage: # needs sage.libs.singular
            sage: A.<x,y,z> = QQ['X,Y,Z'].quotient('X^2+Y^2+Z^2-1')
            sage: B.<t,u,v,w> = QQ['T,U,V,W'].quotient(['T^2+U^2-1', 'V^2+W^2-1'])
            sage: psi = A.hom([v*u, w*u, t], B)
            sage: psi.inverse_image(t^2) == z^2
            True
        """
        graph, from_B, to_A = self._graph_ideal()
        gens_A = graph.ring().gens()[-self.domain().ngens():]
        a = graph.reduce(from_B(b))
        if not all(x in gens_A for x in a.lm().variables()):
            raise ValueError(f"element {b} does not have preimage")
        return to_A(a)

    @cached_method
    def kernel(self):
        """
        Return the kernel ideal of this ring homomorphism.

        EXAMPLES::

            sage: A.<x,y> = QQ[]
            sage: B.<t> = QQ[]
            sage: f = A.hom([t^4, t^3 - t^2], B)
            sage: f.kernel()                                                            # needs sage.libs.singular
            Ideal (y^4 - x^3 + 4*x^2*y - 2*x*y^2 + x^2)
            of Multivariate Polynomial Ring in x, y over Rational Field

        We express a Veronese subring of a polynomial ring as a quotient ring::

            sage: A.<a,b,c,d> = QQ[]
            sage: B.<u,v> = QQ[]
            sage: f = A.hom([u^3, u^2*v, u*v^2, v^3], B)
            sage: f.kernel() == A.ideal(matrix.hankel([a, b, c], [d]).minors(2))        # needs sage.libs.singular
            True
            sage: Q = A.quotient(f.kernel())                                            # needs sage.libs.singular
            sage: Q.hom(f.im_gens(), B).is_injective()                                  # needs sage.libs.singular
            True

        The Steiner-Roman surface::

            sage: R.<x,y,z> = QQ[]
            sage: S = R.quotient(x^2 + y^2 + z^2 - 1)
            sage: f = R.hom([x*y, x*z, y*z], S)                                         # needs sage.libs.singular
            sage: f.kernel()                                                            # needs sage.libs.singular
            Ideal (x^2*y^2 + x^2*z^2 + y^2*z^2 - x*y*z)
             of Multivariate Polynomial Ring in x, y, z over Rational Field

        TESTS:

        The results are cached::

            sage: f.kernel() is f.kernel()                                              # needs sage.libs.singular
            True

        A degenerate case::

            sage: R.<x,y> = QQ[]
            sage: f = R.hom([0, 0], R.quotient(1))                                      # needs sage.libs.singular
            sage: f.kernel().is_one()                                                   # needs sage.libs.singular
            True

        ::

            sage: K.<sqrt2> = QuadraticField(2)                                         # needs sage.rings.number_field
            sage: K.hom([-sqrt2], K).kernel().is_zero()                                 # needs sage.libs.singular sage.rings.number_field
            True

        ::

            sage: # needs sage.rings.number_field
            sage: A.<a> = QuadraticField(2)
            sage: B.<b> = A.extension(A['b']('b^2-3'))
            sage: C.<c> = B.absolute_field()
            sage: A.hom([B(a)], C).kernel().is_zero()                                   # needs sage.libs.singular
            True
            sage: A.hom([a], B).kernel()
            Traceback (most recent call last):
            ...
            NotImplementedError: base rings must be equal
        """
        return self._inverse_image_ideal(self.codomain().zero_ideal())

    def lift(self, x=None):
        """
        Return a lifting map associated to this homomorphism, if
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

    @cached_method
    def _graph_ideal(self):
        """
        Return the ideal corresponding to the graph of this ring homomorphism.

        OUTPUT:

        - the graph as an ideal in the tensor product of codomain and domain
        - a map from the codomain to the ring of the graph ideal
        - a map from the ring of the graph ideal to the domain

        The second map is only meaningful for those elements that involve only
        variables of the domain of ``self``.

        EXAMPLES::

            sage: R.<x,y> = QQ[]
            sage: QQ['t'].hom([x*y^2], R)._graph_ideal()
            (Ideal (x*y^2 - t) of Multivariate Polynomial Ring in x, y, t over
               Rational Field,
             Ring morphism:
               From: Multivariate Polynomial Ring in x, y over Rational Field
               To:   Multivariate Polynomial Ring in x, y, t over Rational Field
               Defn: x |--> x
                     y |--> y,
             Ring morphism:
               From: Multivariate Polynomial Ring in x, y, t over Rational Field
               To:   Univariate Polynomial Ring in t over Rational Field
               Defn: x |--> 0
                     y |--> 0
                     t |--> t)

        TESTS:

        Ideals in quotient rings over ``QQbar`` do not support reduction yet,
        so the graph is constructed in the ambient ring instead::

            sage: # needs sage.libs.singular sage.rings.number_field
            sage: A.<z,w> = QQbar['z,w'].quotient('z*w - 1')
            sage: B.<x,y> = QQbar['x,y'].quotient('2*x^2 + y^2 - 1')
            sage: f = A.hom([QQbar(2).sqrt()*x + QQbar(I)*y,
            ....:            QQbar(2).sqrt()*x - QQbar(I)*y], B)
            sage: f._graph_ideal()[0]
            Ideal (z*w - 1, 2*x^2 + y^2 - 1,
            1.414213562373095?*x + I*y - z,
            1.414213562373095?*x + (-I)*y - w)
            of Multivariate Polynomial Ring in x, y, z, w over Algebraic Field
            sage: f.inverse()(f(z)), f.inverse()(f(w))
            (z, w)

        Non-trivial base maps are not supported::

            sage: # needs sage.rings.number_field
            sage: K.<a> = QuadraticField(2)
            sage: R.<x,y> = K[]
            sage: f = R.hom([x, a*x + y], R, base_map=K.hom([-a], K))
            sage: f._graph_ideal()
            Traceback (most recent call last):
            ...
            NotImplementedError: base map must be trivial

        Non-commutative rings are not supported (:issue:`32824`)::

            sage: A = GradedCommutativeAlgebra(QQ, 'x,y,z')                             # needs sage.combinat sage.modules
            sage: A.hom(A.gens(), A).kernel()                                           # needs sage.combinat sage.modules
            Traceback (most recent call last):
            ...
            NotImplementedError: rings are not commutative
        """
        from sage.rings.quotient_ring import QuotientRing_nc
        from sage.rings.ideal import Ideal_generic
        A = self.domain()
        B = self.codomain()
        Comm = Rings().Commutative()
        if not (A in Comm and B in Comm):
            raise NotImplementedError("rings are not commutative")
        if A.base_ring() != B.base_ring():
            raise NotImplementedError("base rings must be equal")
        try:
            base_map = self.base_map()
        except AttributeError:
            pass
        else:
            if base_map is not None:
                raise NotImplementedError("base map must be trivial")
        Q = _tensor_product_ring(B, A)
        A_to_Q = A.hom(Q.gens()[B.ngens():], Q, check=False)
        B_to_Q = B.hom(Q.gens()[:B.ngens()], Q, check=False)
        graph = Q.ideal([B_to_Q(self(x)) - A_to_Q(x) for x in A.gens()])
        R = Q.cover_ring() if isinstance(Q, QuotientRing_nc) else Q
        R_to_A = R.hom(tuple([0] * B.ngens()) + A.gens(), A, check=False)
        Q_to_A = R_to_A if R is Q else R_to_A * Q.lifting_map()

        # Since we compute normal forms modulo the graph ideal, check that
        # the default `reduce` method has been overwritten
        if graph.reduce.__func__ is Ideal_generic.reduce:
            if Q is not R:
                # Although the graph naturally lives in the quotient Q, we try
                # to lift it to the ambient R as a workaround, since in some
                # cases (e.g. over QQbar) reduction is supported in R
                graph_R = Q.cover()._inverse_image_ideal(graph)
                if graph_R.reduce.__func__ is not Ideal_generic.reduce:
                    return graph_R, (Q.lifting_map() * B_to_Q), R_to_A
            raise NotImplementedError('"reduce" not implemented for %s' % Q)
        return graph, B_to_Q, Q_to_A

    @cached_method
    def inverse(self):
        """
        Return the inverse of this ring homomorphism if it exists.

        Raises a :exc:`ZeroDivisionError` if the inverse does not exist.

        ALGORITHM:

        By default, this computes a Gröbner basis of the ideal corresponding to
        the graph of the ring homomorphism.

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: f = R.hom([2*t - 1], R)
            sage: f.inverse()                                                           # needs sage.libs.singular
            Ring endomorphism of Univariate Polynomial Ring in t over Rational Field
              Defn: t |--> 1/2*t + 1/2

        The following non-linear homomorphism is not invertible, but it induces
        an isomorphism on a quotient ring::

            sage: # needs sage.libs.singular
            sage: R.<x,y,z> = QQ[]
            sage: f = R.hom([y*z, x*z, x*y], R)
            sage: f.inverse()
            Traceback (most recent call last):
            ...
            ZeroDivisionError: ring homomorphism not surjective
            sage: f.is_injective()
            True
            sage: Q.<x,y,z> = R.quotient(x*y*z - 1)
            sage: g = Q.hom([y*z, x*z, x*y], Q)
            sage: g.inverse()
            Ring endomorphism of Quotient of Multivariate Polynomial Ring
            in x, y, z over Rational Field by the ideal (x*y*z - 1)
              Defn: x |--> y*z
                    y |--> x*z
                    z |--> x*y

        Homomorphisms over the integers are supported::

            sage: S.<x,y> = ZZ[]
            sage: f = S.hom([x + 2*y, x + 3*y], S)
            sage: f.inverse()                                                           # needs sage.libs.singular
            Ring endomorphism of Multivariate Polynomial Ring in x, y over Integer Ring
              Defn: x |--> 3*x - 2*y
                    y |--> -x + y
            sage: (f.inverse() * f).is_identity()                                       # needs sage.libs.singular
            True

        The following homomorphism is invertible over the rationals, but not
        over the integers::

            sage: g = S.hom([x + y, x - y - 2], S)
            sage: g.inverse()                                                           # needs sage.libs.singular
            Traceback (most recent call last):
            ...
            ZeroDivisionError: ring homomorphism not surjective
            sage: R.<x,y> = QQ[x,y]
            sage: h = R.hom([x + y, x - y - 2], R)
            sage: (h.inverse() * h).is_identity()                                       # needs sage.libs.singular
            True

        This example by M. Nagata is a wild automorphism::

            sage: R.<x,y,z> = QQ[]
            sage: sigma = R.hom([x - 2*y*(z*x+y^2) - z*(z*x+y^2)^2,
            ....:                y + z*(z*x+y^2), z], R)
            sage: tau = sigma.inverse(); tau                                            # needs sage.libs.singular
            Ring endomorphism of Multivariate Polynomial Ring in x, y, z over
            Rational Field
              Defn: x |--> -y^4*z - 2*x*y^2*z^2 - x^2*z^3 + 2*y^3 + 2*x*y*z + x
                    y |--> -y^2*z - x*z^2 + y
                    z |--> z
            sage: (tau * sigma).is_identity()                                           # needs sage.libs.singular
            True

        We compute the triangular automorphism that converts moments to
        cumulants, as well as its inverse, using the moment generating
        function. The choice of a term ordering can have a great impact on the
        computation time of a Gröbner basis, so here we choose a weighted
        ordering such that the images of the generators are homogeneous
        polynomials.  ::

            sage: d = 12
            sage: T = TermOrder('wdegrevlex', [1..d])
            sage: R = PolynomialRing(QQ, ['x%s' % j for j in (1..d)], order=T)
            sage: S.<t> = PowerSeriesRing(R)
            sage: egf = S([0] + list(R.gens())).ogf_to_egf().exp(prec=d+1)
            sage: phi = R.hom(egf.egf_to_ogf().list()[1:], R)
            sage: phi.im_gens()[:5]
            [x1,
             x1^2 + x2,
             x1^3 + 3*x1*x2 + x3,
             x1^4 + 6*x1^2*x2 + 3*x2^2 + 4*x1*x3 + x4,
             x1^5 + 10*x1^3*x2 + 15*x1*x2^2 + 10*x1^2*x3 + 10*x2*x3 + 5*x1*x4 + x5]
            sage: all(p.is_homogeneous() for p in phi.im_gens())                        # needs sage.libs.singular
            True
            sage: phi.inverse().im_gens()[:5]                                           # needs sage.libs.singular
            [x1,
             -x1^2 + x2,
             2*x1^3 - 3*x1*x2 + x3,
             -6*x1^4 + 12*x1^2*x2 - 3*x2^2 - 4*x1*x3 + x4,
             24*x1^5 - 60*x1^3*x2 + 30*x1*x2^2 + 20*x1^2*x3 - 10*x2*x3 - 5*x1*x4 + x5]
            sage: (phi.inverse() * phi).is_identity()                                   # needs sage.libs.singular
            True

        Automorphisms of number fields as well as Galois fields are supported::

            sage: K.<zeta7> = CyclotomicField(7)                                        # needs sage.rings.number_field
            sage: c = K.hom([1/zeta7])                                                  # needs sage.rings.number_field
            sage: (c.inverse() * c).is_identity()                                       # needs sage.libs.singular sage.rings.number_field
            True

            sage: F.<t> = GF(7^3)                                                       # needs sage.rings.finite_rings
            sage: f = F.hom(t^7, F)                                                     # needs sage.rings.finite_rings
            sage: (f.inverse() * f).is_identity()                                       # needs sage.libs.singular sage.rings.finite_rings
            True

        An isomorphism between the algebraic torus and the circle over a number
        field::

            sage: # needs sage.libs.singular sage.rings.number_field
            sage: K.<i> = QuadraticField(-1)
            sage: A.<z,w> = K['z,w'].quotient('z*w - 1')
            sage: B.<x,y> = K['x,y'].quotient('x^2 + y^2 - 1')
            sage: f = A.hom([x + i*y, x - i*y], B)
            sage: g = f.inverse()
            sage: g.morphism_from_cover().im_gens()
            [1/2*z + 1/2*w, (-1/2*i)*z + (1/2*i)*w]
            sage: all(g(f(z)) == z for z in A.gens())
            True

        TESTS:

        Morphisms involving quotient rings::

            sage: # needs sage.libs.singular
            sage: R.<x,y> = QQ[]
            sage: S.<s,u,t> = QQ['s,u,t'].quotient('u-t^2')
            sage: f = R.hom([s, -t], S)
            sage: (f.inverse() * f).is_identity()
            True
            sage: Q.<v,w> = R.quotient(x - y^2)
            sage: g = Q.hom([v, -w], Q)
            sage: g.inverse()(g(v)) == v and g.inverse()(g(w)) == w
            True
            sage: S.<z> = QQ[]
            sage: h = Q.hom([z^2, -z], S)
            sage: h.inverse()(h(v)) == v and h.inverse()(h(w)) == w
            True

        Morphisms between number fields and quotient rings::

            sage: # needs sage.rings.number_field
            sage: K.<sqrt2> = QuadraticField(2)
            sage: f = K.hom([-sqrt2], K.polynomial_quotient_ring())
            sage: (f.inverse() * f).is_identity()
            True
            sage: g = K.polynomial_quotient_ring().hom([-sqrt2], K)
            sage: (g.inverse() * g).is_identity()
            True

        Morphisms involving Galois fields::

            sage: # needs sage.rings.finite_rings
            sage: A.<t> = GF(7^3)
            sage: R = A.polynomial_ring().quotient(A.polynomial())
            sage: g = A.hom(R.gens(), R)
            sage: (g.inverse() * g).is_identity()                                       # needs sage.libs.singular
            True
            sage: B.<T>, f = A.extension(3, map=True)
            sage: f.inverse()
            Traceback (most recent call last):
            ...
            ZeroDivisionError: ring homomorphism not surjective
            sage: B.<T>, f = A.extension(1, map=True)
            sage: f.inverse()
            Ring morphism:
              From: Finite Field in T of size 7^3
              To:   Finite Field in t of size 7^3
              Defn: T |--> t

        Non-injective homomorphisms::

            sage: # needs sage.libs.singular
            sage: R.<x,y> = QQ[]
            sage: S.<a,b,c> = QQ[]
            sage: S.hom([x, y, 0], R).inverse()
            Traceback (most recent call last):
            ...
            ZeroDivisionError: ring homomorphism not injective
            sage: T.<z> = QQ[]
            sage: R.hom([2*z, 3*z], T).inverse()
            Traceback (most recent call last):
            ...
            ZeroDivisionError: ring homomorphism not injective
            sage: Q.<u,v> = R.quotient([x^5, y^4])
            sage: R.hom([u, v], Q).inverse()
            Traceback (most recent call last):
            ...
            ZeroDivisionError: ring homomorphism not injective
            sage: Q.cover().inverse()
            Traceback (most recent call last):
            ...
            ZeroDivisionError: ring homomorphism not injective

        Univariate quotient rings::

            sage: R.<t> = QQ['t'].quotient('t^5')                                       # needs sage.libs.pari
            sage: f = R.hom([2*t], R)
            sage: (f.inverse() * f).is_identity()                                       # needs sage.libs.singular
            True

        A homomorphism over ``QQbar``::

            sage: R.<x,y> = QQbar[]                                                     # needs sage.rings.number_field
            sage: f = R.hom([x + QQbar(I)*y^2, -y], R)                                  # needs sage.rings.number_field
            sage: (f.inverse() * f).is_identity()                                       # needs sage.rings.number_field
            True

        Check that results are cached::

            sage: R.<x,y> = GF(823)[]
            sage: f = R.hom([x, y + x^2], R)
            sage: f.inverse() is f.inverse()                                            # needs sage.libs.singular
            True

        Some subclasses of ring homomorphisms are not supported::

            sage: from sage.rings.morphism import FrobeniusEndomorphism_generic
            sage: K.<u> = PowerSeriesRing(GF(5))
            sage: FrobeniusEndomorphism_generic(K).inverse()
            Traceback (most recent call last):
            ...
            NotImplementedError

        ::

            sage: R.<x,y> = LaurentPolynomialRing(QQ)                                   # needs sage.modules
            sage: R.hom([y, x], R).inverse()                                            # needs sage.libs.singular sage.modules
            Traceback (most recent call last):
            ...
            NotImplementedError

        ::

            sage: K.<x> = FunctionField(QQ)
            sage: K.hom(1/x).inverse()
            Traceback (most recent call last):
            ...
            NotImplementedError: inverse not implemented...

        The implementation performs several computations that require a Gröbner
        basis of the graph ideal, so we check that the Gröbner basis is cached
        after the first such computation::

            sage: R.<x,y> = QQ[]
            sage: f = R.hom([x + 123*y^2, y], R)
            sage: f._graph_ideal()[0].groebner_basis.is_in_cache()
            False
            sage: f.is_injective()                                                      # needs sage.libs.singular
            True
            sage: f._graph_ideal()[0].groebner_basis.is_in_cache()                      # needs sage.libs.singular
            True
        """
        if not self.is_injective():
            raise ZeroDivisionError("ring homomorphism not injective")
        ys = self.codomain().gens()
        try:
            preimages = [self._inverse_image_element(y) for y in ys]
        except ValueError:
            raise ZeroDivisionError("ring homomorphism not surjective")
        return self.parent().reversed()(preimages, check=False)

    def __invert__(self):
        """
        Return the inverse of this ring homomorphism if it exists.

        This simply calls :meth:`inverse`.

        EXAMPLES::

            sage: R.<x,y> = GF(17)[]
            sage: f = R.hom([3*x, y + x^2 + x^3], R)
            sage: (f * ~f).is_identity()                                                # needs sage.libs.singular
            True
        """
        return self.inverse()

    def is_surjective(self):
        """
        Return whether this ring homomorphism is surjective.

        EXAMPLES::

            sage: R.<x,y,z> = QQ[]
            sage: R.hom([y*z, x*z, x*y], R).is_surjective()                             # needs sage.libs.singular
            False
            sage: Q.<x,y,z> = R.quotient(x*y*z - 1)                                     # needs sage.libs.singular
            sage: R.hom([y*z, x*z, x*y], Q).is_surjective()                             # needs sage.libs.singular
            True

        ALGORITHM:

        By default, this requires the computation of a Gröbner basis.
        """
        for y in self.codomain().gens():
            try:
                self._inverse_image_element(y)
            except ValueError:
                return False
        return True

    def is_invertible(self):
        """
        Return whether this ring homomorphism is bijective.

        EXAMPLES::

            sage: R.<x,y,z> = QQ[]
            sage: R.hom([y*z, x*z, x*y], R).is_invertible()                             # needs sage.libs.singular
            False
            sage: Q.<x,y,z> = R.quotient(x*y*z - 1)                                     # needs sage.libs.singular
            sage: Q.hom([y*z, x*z, x*y], Q).is_invertible()                             # needs sage.libs.singular
            True

        ALGORITHM:

        By default, this requires the computation of a Gröbner basis.
        """
        return self.is_injective() and self.is_surjective()


cdef class RingHomomorphism_im_gens(RingHomomorphism):
    """
    A ring homomorphism determined by the images of generators.
    """
    def __init__(self, parent, im_gens, check=True, base_map=None):
        """
        EXAMPLES::

            sage: R.<x,y> = QQ[]
            sage: phi = R.hom([x, x + y]); phi
            Ring endomorphism of Multivariate Polynomial Ring in x, y over Rational Field
              Defn: x |--> x
                    y |--> x + y
            sage: type(phi)
            <class 'sage.rings.morphism.RingHomomorphism_im_gens'>

        Here's another example where the domain isn't free::

            sage: S.<xx,yy> = R.quotient(x - y)                                         # needs sage.libs.singular
            sage: phi = S.hom([xx + 1, xx + 1])                                         # needs sage.libs.singular

        Note that one has to specify valid images::

            sage: phi = S.hom([xx + 1, xx - 1])                                         # needs sage.libs.singular
            Traceback (most recent call last):
            ...
            ValueError: relations do not all (canonically) map to 0
             under map determined by images of generators

        You can give a map of the base ring::

            sage: # needs sage.rings.number_field
            sage: Zx.<x> = ZZ[]
            sage: K.<i> = NumberField(x^2 + 1)
            sage: cc = K.hom([-i])
            sage: R.<t> = K[]
            sage: z = 1 + i*t + (3+4*i)*t^2
            sage: z._im_gens_(R, [t^2], base_map=cc)
            (-4*i + 3)*t^4 - i*t^2 + 1

        The base map's codomain is extended to the whole codomain::

            sage: S.<x> = QQ[]
            sage: T.<y> = S[]
            sage: cc = S.hom([x + 1])
            sage: f = T.hom([x - y], base_map=cc)
            sage: g = T.hom([x - y], base_map=cc.extend_codomain(T))
            sage: f == g
            True
            sage: f.base_map() == cc.extend_codomain(T)
            True

        There is a check option, but it may be ignored in some cases
        -- it's purpose isn't so you can lie to Sage, but to sometimes
        speed up creation of a homomorphism::

            sage: R.<x,y> = QQ[]
            sage: S.<xx,yy> = R.quotient(x - y)                                         # needs sage.libs.singular
            sage: phi = S.hom([xx + 1, xx - 1], check=False)                            # needs sage.libs.singular
            Traceback (most recent call last):
            ...
            ValueError: relations do not all (canonically) map to 0
             under map determined by images of generators
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
            if base_map is None:
                tkwds = {}
            else:
                if base_map.codomain() is not self.codomain():
                    base_map = base_map.extend_codomain(self.codomain())
                tkwds = {'base_map': base_map}
            t = parent.domain()._is_valid_homomorphism_(parent.codomain(), im_gens, **tkwds)
            if not t:
                raise ValueError("relations do not all (canonically) map to 0 under map determined by images of generators")
        if not im_gens.is_immutable():
            import copy
            im_gens = copy.copy(im_gens)
            im_gens.set_immutable()
        self._im_gens = im_gens
        self._base_map = base_map

    def im_gens(self):
        """
        Return the images of the generators of the domain.

        OUTPUT:

        - ``list`` -- a copy of the list of gens (it is safe to change this)

        EXAMPLES::

            sage: R.<x,y> = QQ[]
            sage: f = R.hom([x, x + y])
            sage: f.im_gens()
            [x, x + y]

        We verify that the returned list of images of gens is a copy,
        so changing it doesn't change ``f``::

            sage: f.im_gens()[0] = 5
            sage: f.im_gens()
            [x, x + y]
        """
        return list(self._im_gens)

    def base_map(self):
        """
        Return the map on the base ring that is part of the defining
        data for this morphism.  May return ``None`` if a coercion is used.

        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: R.<x> = ZZ[]
            sage: K.<i> = NumberField(x^2 + 1)
            sage: cc = K.hom([-i])
            sage: S.<y> = K[]
            sage: phi = S.hom([y^2], base_map=cc)
            sage: phi
            Ring endomorphism of Univariate Polynomial Ring in y
             over Number Field in i with defining polynomial x^2 + 1
              Defn: y |--> y^2
                    with map of base ring
            sage: phi(y)
            y^2
            sage: phi(i*y)
            -i*y^2
            sage: phi.base_map()
            Composite map:
              From: Number Field in i with defining polynomial x^2 + 1
              To:   Univariate Polynomial Ring in y over Number Field in i
                    with defining polynomial x^2 + 1
              Defn:   Ring endomorphism of Number Field in i with defining polynomial x^2 + 1
                      Defn: i |--> -i
                    then
                      Polynomial base injection morphism:
                      From: Number Field in i with defining polynomial x^2 + 1
                      To:   Univariate Polynomial Ring in y over Number Field in i
                            with defining polynomial x^2 + 1
        """
        return self._base_map

    cdef _update_slots(self, dict _slots):
        """
        Helper for copying and pickling.

        EXAMPLES::

            sage: R.<x,y> = QQ[]
            sage: f = R.hom([x, x + y])
            sage: g = copy(f)   # indirect doctest
            sage: g == f
            True
            sage: g is f
            False
            sage: g(y)
            x + y
        """
        self._im_gens = _slots['__im_gens'] # double underscores for legacy pickles
        self._base_map = _slots.get('_base_map')
        RingHomomorphism._update_slots(self, _slots)

    cdef dict _extra_slots(self):
        """
        Helper for copying and pickling.

        EXAMPLES::

            sage: R.<x,y> = QQ[]
            sage: f = R.hom([x, x + y])
            sage: g = copy(f)   # indirect doctest
            sage: g == f
            True
            sage: g is f
            False
            sage: g(y)
            x + y
        """
        slots = RingHomomorphism._extra_slots(self)
        slots['__im_gens'] = self._im_gens
        slots['_base_map'] = self._base_map
        return slots

    cpdef _richcmp_(self, other, int op):
        r"""
        EXAMPLES:

        A single variate quotient over `\QQ`::

            sage: # needs sage.libs.pari
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

            sage: loads(dumps(f2)) == f2                                                # needs sage.libs.pari
            True

        ::

            sage: R.<x,y> = QQ[]; f = R.hom([x, x + y]); g = R.hom([y, x])
            sage: f == g             # indirect doctest
            False

        EXAMPLES:

        A multivariate quotient over a finite field::

            sage: # needs sage.libs.singular
            sage: R.<x,y> = GF(7)[]
            sage: Q.<a,b> = R.quotient([x^2 + x + 1, y^2 + y + 1])
            sage: f1 = R.hom([a, b])
            sage: f2 = R.hom([a + a^2 + a + 1, b + b^2 + b + 1])
            sage: f1 == f2
            True
            sage: f1 == R.hom([b, a])
            False
            sage: x^3 + x + y^2
            x^3 + y^2 + x
            sage: f1(x^3 + x + y^2)
            a - b
            sage: f2(x^3 + x + y^2)
            a - b

        TESTS::

            sage: loads(dumps(f2)) == f2                                                # needs sage.libs.pari
            True

        This was fixed in :issue:`24277`::

            sage: H = End(QQ)
            sage: H(1) == H.identity()
            True
        """
        if not isinstance(other, RingHomomorphism_im_gens):
            # Generic comparison
            return RingMap._richcmp_(self, other, op)
        # Check equality using the images of the generators.
        self_im = self._im_gens
        other_im = (<RingHomomorphism_im_gens>other)._im_gens
        self_hom = self._base_map
        other_hom = (<RingHomomorphism_im_gens>other)._base_map
        return richcmp((self_im, self_hom), (other_im, other_hom), op)

    def __hash__(self):
        """
        Return the hash of this morphism.

        EXAMPLES::

            sage: R.<x> = ZZ[]
            sage: s = R.hom([x+1])
            sage: type(s)
            <class 'sage.rings.morphism.RingHomomorphism_im_gens'>
            sage: hash(s) == hash(s)
            True
            sage: {s: 1}[s]
            1
        """
        return hash((self._im_gens, self._base_map))

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
        if D in FacadeSets():
            D, = D.facade_for()
        ig = self._im_gens
        s = '\n'.join('{} |--> {}'.format(D.gen(i), ig[i])
                       for i in range(D.ngens()))
        if s and self._base_map is not None:
            s += '\nwith map of base ring'
        return s

    cpdef Element _call_(self, x):
        """
        Evaluate this homomorphism at ``x``.

        EXAMPLES::

            sage: R.<x,y,z> = ZZ[]; f = R.hom([2*x,z,y])
            sage: f(x+2*y+3*z)             # indirect doctest
            2*x + 3*y + 2*z
        """
        return x._im_gens_(self.codomain(), self.im_gens(), base_map=self.base_map())


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
          From: Univariate Polynomial Ring in t
                over Multivariate Polynomial Ring in x, y over Rational Field
          To:   Univariate Polynomial Ring in t
                over Univariate Polynomial Ring in z over Rational Field
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

        sage: # needs sage.modules
        sage: MR = MatrixSpace(R, 2, 2)
        sage: MS = MatrixSpace(S, 2, 2)
        sage: M = MR([x^2 + 1/7*x*y - y^2, -1/2*y^2 + 2*y + 1/6,
        ....:         4*x^2 - 14*x, 1/2*y^2 + 13/4*x - 2/11*y])
        sage: Mf = MR.hom(f, MS)
        sage: Mf
        Ring morphism:
          From: Full MatrixSpace of 2 by 2 dense matrices
                over Multivariate Polynomial Ring in x, y over Rational Field
          To:   Full MatrixSpace of 2 by 2 dense matrices
                over Univariate Polynomial Ring in z over Rational Field
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

        sage: # needs sage.modules
        sage: MPR = MatrixSpace(PR, 2)
        sage: MPS = MatrixSpace(PS, 2)
        sage: M = MPR([(-x + y)*t^2 + 58*t - 3*x^2 + x*y,
        ....:          (- 1/7*x*y - 1/40*x)*t^2 + (5*x^2 + y^2)*t + 2*y,
        ....:          (- 1/3*y + 1)*t^2 + 1/3*x*y + y^2 + 5/2*y + 1/4,
        ....:          (x + 6*y + 1)*t^2])
        sage: MPf = MPR.hom(f, MPS); MPf
        Ring morphism:
          From: Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial
                Ring in t over Multivariate Polynomial Ring in x, y over Rational Field
          To:   Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial
                Ring in t over Univariate Polynomial Ring in z over Rational Field
          Defn: Induced from base ring by
                Ring morphism:
                  From: Univariate Polynomial Ring in t
                        over Multivariate Polynomial Ring in x, y over Rational Field
                  To:   Univariate Polynomial Ring in t
                        over Univariate Polynomial Ring in z over Rational Field
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
            sage: f = R.hom([2*x], R)
            sage: P = MatrixSpace(R, 2).Hom(MatrixSpace(R, 2))                          # needs sage.modules
            sage: g = RingHomomorphism_from_base(P, f)                                  # needs sage.modules
            sage: g                                                                     # needs sage.modules
            Ring endomorphism of Full MatrixSpace of 2 by 2 dense matrices
             over Univariate Polynomial Ring in x over Integer Ring
              Defn: Induced from base ring by
                    Ring endomorphism of Univariate Polynomial Ring in x over Integer Ring
                      Defn: x |--> 2*x

        Note that an induced homomorphism only makes sense if domain and
        codomain are constructed in a compatible way. So, the following
        results in an error::

            sage: P = MatrixSpace(R, 2).Hom(R['t'])                                     # needs sage.modules
            sage: g = RingHomomorphism_from_base(P, f)                                  # needs sage.modules
            Traceback (most recent call last):
            ...
            ValueError: domain (Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in x over Integer Ring)
             and codomain (Univariate Polynomial Ring in t over Univariate Polynomial Ring in x over Integer Ring)
             must have the same functorial construction over their base rings
        """
        RingHomomorphism.__init__(self, parent)
        if underlying.domain() != parent.domain().base():
            raise ValueError("The given homomorphism has to have the domain %s"%parent.domain().base())
        if underlying.codomain() != parent.codomain().base():
            raise ValueError("The given homomorphism has to have the codomain %s"%parent.codomain().base())
        if parent.domain().construction()[0] != parent.codomain().construction()[0]:
            raise ValueError(f"domain ({parent.domain()}) and codomain ({parent.codomain()}) must have the same functorial construction over their base rings")
        self._underlying = underlying

    def underlying_map(self):
        """
        Return the underlying homomorphism of the base ring.

        EXAMPLES::

            sage: # needs sage.modules
            sage: R.<x,y> = QQ[]
            sage: S.<z> = QQ[]
            sage: f = R.hom([2*z, 3*z], S)
            sage: MR = MatrixSpace(R, 2)
            sage: MS = MatrixSpace(S, 2)
            sage: g = MR.hom(f, MS)
            sage: g.underlying_map() == f
            True
        """
        return self._underlying

    cdef _update_slots(self, dict _slots):
        """
        Helper for copying and pickling.

        EXAMPLES::

            sage: R.<x,y> = QQ[]
            sage: S.<z> = QQ[]
            sage: f = R.hom([2*z, 3*z],S)
            sage: PR.<t> = R[]
            sage: PS = S['t']
            sage: phi = PR.hom(f, PS)
            sage: type(phi)
            <class 'sage.rings.morphism.RingHomomorphism_from_base'>
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
        self._underlying = _slots['__underlying'] # double underscore for legacy pickles
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
            <class 'sage.rings.morphism.RingHomomorphism_from_base'>
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
        slots['__underlying'] = self._underlying
        return slots

    cpdef _richcmp_(self, other, int op):
        r"""
        EXAMPLES:

        A multivariate polynomial ring over a single variate quotient over
        `\QQ`::

            sage: # needs sage.libs.pari sage.libs.singular sage.modules
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

            sage: f1P == loads(dumps(f1P))                                              # needs sage.libs.pari sage.libs.singular sage.modules
            True

            sage: R.<x,y> = QQ[]; f = R.hom([x, x + y]); g = R.hom([y, x])
            sage: S.<z> = R[]
            sage: fS = S.hom(f, S); gS = S.hom(g, S)
            sage: fS != gS   # indirect doctest
            True

        EXAMPLES:

        A matrix ring over a multivariate quotient over a finite field::

            sage: # needs sage.libs.singular sage.modules
            sage: R.<x,y> = GF(7)[]
            sage: Q.<a,b> = R.quotient([x^2 + x + 1, y^2 + y + 1])
            sage: f1 = R.hom([a, b])
            sage: f2 = R.hom([a + a^2 + a + 1, b + b^2 + b + 1])
            sage: MR = MatrixSpace(R, 2)
            sage: MQ = MatrixSpace(Q, 2)
            sage: f1M = MR.hom(f1, MQ)
            sage: f2M = MR.hom(f2, MQ)
            sage: f1M == f2M
            True

        TESTS::

            sage: f1M == loads(dumps(f1M))                                              # needs sage.libs.singular sage.modules
            True
        """
        if not isinstance(other, RingHomomorphism_from_base):
            # Generic comparison
            return RingMap._richcmp_(self, other, op)
        self_underlying = self._underlying
        other_underlying = (<RingHomomorphism_from_base>other)._underlying
        return richcmp(self_underlying, other_underlying, op)

    def _repr_defn(self):
        """
        Used in constructing string representation of ``self``.

        EXAMPLES:

        We use a matrix ring over univariate polynomial ring over the fraction field
        over a multivariate polynomial ring::

            sage: R1.<x,y> = ZZ[]
            sage: f = R1.hom([x + y, x - y])
            sage: R2 = MatrixSpace(FractionField(R1)['t'], 2)                           # needs sage.modules
            sage: g = R2.hom(f, R2)                                                     # needs sage.modules
            sage: g         # indirect doctest                                          # needs sage.modules
            Ring endomorphism of Full MatrixSpace of 2 by 2 dense matrices
             over Univariate Polynomial Ring in t over Fraction Field of Multivariate Polynomial Ring in x, y over Integer Ring
              Defn: Induced from base ring by
                    Ring endomorphism of Univariate Polynomial Ring in t over Fraction Field of Multivariate Polynomial Ring in x, y over Integer Ring
                      Defn: Induced from base ring by
                            Ring endomorphism of Fraction Field of Multivariate Polynomial Ring in x, y over Integer Ring
                              Defn: x |--> x + y
                                    y |--> x - y
        """
        U = repr(self._underlying).split('\n')
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
            return P({a: self._underlying(b) for a, b in x.monomial_coefficients().items()})
        except Exception:
            pass
        try:
            return P([self._underlying(b) for b in x])
        except Exception:
            pass
        try:
            return P(self._underlying(x.numerator()))/P(self._underlying(x.denominator()))
        except Exception:
            raise TypeError("invalid argument %s" % repr(x))

    @cached_method
    def inverse(self):
        """
        Return the inverse of this ring homomorphism if the underlying
        homomorphism of the base ring is invertible.

        EXAMPLES::

            sage: R.<x,y> = QQ[]
            sage: S.<a,b> = QQ[]
            sage: f = R.hom([a + b, a - b], S)
            sage: PR.<t> = R[]
            sage: PS = S['t']
            sage: Pf = PR.hom(f, PS)
            sage: Pf.inverse()                                                          # needs sage.libs.singular
            Ring morphism:
              From: Univariate Polynomial Ring in t over Multivariate
                    Polynomial Ring in a, b over Rational Field
              To:   Univariate Polynomial Ring in t over Multivariate
                    Polynomial Ring in x, y over Rational Field
              Defn: Induced from base ring by
                    Ring morphism:
                      From: Multivariate Polynomial Ring in a, b over Rational Field
                      To:   Multivariate Polynomial Ring in x, y over Rational Field
                      Defn: a |--> 1/2*x + 1/2*y
                            b |--> 1/2*x - 1/2*y
            sage: Pf.inverse()(Pf(x*t^2 + y*t))                                         # needs sage.libs.singular
            x*t^2 + y*t
        """
        return self.parent().reversed()(self._underlying.inverse())


cdef class RingHomomorphism_from_fraction_field(RingHomomorphism):
    r"""
    Morphisms between fraction fields.

    TESTS::

        sage: S.<x> = QQ[]
        sage: f = S.hom([x^2])
        sage: g = f.extend_to_fraction_field()                                          # needs sage.libs.singular
        sage: type(g)                                                                   # needs sage.libs.singular
        <class 'sage.rings.morphism.RingHomomorphism_from_fraction_field'>
    """
    def __init__(self, parent, morphism):
        r"""
        Initialize this morphism.

        TESTS::

            sage: # needs sage.rings.number_field
            sage: x = polygen(ZZ, 'x')
            sage: A.<a> = ZZ.extension(x^2 - 2)
            sage: f = A.coerce_map_from(ZZ)
            sage: g = f.extend_to_fraction_field()   # indirect doctest
            sage: g
            Ring morphism:
              From: Rational Field
              To:   Number Field in a with defining polynomial x^2 - 2
        """
        RingHomomorphism.__init__(self, parent)
        self._morphism = morphism

    def _repr_defn(self):
        r"""
        Return a string definition of this morphism.

        EXAMPLES::

            sage: S.<x> = QQ[]
            sage: f = S.hom([x^2]).extend_to_fraction_field()                           # needs sage.libs.singular
            sage: f                                                                     # needs sage.libs.singular
            Ring endomorphism of Fraction Field of Univariate Polynomial Ring in x over Rational Field
              Defn: x |--> x^2
            sage: f._repr_defn()                                                        # needs sage.libs.singular
            'x |--> x^2'
        """
        return self._morphism._repr_defn()

    cpdef Element _call_(self, x):
        r"""
        Return the value of this morphism at ``x``.

        INPUT:

        - ``x`` -- an element in the domain of this morphism

        EXAMPLES::

            sage: S.<x> = QQ[]
            sage: f = S.hom([x + 1]).extend_to_fraction_field()                         # needs sage.libs.singular
            sage: f(1/x)                                                                # needs sage.libs.singular
            1/(x + 1)
            sage: f(1/(x-1))                                                            # needs sage.libs.singular
            1/x
        """
        return self._morphism(x.numerator()) / self._morphism(x.denominator())

    cdef _update_slots(self, dict _slots):
        """
        Helper function for copying and pickling.

        TESTS::

            sage: # needs sage.libs.singular
            sage: S.<x> = QQ[]
            sage: f = S.hom([x + 1]).extend_to_fraction_field()

            sage: g = copy(f)    # indirect doctest                                     # needs sage.libs.singular
            sage: f == g                                                                # needs sage.libs.singular
            True
            sage: f is g                                                                # needs sage.libs.singular
            False
        """
        self._morphism = _slots['_morphism']
        RingHomomorphism._update_slots(self, _slots)

    cdef dict _extra_slots(self):
        """
        Helper function for copying and pickling.

        TESTS::

            sage: S.<x> = QQ[]
            sage: f = S.hom([x + 1]).extend_to_fraction_field()                         # needs sage.libs.singular
            sage: loads(dumps(f)) == f                                                  # needs sage.libs.singular
            True
        """
        slots = RingHomomorphism._extra_slots(self)
        slots['_morphism'] = self._morphism
        return slots

    @cached_method
    def inverse(self):
        """
        Return the inverse of this ring homomorphism if it exists.

        EXAMPLES::

            sage: S.<x> = QQ[]
            sage: f = S.hom([2*x - 1])
            sage: g = f.extend_to_fraction_field()                                      # needs sage.libs.singular
            sage: g.inverse()                                                           # needs sage.libs.singular
            Ring endomorphism of Fraction Field of Univariate Polynomial Ring
             in x over Rational Field
              Defn: x |--> 1/2*x + 1/2
        """
        return self.parent().reversed()(self._morphism.inverse())


cdef class RingHomomorphism_cover(RingHomomorphism):
    r"""
    A homomorphism induced by quotienting a ring out by an ideal.

    EXAMPLES::

        sage: R.<x,y> = PolynomialRing(QQ, 2)
        sage: S.<a,b> = R.quo(x^2 + y^2)                                                # needs sage.libs.singular
        sage: phi = S.cover(); phi                                                      # needs sage.libs.singular
        Ring morphism:
          From: Multivariate Polynomial Ring in x, y over Rational Field
          To:   Quotient of Multivariate Polynomial Ring in x, y over Rational Field
                by the ideal (x^2 + y^2)
          Defn: Natural quotient map
        sage: phi(x + y)                                                                # needs sage.libs.singular
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
            <class 'sage.rings.morphism.RingHomomorphism_cover'>
        """
        RingHomomorphism.__init__(self, parent)

    cpdef Element _call_(self, x):
        """
        Evaluate this covering homomorphism at ``x``, which just involves
        coercing ``x`` into the domain, then codomain.

        EXAMPLES::

            sage: f = Zmod(6).cover()
            sage: type(f)
            <class 'sage.rings.morphism.RingHomomorphism_cover'>
            sage: f(-5)                 # indirect doctest
            1

        TESTS:

        We verify that calling directly raises the expected error
        (just coercing into the codomain), but calling with __call__
        (the second call below) gives a :exc:`TypeError` since 1/2 cannot be
        coerced into the domain. ::

            sage: f._call_(1/2)
            Traceback (most recent call last):
            ...
            ZeroDivisionError: inverse of Mod(2, 6) does not exist
            sage: f(1/2)
            Traceback (most recent call last):
            ...
            TypeError: 1/2 fails to convert into the map's domain Integer Ring,
            but a `pushforward` method is not properly implemented
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
            <class 'sage.rings.morphism.RingHomomorphism_cover'>
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

            sage: # needs sage.libs.singular
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

            sage: # needs sage.libs.singular
            sage: R.<x,y> = PolynomialRing(QQ, 2)
            sage: S.<a,b> = R.quo(x^2 + y^2)
            sage: phi = S.cover()
            sage: type(phi)
            <class 'sage.rings.morphism.RingHomomorphism_cover'>
            sage: hash(phi) == hash(phi)
            True
            sage: {phi: 1}[phi]
            1
        """
        return hash((self.domain(), self.codomain()))

    def _inverse_image_ideal(self, I):
        """
        Return the inverse image of the ideal `I` under this covering morphism.

        INPUT:

        - ``I`` -- an ideal in the quotient ring

        EXAMPLES::

            sage: # needs sage.libs.singular
            sage: R.<x,y> = QQ['x,y'].quotient('x^2 * y^2')
            sage: R.cover().inverse_image(R.ideal(x^3, y^3 + 1))
            Ideal (x^2*y^2, x^3, y^3 + 1) of Multivariate Polynomial Ring
            in x, y over Rational Field
            sage: S.<u,v> = QQbar['u,v'].quotient('u^4 - 1')
            sage: S.cover().inverse_image(S.ideal(u^2 - 1))
            Ideal (u^4 - 1, u^2 - 1) of Multivariate Polynomial Ring in u, v
            over Algebraic Field
        """
        if I.is_zero():
            return self.kernel()
        return self.kernel() + [f.lift() for f in I.gens()]

    def _inverse_image_element(self, b):
        """
        Lift an element from the quotient to the cover ring of this ring
        homomorphism.

        EXAMPLES::

            sage: Q.<u,v> = QQ['x,y'].quotient('x + y')                                 # needs sage.libs.singular
            sage: Q.cover().inverse_image(u)                                            # needs sage.libs.singular
            -y
        """
        return b.lift()


cdef class RingHomomorphism_from_quotient(RingHomomorphism):
    r"""
    A ring homomorphism with domain a generic quotient ring.

    INPUT:

    - ``parent`` -- a ring homset ``Hom(R,S)``

    - ``phi`` -- a ring homomorphism ``C --> S``, where ``C`` is the
      domain of ``R.cover()``

    OUTPUT: a ring homomorphism

    The domain `R` is a quotient object `C \to R`, and
    ``R.cover()`` is the ring homomorphism
    `\varphi: C \to R`. The condition on the elements
    ``im_gens`` of `S` is that they define a
    homomorphism `C \to S` such that each generator of the
    kernel of `\varphi` maps to `0`.

    EXAMPLES::

        sage: # needs sage.libs.singular
        sage: R.<x, y, z> = PolynomialRing(QQ, 3)
        sage: S.<a, b, c> = R.quo(x^3 + y^3 + z^3)
        sage: phi = S.hom([b, c, a]); phi
        Ring endomorphism of Quotient of Multivariate Polynomial Ring in x, y, z
         over Rational Field by the ideal (x^3 + y^3 + z^3)
          Defn: a |--> b
                b |--> c
                c |--> a
        sage: phi(a + b + c)
        a + b + c
        sage: loads(dumps(phi)) == phi
        True

    Validity of the homomorphism is determined, when possible, and a
    :exc:`TypeError` is raised if there is no homomorphism sending the
    generators to the given images::

        sage: S.hom([b^2, c^2, a^2])                                                    # needs sage.libs.singular
        Traceback (most recent call last):
        ...
        ValueError: relations do not all (canonically) map to 0
        under map determined by images of generators
    """
    def __init__(self, parent, phi):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: R.<x,y> = QQ[]; S.<xx,yy> = R.quo([x^2,y^2]); S.hom([yy,xx])          # needs sage.libs.singular
            Ring endomorphism of Quotient of Multivariate Polynomial Ring in x, y
             over Rational Field by the ideal (x^2, y^2)
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

            sage: # needs sage.libs.singular
            sage: R.<x, y, z> = PolynomialRing(QQ, 3)
            sage: S.<a, b, c> = R.quo(x^3 + y^3 + z^3)
            sage: phi = S.hom([b, c, a]); phi
            Ring endomorphism of Quotient of Multivariate Polynomial Ring in x, y, z
             over Rational Field by the ideal (x^3 + y^3 + z^3)
              Defn: a |--> b
                    b |--> c
                    c |--> a
            sage: phi(a + b + c)
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

            sage: # needs sage.libs.singular
            sage: R.<x, y, z> = PolynomialRing(QQ, 3)
            sage: S.<a, b, c> = R.quo(x^3 + y^3 + z^3)
            sage: phi = S.hom([b, c, a]); phi
            Ring endomorphism of Quotient of Multivariate Polynomial Ring in x, y, z
             over Rational Field by the ideal (x^3 + y^3 + z^3)
              Defn: a |--> b
                    b |--> c
                    c |--> a
            sage: phi(a + b + c)
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

            sage: R.<x,y> = QQ[]; S.<xx,yy> = R.quo([x^2, y^2]); f = S.hom([yy,xx])     # needs sage.libs.singular
            sage: f._phi()                                                              # needs sage.libs.singular
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

            sage: R.<x,y> = QQ[]; S.<xx,yy> = R.quo([x^2, y^2])                         # needs sage.libs.singular
            sage: S.hom([yy,xx]).morphism_from_cover()                                  # needs sage.libs.singular
            Ring morphism:
              From: Multivariate Polynomial Ring in x, y over Rational Field
              To:   Quotient of Multivariate Polynomial Ring in x, y
                    over Rational Field by the ideal (x^2, y^2)
              Defn: x |--> yy
                    y |--> xx
        """
        return self.phi

    cpdef _richcmp_(self, other, int op):
        """
        Compare ``self`` to ``other``.

        EXAMPLES::

            sage: # needs sage.libs.singular
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

            sage: # needs sage.libs.singular
            sage: R.<x, y, z> = PolynomialRing(GF(19), 3)
            sage: S.<a, b, c> = R.quo(x^3 + y^3 + z^3)
            sage: phi = S.hom([b, c, a])
            sage: type(phi)
            <class 'sage.rings.morphism.RingHomomorphism_from_quotient'>
            sage: hash(phi) == hash(phi)
            True
            sage: {phi: 1}[phi]
            1
        """
        return hash(self.phi)

    def _repr_defn(self) -> str:
        """
        Used internally for printing this function.

        EXAMPLES::

            sage: R.<x,y> = QQ[]; S.<xx,yy> = R.quo([x^2,y^2])                          # needs sage.libs.singular
            sage: f = S.hom([yy, xx])                                                   # needs sage.libs.singular
            sage: print(f._repr_defn())                                                 # needs sage.libs.singular
            xx |--> yy
            yy |--> xx
        """
        D = self.domain()
        ig = self.phi.im_gens()
        return '\n'.join('{} |--> {}'.format(D.gen(i), ig[i])
                         for i in range(D.ngens()))

    cpdef Element _call_(self, x):
        """
        Evaluate this function at ``x``.

        EXAMPLES::

            sage: R.<x,y> = QQ[]; S.<xx,yy> = R.quo([x^2, y^2]); f = S.hom([yy, xx])    # needs sage.libs.singular
            sage: f(3*x + (1/2)*y)   # indirect doctest                                 # needs sage.libs.singular
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

        - ``domain`` -- a ring

        - ``n`` -- nonnegative integer (default: 1)

        OUTPUT:

        The `n`-th power of the absolute (arithmetic) Frobenius
        endomorphism on ``domain``

        TESTS::

            sage: from sage.rings.morphism import FrobeniusEndomorphism_generic
            sage: K.<u> = PowerSeriesRing(GF(5))
            sage: FrobeniusEndomorphism_generic(K)
            Frobenius endomorphism x |--> x^5 of Power Series Ring in u
             over Finite Field of size 5
            sage: FrobeniusEndomorphism_generic(K, 2)
            Frobenius endomorphism x |--> x^(5^2) of Power Series Ring in u
             over Finite Field of size 5
        """
        from sage.rings.ring import CommutativeRing
        from sage.categories.commutative_rings import CommutativeRings
        from sage.categories.homset import Hom
        if not (domain in CommutativeRings() or
                isinstance(domain, CommutativeRing)):  # TODO: remove this line
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

    cdef _update_slots(self, dict _slots):
        """
        Update information with the given slots.

        Helper function for copying or pickling.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K = Frac(GF(5)['T'])
            sage: phi = K.frobenius_endomorphism()
            sage: psi = copy(phi)
            sage: phi == psi
            True
        """
        self._p = _slots['_domain'].characteristic()
        self._power = _slots['_power']
        self._q = self._p ** self._power
        RingHomomorphism._update_slots(self, _slots)

    cdef dict _extra_slots(self):
        """
        Return additional information about this morphism
        as a dictionary.

        Helper function for copying or pickling.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K = Frac(GF(25)['T'])
            sage: phi = K.frobenius_endomorphism(2)
            sage: phi
            Frobenius endomorphism x |--> x^(5^2) of Fraction Field of
             Univariate Polynomial Ring in T over Finite Field in z2 of size 5^2
            sage: psi = loads(dumps(phi)); psi
            Frobenius endomorphism x |--> x^(5^2) of Fraction Field of
             Univariate Polynomial Ring in T over Finite Field in z2 of size 5^2
            sage: phi == psi
            True
        """
        slots = RingHomomorphism._extra_slots(self)
        slots['_power'] = self._power
        return slots

    def _repr_(self):
        """
        Return a string representation of this endomorphism.

        EXAMPLES::

            sage: K.<u> = PowerSeriesRing(GF(5))
            sage: Frob = K.frobenius_endomorphism(); Frob
            Frobenius endomorphism x |--> x^5 of Power Series Ring in u
             over Finite Field of size 5

            sage: Frob._repr_()
            'Frobenius endomorphism x |--> x^5 of Power Series Ring in u
             over Finite Field of size 5'
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

            sage: # needs sage.rings.finite_rings
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

            sage: # needs sage.rings.finite_rings
            sage: K.<u> = PowerSeriesRing(GF(5))
            sage: Frob = K.frobenius_endomorphism()
            sage: Frob(u)
            u^5
            sage: (Frob^2)(1 + u)
            1 + u^25
        """
        return x ** self._q

    def power(self):
        """
        Return an integer `n` such that this endomorphism
        is the `n`-th power of the absolute (arithmetic)
        Frobenius.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
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
            Frobenius endomorphism x |--> x^5 of Power Series Ring in u
             over Finite Field of size 5
            sage: Frob^2
            Frobenius endomorphism x |--> x^(5^2) of Power Series Ring in u
             over Finite Field of size 5
        """
        return self.__class__(self.domain(), self.power()*n)

    def _composition(self, right):
        """
        Return ``self`` o ``right``.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<u> = PowerSeriesRing(GF(5))
            sage: f = K.frobenius_endomorphism(); f
            Frobenius endomorphism x |--> x^5 of Power Series Ring in u
             over Finite Field of size 5
            sage: g = K.frobenius_endomorphism(2); g
            Frobenius endomorphism x |--> x^(5^2) of Power Series Ring in u
             over Finite Field of size 5
            sage: f * g
            Frobenius endomorphism x |--> x^(5^3) of Power Series Ring in u
             over Finite Field of size 5
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


def _tensor_product_ring(B, A):
    """
    Construct a quotient ring representing the tensor product of two rings
    over a common base ring.

    Allowed arguments are polynomial rings, quotient rings, number fields and
    finite fields.

    EXAMPLES::

        sage: # needs sage.libs.singular
        sage: from sage.rings.morphism import _tensor_product_ring
        sage: R.<x,y> = QQ[]
        sage: S.<u,v> = R.quotient(x^2 + y^2)
        sage: Q = _tensor_product_ring(S, R); Q
        Quotient of Multivariate Polynomial Ring in u, v, x, y over
        Rational Field by the ideal (u^2 + v^2)
        sage: Q.term_order()
        Block term order with blocks:
        (Degree reverse lexicographic term order of length 2,
         Degree reverse lexicographic term order of length 2)
        sage: _tensor_product_ring(R, R)
        Multivariate Polynomial Ring in y0, y1, x0, x1 over Rational Field

    TESTS:

    Local orderings are not supported::

        sage: R = PolynomialRing(QQ, 'x,y', order='negdeglex')
        sage: _tensor_product_ring(R, R)                                                # needs sage.libs.singular
        Traceback (most recent call last):
        ...
        ValueError: term ordering must be global
    """
    from sage.rings.finite_rings.finite_field_base import FiniteField
    from sage.rings.number_field.number_field_base import NumberField
    from sage.rings.polynomial.multi_polynomial_ring import MPolynomialRing_base
    from sage.rings.polynomial.polynomial_quotient_ring import PolynomialQuotientRing_generic
    from sage.rings.polynomial.polynomial_ring import PolynomialRing_generic
    from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
    from sage.rings.polynomial.term_order import TermOrder
    from sage.rings.quotient_ring import QuotientRing_nc

    if set(B.variable_names()).isdisjoint(A.variable_names()):
        names = B.variable_names() + A.variable_names()
    else:
        names = (['y%d' % d for d in range(B.ngens())] +
                 ['x%d' % d for d in range(A.ngens())])

    def term_order(A):
        # univariate rings do not have a term order
        if (isinstance(A, (PolynomialRing_generic, PolynomialQuotientRing_generic))
            or (isinstance(A, (NumberField, FiniteField))
                and not A.is_prime_field())):
            return TermOrder('lex', 1)
        try:
            t = A.term_order()
        except AttributeError:
            raise NotImplementedError("inverse not implemented for "
                                      "morphisms of %s" % A)
        if not t.is_global():
            raise ValueError("term ordering must be global")
        return t
    R = PolynomialRing(A.base_ring(), names=names,
                       order=term_order(B) + term_order(A))

    def relations(A, R_gens_A):
        if isinstance(A, (MPolynomialRing_base, PolynomialRing_generic)):
            return []
        elif isinstance(A, PolynomialQuotientRing_generic):
            to_R = A.ambient().hom(R_gens_A, R, check=False)
            return [to_R(A.modulus())]
        elif isinstance(A, QuotientRing_nc):
            to_R = A.ambient().hom(R_gens_A, R, check=False)
            return list(to_R(A.defining_ideal()).gens())
        elif (isinstance(A, (NumberField, FiniteField))
              and not A.is_prime_field()):
            to_R = A.polynomial_ring().hom(R_gens_A, R, check=False)
            return [to_R(A.polynomial())]
        else:
            raise NotImplementedError("inverse not implemented for "
                                      "morphisms of %s" % A)
    rels_A = relations(A, R.gens()[B.ngens():])
    rels_B = relations(B, R.gens()[:B.ngens()])
    return R.quotient(rels_A + rels_B, names=R.variable_names())
