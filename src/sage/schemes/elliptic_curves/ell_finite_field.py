"""
Elliptic curves over finite fields

AUTHORS:

- William Stein (2005): Initial version

- Robert Bradshaw et al....

- John Cremona (2008-02): Point counting and group structure for
  non-prime fields, Frobenius endomorphism and order, elliptic logs

- Mariah Lenox (2011-03): Added set_order method
"""

# ****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************


from sage.schemes.curves.projective_curve import Hasse_bounds
from .ell_field import EllipticCurve_field
from .constructor import EllipticCurve
from sage.schemes.hyperelliptic_curves.hyperelliptic_finite_field import HyperellipticCurve_finite_field
from sage.rings.all import Integer, ZZ, PolynomialRing, GF, polygen
from sage.rings.finite_rings.element_base import is_FiniteFieldElement
import sage.groups.generic as generic
from . import ell_point
from sage.arith.all import gcd, lcm, binomial
from sage.misc.cachefunc import cached_method
from sage.groups.additive_abelian.additive_abelian_wrapper import AdditiveAbelianGroupWrapper

import sage.plot.all as plot


class EllipticCurve_finite_field(EllipticCurve_field, HyperellipticCurve_finite_field):
    r"""
    Elliptic curve over a finite field.

    EXAMPLES::

        sage: EllipticCurve(GF(101),[2,3])
        Elliptic Curve defined by y^2  = x^3 + 2*x + 3 over Finite Field of size 101

        sage: F = GF(101^2, 'a')
        sage: EllipticCurve([F(2),F(3)])
        Elliptic Curve defined by y^2  = x^3 + 2*x + 3 over Finite Field in a of size 101^2

    Elliptic curves over `\ZZ/N\ZZ` with `N` prime are of type
    "elliptic curve over a finite field"::

        sage: F = Zmod(101)
        sage: EllipticCurve(F, [2, 3])
        Elliptic Curve defined by y^2 = x^3 + 2*x + 3 over Ring of integers modulo 101
        sage: E = EllipticCurve([F(2), F(3)])
        sage: type(E)
        <class 'sage.schemes.elliptic_curves.ell_finite_field.EllipticCurve_finite_field_with_category'>
        sage: E.category()
        Category of schemes over Ring of integers modulo 101

    Elliptic curves over `\ZZ/N\ZZ` with `N` composite are of type
    "generic elliptic curve"::

        sage: F = Zmod(95)
        sage: EllipticCurve(F, [2, 3])
        Elliptic Curve defined by y^2 = x^3 + 2*x + 3 over Ring of integers modulo 95
        sage: E = EllipticCurve([F(2), F(3)])
        sage: type(E)
        <class 'sage.schemes.elliptic_curves.ell_generic.EllipticCurve_generic_with_category'>
        sage: E.category()
        Category of schemes over Ring of integers modulo 95
        sage: TestSuite(E).run(skip=["_test_elements"])
    """

    _point = ell_point.EllipticCurvePoint_finite_field

    def plot(self, *args, **kwds):
        """
        Draw a graph of this elliptic curve over a prime finite field.

        INPUT:

        -  ``*args, **kwds`` - all other options are passed
           to the circle graphing primitive.

        EXAMPLES::

            sage: E = EllipticCurve(FiniteField(17), [0,1])
            sage: P = plot(E, rgbcolor=(0,0,1))
        """
        R = self.base_ring()
        if not R.is_prime_field():
            raise NotImplementedError

        G = plot.Graphics()
        G += plot.points([P[0:2] for P in self.points() if not P.is_zero()], *args, **kwds)

        return G

    def _points_via_group_structure(self):
        """
        Return a list of all the points on the curve, for prime fields only.

        See points() for the general case.

        EXAMPLES::

            sage: S = EllipticCurve(GF(97),[2,3])._points_via_group_structure()
            sage: len(S)
            100

        See :trac:`4687`, where the following example did not work::

            sage: E = EllipticCurve(GF(2),[0, 0, 1, 1, 1])
            sage: E.points()
            [(0 : 1 : 0)]

        ::

            sage: E = EllipticCurve(GF(2),[0, 0, 1, 0, 1])
            sage: E.points()
            [(0 : 1 : 0), (1 : 0 : 1), (1 : 1 : 1)]

        ::

            sage: E = EllipticCurve(GF(4,'a'),[0, 0, 1, 0, 1])
            sage: E.points()
            [(0 : 1 : 0), (0 : a : 1), (0 : a + 1 : 1), (1 : 0 : 1), (1 : 1 : 1), (a : 0 : 1), (a : 1 : 1), (a + 1 : 0 : 1), (a + 1 : 1 : 1)]
        """
        # TODO, eliminate when polynomial calling is fast
        G = self.abelian_group()
        pts = [x.element() for x in G.gens()]

        ni = G.generator_orders()
        ngens = G.ngens()

        H0=[self(0)]
        if ngens == 0:    # trivial group
            return H0
        for m in range(1,ni[0]):
            H0.append(H0[-1]+pts[0])
        if ngens == 1:    # cyclic group
            return H0

        # else noncyclic group
        H1=[self(0)]
        for m in range(1,ni[1]):
            H1.append(H1[-1]+pts[1])
        return [P+Q for P in H0 for Q in H1]

    def points(self):
        r"""
        All the points on this elliptic curve. The list of points is cached
        so subsequent calls are free.

        EXAMPLES::

            sage: p = 5
            sage: F = GF(p)
            sage: E = EllipticCurve(F, [1, 3])
            sage: a_sub_p = E.change_ring(QQ).ap(p); a_sub_p
            2

        ::

            sage: len(E.points())
            4
            sage: p + 1 - a_sub_p
            4
            sage: E.points()
            [(0 : 1 : 0), (1 : 0 : 1), (4 : 1 : 1), (4 : 4 : 1)]

        ::

            sage: K = GF((p, 2),'a')
            sage: E = E.change_ring(K)
            sage: len(E.points())
            32
            sage: (p + 1)**2 - a_sub_p**2
            32
            sage: w = E.points(); w
            [(0 : 1 : 0), (0 : 2*a + 4 : 1), (0 : 3*a + 1 : 1), (1 : 0 : 1), (2 : 2*a + 4 : 1), (2 : 3*a + 1 : 1), (3 : 2*a + 4 : 1), (3 : 3*a + 1 : 1), (4 : 1 : 1), (4 : 4 : 1), (a : 1 : 1), (a : 4 : 1), (a + 2 : a + 1 : 1), (a + 2 : 4*a + 4 : 1), (a + 3 : a : 1), (a + 3 : 4*a : 1), (a + 4 : 0 : 1), (2*a : 2*a : 1), (2*a : 3*a : 1), (2*a + 4 : a + 1 : 1), (2*a + 4 : 4*a + 4 : 1), (3*a + 1 : a + 3 : 1), (3*a + 1 : 4*a + 2 : 1), (3*a + 2 : 2*a + 3 : 1), (3*a + 2 : 3*a + 2 : 1), (4*a : 0 : 1), (4*a + 1 : 1 : 1), (4*a + 1 : 4 : 1), (4*a + 3 : a + 3 : 1), (4*a + 3 : 4*a + 2 : 1), (4*a + 4 : a + 4 : 1), (4*a + 4 : 4*a + 1 : 1)]

        Note that the returned list is an immutable sorted Sequence::

            sage: w[0] = 9
            Traceback (most recent call last):
            ...
            ValueError: object is immutable; please change a copy instead.
        """
        try:
            return self.__points
        except AttributeError:
            pass

        from sage.structure.sequence import Sequence
        k = self.base_ring()
        if k.is_prime_field() and k.order()>50:
            v = self._points_via_group_structure()
        else:
            v =self._points_fast_sqrt()
        v.sort()
        self.__points = Sequence(v, immutable=True)
        return self.__points

    rational_points = points

    def count_points(self, n=1):
        """
        Return the cardinality of this elliptic curve over the base field or extensions.

        INPUT:

        - ``n`` (int) -- a positive integer

        OUTPUT:

        If `n=1`, returns the cardinality of the curve over its base field.

        If `n>1`, returns a list `[c_1, c_2, ..., c_n]` where `c_d` is
        the cardinality of the curve over the extension of degree `d`
        of its base field.

        EXAMPLES::

            sage: p = 101
            sage: F = GF(p)
            sage: E = EllipticCurve(F, [2,3])
            sage: E.count_points(1)
            96
            sage: E.count_points(5)
            [96, 10368, 1031904, 104053248, 10509895776]

        ::

            sage: F.<a> = GF(p^2)
            sage: E = EllipticCurve(F, [a,a])
            sage: E.cardinality()
            10295
            sage: E.count_points()
            10295
            sage: E.count_points(1)
            10295
            sage: E.count_points(5)
            [10295, 104072155, 1061518108880, 10828567126268595, 110462212555439192375]

        """
        try:
            n = Integer(n)
        except TypeError:
            raise TypeError("n must be a positive integer")

        if n < 1:
            raise ValueError("n must be a positive integer")

        if n == 1:
            return self.cardinality()

        return [self.cardinality(extension_degree=i) for i in range(1, n + 1)]

    def random_element(self):
        """
        Return a random point on this elliptic curve, uniformly chosen
        among all rational points.

        ALGORITHM:

        Choose the point at infinity with probability `1/(2q + 1)`.
        Otherwise, take a random element from the field as x-coordinate
        and compute the possible y-coordinates. Return the i'th
        possible y-coordinate, where i is randomly chosen to be 0 or 1.
        If the i'th y-coordinate does not exist (either there is no
        point with the given x-coordinate or we hit a 2-torsion point
        with i == 1), try again.

        This gives a uniform distribution because you can imagine
        `2q + 1` buckets, one for the point at infinity and 2 for each
        element of the field (representing the x-coordinates). This
        gives a 1-to-1 map of elliptic curve points into buckets. At
        every iteration, we simply choose a random bucket until we find
        a bucket containing a point.

        AUTHOR:

        - Jeroen Demeyer (2014-09-09): choose points uniformly random,
          see :trac:`16951`.

        EXAMPLES::

            sage: k = GF(next_prime(7^5))
            sage: E = EllipticCurve(k,[2,4])
            sage: P = E.random_element(); P  # random
            (16740 : 12486 : 1)
            sage: type(P)
            <class 'sage.schemes.elliptic_curves.ell_point.EllipticCurvePoint_finite_field'>
            sage: P in E
            True

        ::

            sage: k.<a> = GF(7^5)
            sage: E = EllipticCurve(k,[2,4])
            sage: P = E.random_element(); P  # random
            (5*a^4 + 3*a^3 + 2*a^2 + a + 4 : 2*a^4 + 3*a^3 + 4*a^2 + a + 5 : 1)
            sage: type(P)
            <class 'sage.schemes.elliptic_curves.ell_point.EllipticCurvePoint_finite_field'>
            sage: P in E
            True

        ::

            sage: k.<a> = GF(2^5)
            sage: E = EllipticCurve(k,[a^2,a,1,a+1,1])
            sage: P = E.random_element();P  # random
            (a^4 + a : a^4 + a^3 + a^2 : 1)
            sage: type(P)
            <class 'sage.schemes.elliptic_curves.ell_point.EllipticCurvePoint_finite_field'>
            sage: P in E
            True

        Ensure that the entire point set is reachable::

            sage: E = EllipticCurve(GF(11), [2,1])
            sage: S = set()
            sage: while len(S) < E.cardinality():
            ....:     S.add(E.random_element())

        TESTS:

        See :trac:`8311`::

            sage: E = EllipticCurve(GF(3), [0,0,0,2,2])
            sage: E.random_element()
            (0 : 1 : 0)
            sage: E.cardinality()
            1

            sage: E = EllipticCurve(GF(2), [0,0,1,1,1])
            sage: E.random_point()
            (0 : 1 : 0)
            sage: E.cardinality()
            1

            sage: F.<a> = GF(4)
            sage: E = EllipticCurve(F, [0, 0, 1, 0, a])
            sage: E.random_point()
            (0 : 1 : 0)
            sage: E.cardinality()
            1

        """
        k = self.base_field()
        n = 2 * k.order() + 1

        while True:
            # Choose the point at infinity with probability 1/(2q + 1)
            i = ZZ.random_element(n)
            if not i:
                return self.point(0)

            v = self.lift_x(k.random_element(), all=True)
            try:
                return v[i % 2]
            except IndexError:
                pass

    random_point = random_element

    def trace_of_frobenius(self):
        r"""
        Return the trace of Frobenius acting on this elliptic curve.

        .. NOTE::

            This computes the curve cardinality, which may be
            time-consuming.

        EXAMPLES::

            sage: E = EllipticCurve(GF(101),[2,3])
            sage: E.trace_of_frobenius()
            6
            sage: E = EllipticCurve(GF(11^5,'a'),[2,5])
            sage: E.trace_of_frobenius()
            802

        The following shows that the issue from :trac:`2849` is fixed::

            sage: E = EllipticCurve(GF(3^5,'a'),[-1,-1])
            sage: E.trace_of_frobenius()
            -27
        """
        return 1 + self.base_field().order() - self.cardinality()

    def cardinality(self, algorithm=None, extension_degree=1):
        r"""
        Return the number of points on this elliptic curve.

        INPUT:

        - ``algorithm`` -- (optional) string:

          - ``'pari'`` -- use the PARI C-library function ``ellcard``.

          - ``'bsgs'`` -- use the baby-step giant-step method as
             implemented in Sage, with the Cremona-Sutherland version
             of Mestre's trick.

          - ``'exhaustive'`` -- naive point counting.

          - ``'subfield'`` -- reduce to a smaller field, provided that
            the j-invariant lies in a subfield.

          - ``'all'`` -- compute cardinality with both ``'pari'`` and
            ``'bsgs'``; return result if they agree or raise a
            ``AssertionError`` if they do not

        - ``extension_degree`` -- an integer `d` (default: 1): if the
          base field is `\GF{q}`, return the cardinality of ``self``
          over the extension `\GF{q^d}` of degree `d`.

        OUTPUT:

        The order of the group of rational points of ``self`` over its
        base field, or over an extension field of degree `d` as above.
        The result is cached.

        EXAMPLES::

            sage: EllipticCurve(GF(4, 'a'), [1,2,3,4,5]).cardinality()
            8
            sage: k.<a> = GF(3^3)
            sage: l = [a^2 + 1, 2*a^2 + 2*a + 1, a^2 + a + 1, 2, 2*a]
            sage: EllipticCurve(k,l).cardinality()
            29

        ::

            sage: l = [1, 1, 0, 2, 0]
            sage: EllipticCurve(k, l).cardinality()
            38

        An even bigger extension (which we check against Magma)::

            sage: EllipticCurve(GF(3^100, 'a'), [1,2,3,4,5]).cardinality()
            515377520732011331036459693969645888996929981504
            sage: magma.eval("Order(EllipticCurve([GF(3^100)|1,2,3,4,5]))")    # optional - magma
            '515377520732011331036459693969645888996929981504'

        ::

            sage: EllipticCurve(GF(10007), [1,2,3,4,5]).cardinality()
            10076
            sage: EllipticCurve(GF(10007), [1,2,3,4,5]).cardinality(algorithm='pari')
            10076
            sage: EllipticCurve(GF(next_prime(10**20)), [1,2,3,4,5]).cardinality()
            100000000011093199520

        The cardinality is cached::

            sage: E = EllipticCurve(GF(3^100, 'a'), [1,2,3,4,5])
            sage: E.cardinality() is E.cardinality()
            True

        The following is very fast since the curve is actually defined
        over the prime field::

            sage: k.<a> = GF(11^100)
            sage: E1 = EllipticCurve(k, [3,3])
            sage: N1 = E1.cardinality(algorithm="subfield"); N1
            137806123398222701841183371720896367762643312000384671846835266941791510341065565176497846502742959856128
            sage: E1.cardinality_pari() == N1
            True
            sage: E2 = E1.quadratic_twist()
            sage: N2 = E2.cardinality(algorithm="subfield"); N2
            137806123398222701841183371720896367762643312000384656816094284101308193849980588362304472492174093035876
            sage: E2.cardinality_pari() == N2
            True
            sage: N1 + N2 == 2*(k.cardinality() + 1)
            True

        We can count points over curves defined as a reduction::

            sage: x = polygen(QQ)
            sage: K.<w> = NumberField(x^2 + x + 1)
            sage: EK = EllipticCurve(K, [0, 0, w, 2, 1])
            sage: E = EK.base_extend(K.residue_field(2))
            sage: E
            Elliptic Curve defined by y^2 + wbar*y = x^3 + 1 over Residue field in wbar of Fractional ideal (2)
            sage: E.cardinality()
            7
            sage: E = EK.base_extend(K.residue_field(w - 1))
            sage: E.abelian_group()
            Trivial group embedded in Abelian group of points on Elliptic Curve defined by y^2 + y = x^3 + 2*x + 1 over Residue field of Fractional ideal (w - 1)

        ::

            sage: R.<x> = GF(17)[]
            sage: pol = R.irreducible_element(5)
            sage: k.<a> = R.residue_field(pol)
            sage: E = EllipticCurve(R, [1, x]).base_extend(k)
            sage: E
            Elliptic Curve defined by y^2 = x^3 + x + a over Residue field in a of Principal ideal (x^5 + x + 14) of Univariate Polynomial Ring in x over Finite Field of size 17
            sage: E.cardinality()
            1421004

        TESTS::

            sage: EllipticCurve(GF(10009), [1,2,3,4,5]).cardinality(algorithm='foobar')
            Traceback (most recent call last):
            ...
            ValueError: algorithm 'foobar' is not known

        If the cardinality has already been computed, then the ``algorithm``
        keyword is ignored::

            sage: E = EllipticCurve(GF(10007), [1,2,3,4,5])
            sage: E.cardinality(algorithm='pari')
            10076
            sage: E.cardinality(algorithm='foobar')
            10076

        Check that a bug noted at :trac:`15667` is fixed::

            sage: F.<a> = GF(3^6)
            sage: EllipticCurve([a^5 + 2*a^3 + 2*a^2 + 2*a, a^4 + a^3 + 2*a + 1]).cardinality()
            784
        """
        if extension_degree > 1:
            # A recursive call to cardinality() with
            # extension_degree=1, which will cache the cardinality, is
            # made by the call to frobenius_order() here:
            frob = self.frobenius() ** extension_degree - 1
            R = self.frobenius_order()
            if R.degree() == 1:
                return frob * frob
            else:
                return frob.norm()

        # We need manual caching (not @cached_method) since various
        # other methods refer to this _order attribute, in particular
        # self.set_order().
        try:
            return self._order
        except AttributeError:
            pass

        jpol = None
        if algorithm is None:
            # Check for j in subfield
            jpol = self.j_invariant().minimal_polynomial()
            if jpol.degree() < self.base_field().degree():
                algorithm = "subfield"
            else:
                algorithm = "pari"

        if algorithm == "pari":
            N = self.cardinality_pari()
        elif algorithm == "subfield":
            if jpol is None:
                jpol = self.j_invariant().minimal_polynomial()
            N = self._cardinality_subfield(jpol)
        elif algorithm == "bsgs":
            N = self.cardinality_bsgs()
        elif algorithm == "exhaustive":
            N = self.cardinality_exhaustive()
        elif algorithm == "all":
            N = self.cardinality_pari()
            N2 = self.cardinality_bsgs()
            if N != N2:
                raise AssertionError("cardinality with pari=%s but with bsgs=%s" % (N, N2))
        else:
            raise ValueError("algorithm {!r} is not known".format(algorithm))

        self._order = N
        return N

    from .cardinality import (cardinality_bsgs,
                              cardinality_exhaustive, _cardinality_subfield)

    order = cardinality  # alias

    def frobenius_polynomial(self):
        r"""
        Return the characteristic polynomial of Frobenius.

        The Frobenius endomorphism of the elliptic curve has quadratic
        characteristic polynomial. In most cases this is irreducible and
        defines an imaginary quadratic order; for some supersingular
        curves, Frobenius is an integer a and the polynomial is
        `(x-a)^2`.

        .. NOTE::

            This computes the curve cardinality, which may be
            time-consuming.

        EXAMPLES::

            sage: E = EllipticCurve(GF(11),[3,3])
            sage: E.frobenius_polynomial()
            x^2 - 4*x + 11

        For some supersingular curves, Frobenius is in Z and the polynomial
        is a square::

            sage: E = EllipticCurve(GF(25,'a'),[0,0,0,0,1])
            sage: E.frobenius_polynomial().factor()
            (x + 5)^2
        """
        x = polygen(ZZ)
        return x**2-self.trace_of_frobenius()*x+self.base_field().cardinality()

    def frobenius_order(self):
        r"""
        Return the quadratic order Z[phi] where phi is the Frobenius
        endomorphism of the elliptic curve.

        .. NOTE::

            This computes the curve cardinality, which may be
            time-consuming.

        EXAMPLES::

            sage: E = EllipticCurve(GF(11),[3,3])
            sage: E.frobenius_order()
            Order in Number Field in phi with defining polynomial x^2 - 4*x + 11

        For some supersingular curves, Frobenius is in Z and the Frobenius
        order is Z::

            sage: E = EllipticCurve(GF(25,'a'),[0,0,0,0,1])
            sage: R = E.frobenius_order()
            sage: R
            Order in Number Field in phi with defining polynomial x + 5
            sage: R.degree()
            1
        """
        f = self.frobenius_polynomial().factor()[0][0]
        return ZZ.extension(f,names='phi')

    def frobenius(self):
        r"""
        Return the frobenius of ``self`` as an element of a quadratic order.

        .. NOTE::

            This computes the curve cardinality, which may be
            time-consuming.

        Frobenius is only determined up to conjugacy.

        EXAMPLES::

            sage: E = EllipticCurve(GF(11),[3,3])
            sage: E.frobenius()
            phi
            sage: E.frobenius().minpoly()
            x^2 - 4*x + 11

        For some supersingular curves, Frobenius is in Z::

            sage: E = EllipticCurve(GF(25,'a'),[0,0,0,0,1])
            sage: E.frobenius()
            -5
        """
        R = self.frobenius_order()
        if R.degree() == 1:
            return self.frobenius_polynomial().roots(multiplicities=False)[0]
        else:
            return R.gen(1)

    def cardinality_pari(self):
        r"""
        Return the cardinality of ``self`` using PARI.

        This uses :pari:`ellcard`.

        EXAMPLES::

            sage: p = next_prime(10^3)
            sage: E = EllipticCurve(GF(p),[3,4])
            sage: E.cardinality_pari()
            1020
            sage: K = GF(next_prime(10^6))
            sage: E = EllipticCurve(K,[1,0,0,1,1])
            sage: E.cardinality_pari()
            999945

        Since :trac:`16931`, this now works over finite fields which
        are not prime fields::

            sage: k.<a> = GF(7^3)
            sage: E = EllipticCurve_from_j(a)
            sage: E.cardinality_pari()
            318
            sage: K.<a> = GF(3^20)
            sage: E = EllipticCurve(K,[1,0,0,1,a])
            sage: E.cardinality_pari()
            3486794310

        TESTS::

            sage: E.cardinality_pari().parent()
            Integer Ring
        """
        return Integer(self.__pari__().ellcard())

    @cached_method
    def gens(self):
        r"""
        Return points which generate the abelian group of points on
        this elliptic curve.

        The algorithm involves factoring the group order of ``self``,
        but is otherwise (randomized) polynomial-time.

        (The points returned by this function are not guaranteed to be
        the same each time, although they should remain fixed within a
        single run of Sage unless :meth:`abelian_group` is called.)

        OUTPUT: a tuple of points on the curve.

        - if the group is trivial: an empty tuple.

        - if the group is cyclic: a tuple with 1 point, a generator.

        - if the group is not cyclic: a tuple with 2 points, where the
          order of the first point equals the exponent of the group.

        .. WARNING::

            In the case of 2 generators `P` and `Q`, it is not
            guaranteed that the group is the cartesian product of the 2
            cyclic groups `\langle P \rangle` and `\langle Q \rangle`.
            In other words, the order of `Q` is not as small as possible.
            If you really need a basis (rather than just a generating set)
            of the group, use :meth:`abelian_group`.

        EXAMPLES::

            sage: E = EllipticCurve(GF(11),[2,5])
            sage: P = E.gens()[0]; P # random
            (0 : 7 : 1)
            sage: E.cardinality(), P.order()
            (10, 10)
            sage: E = EllipticCurve(GF(41),[2,5])
            sage: E.gens()  # random
            ((20 : 38 : 1), (25 : 31 : 1))
            sage: E.cardinality()
            44

        If the abelian group has been computed, return those generators
        instead::

            sage: E.abelian_group()
            Additive abelian group isomorphic to Z/22 + Z/2 embedded in Abelian group of points on Elliptic Curve defined by y^2 = x^3 + 2*x + 5 over Finite Field of size 41
            sage: ab_gens = E.abelian_group().gens()
            sage: ab_gens == E.gens()
            True
            sage: E.gens()[0].order()
            22
            sage: E.gens()[1].order()
            2

        Examples with 1 and 0 generators::

            sage: F.<a> = GF(3^6)
            sage: E = EllipticCurve([a, a+1])
            sage: pts = E.gens()
            sage: len(pts)
            1
            sage: pts[0].order() == E.cardinality()
            True
            sage: E = EllipticCurve(GF(2), [0,0,1,1,1])
            sage: E.gens()
            ()

        This works over larger finite fields where :meth:`abelian_group`
        may be too expensive::

            sage: k.<a> = GF(5^60)
            sage: E = EllipticCurve([a, a])
            sage: len(E.gens())
            2
            sage: E.cardinality()
            867361737988403547206134229616487867594472
            sage: a = E.gens()[0].order(); a # random
            433680868994201773603067114808243933797236
            sage: b = E.gens()[1].order(); b # random
            30977204928157269543076222486303138128374
            sage: lcm(a,b)
            433680868994201773603067114808243933797236
        """
        card, ords, pts = self.__pari__().ellgroup(flag=1)
        if not hasattr(self, '_order'):
            self._order = ZZ(card)
        pts = tuple(self.point(list(P)) for P in pts)
        if len(pts) >= 1:
            pts[0]._order = ZZ(ords[0]) # PARI documentation: "P is of order d_1"
        return pts

    def __iter__(self):
        """
        Return an iterator through the points of this elliptic curve.

        EXAMPLES::

            sage: E = EllipticCurve(GF(11), [1,2])
            sage: for P in E:  print("{} {}".format(P, P.order()))
            (0 : 1 : 0) 1
            (1 : 2 : 1) 4
            (1 : 9 : 1) 4
            (2 : 1 : 1) 8
            ...
            (10 : 0 : 1) 2
        """
        for P in self.points():
            yield P

    def __getitem__(self, n):
        """
        Return the n'th point in self's __points list.

        This enables users to iterate over the curve's point set.

        EXAMPLES::

            sage: E = EllipticCurve(GF(97),[2,3])
            sage: S = E.points()
            sage: E[10]
            (10 : 76 : 1)
            sage: E[15]
            (17 : 10 : 1)
            sage: for P in E: print(P.order())
            1
            50
            50
            50
            50
            5
            5
            50
            ...
        """
        return self.points()[n]

    @cached_method
    def abelian_group(self):
        r"""
        Return the abelian group structure of the group of points on this
        elliptic curve.

        .. SEEALSO::

            If you do not need the complete abelian group structure but
            only generators of the group, use :meth:`gens` which can
            be much faster in some cases.

        This method relies on :meth:`gens`, which uses random points on the
        curve and hence the generators are likely to differ from one run to
        another. However, the group is cached, so the generators will not
        change in any one run of Sage.

        OUTPUT:

        - an :class:`AdditiveAbelianGroupWrapper` object encapsulating the
          abelian group of rational points on this elliptic curve

        ALGORITHM:

        We first call :meth:`gens` to obtain a generating set `(P,Q)`.
        Letting `P` denote the point of larger order `n_1`, we extend `P`
        to a basis `(P,Q')` by computing a scalar `x` such that `Q'=Q-[x]P`
        has order `n_2=\#E/n_1`. Finding `x` involves a (typically easy)
        discrete-logarithm computation.

        The complexity of the algorithm is the cost of factoring the group
        order, plus `\Theta(\sqrt{\ell})` for each prime `\ell` such that
        the rational `\ell^\infty`-torsion of ``self`` is isomorphic to
        `\ZZ/\ell^r\times\ZZ/\ell^s` with `r>s>0`, times a polynomial in
        the logarithm of the base-field size.

        AUTHORS:

        - John Cremona: original implementation
        - Lorenz Panny (2021): current implementation

        EXAMPLES::

            sage: E = EllipticCurve(GF(11),[2,5])
            sage: E.abelian_group()
            Additive abelian group isomorphic to Z/10 embedded in Abelian group of points on Elliptic Curve defined by y^2 = x^3 + 2*x + 5 over Finite Field of size 11

        ::

            sage: E = EllipticCurve(GF(41),[2,5])
            sage: E.abelian_group()
            Additive abelian group isomorphic to Z/22 + Z/2 ...

        ::

            sage: F.<a> = GF(3^6,'a')
            sage: E = EllipticCurve([a^4 + a^3 + 2*a^2 + 2*a, 2*a^5 + 2*a^3 + 2*a^2 + 1])
            sage: E.abelian_group()
            Additive abelian group isomorphic to Z/26 + Z/26 ...

        ::

            sage: F.<a> = GF(101^3,'a')
            sage: E = EllipticCurve([2*a^2 + 48*a + 27, 89*a^2 + 76*a + 24])
            sage: E.abelian_group()
            Additive abelian group isomorphic to Z/1031352 ...

        The group can be trivial::

            sage: E = EllipticCurve(GF(2),[0,0,1,1,1])
            sage: E.abelian_group()
            Trivial group embedded in Abelian group of points on Elliptic Curve defined by y^2 + y = x^3 + x + 1 over Finite Field of size 2

        Of course, there are plenty of points if we extend the field::

            sage: E.cardinality(extension_degree=100)
            1267650600228231653296516890625

        This tests the patch for :trac:`3111`, using 10 primes randomly
        selected::

            sage: E = EllipticCurve('389a')
            sage: for p in [5927, 2297, 1571, 1709, 3851, 127, 3253, 5783, 3499, 4817]:
            ....:     G = E.change_ring(GF(p)).abelian_group()
            sage: for p in prime_range(10000):  # long time (19s on sage.math, 2011)
            ....:     if p != 389:
            ....:         G = E.change_ring(GF(p)).abelian_group()

        This tests that the bug reported in :trac:`3926` has been fixed::

            sage: K.<i> = QuadraticField(-1)
            sage: OK = K.ring_of_integers()
            sage: P=K.factor(10007)[0][0]
            sage: OKmodP = OK.residue_field(P)
            sage: E = EllipticCurve([0,0,0,i,i+3])
            sage: Emod = E.change_ring(OKmodP); Emod
            Elliptic Curve defined by y^2  = x^3 + ibar*x + (ibar+3) over Residue field in ibar of Fractional ideal (10007)
            sage: Emod.abelian_group() #random generators
            (Multiplicative Abelian group isomorphic to C50067594 x C2,
            ((3152*ibar + 7679 : 7330*ibar + 7913 : 1), (8466*ibar + 1770 : 0 : 1)))
        """

        gens = self.gens()
        assert len(gens) <= 2

        if len(gens) == 2:

            P, Q = gens
            n = self.cardinality()      # cached
            n1 = P.order()              # cached
            n2 = n//n1
            assert not n1 * Q           # PARI should guarantee this

            k = n1.prime_to_m_part(n2)
            Q *= k                      # don't need; kill that part
            nQ = n2 * generic.order_from_multiple(n2*Q, n1//k//n2)

            S = n//nQ * P
            T = n2 * Q
            S.set_order(nQ//n2)         # for .discrete_log()
            x = S.discrete_log(T)
            Q -= x * n1//nQ * P

            Q.set_order(n2)             # verifies n2*Q == 0
            gens = P, Q

        orders = [T.order() for T in gens]  # cached

        self.gens.set_cache(gens)
        return AdditiveAbelianGroupWrapper(self.point_homset(), gens, orders)

    def is_isogenous(self, other, field=None, proof=True):
        """
        Return whether or not self is isogenous to other.

        INPUT:

        - ``other`` -- another elliptic curve.

        - ``field`` (default None) -- a field containing the base
          fields of the two elliptic curves into which the two curves
          may be extended to test if they are isogenous over this
          field. By default is_isogenous will not try to find this
          field unless one of the curves can be extended into the base
          field of the other, in which case it will test over the
          larger base field.

        - ``proof`` (default True) -- this parameter is here only to
          be consistent with versions for other types of elliptic
          curves.

        OUTPUT:

        (bool) True if there is an isogeny from curve ``self`` to
        curve ``other`` defined over ``field``.

        EXAMPLES::

            sage: E1 = EllipticCurve(GF(11^2,'a'),[2,7]); E1
            Elliptic Curve defined by y^2 = x^3 + 2*x + 7 over Finite Field in a of size 11^2
            sage: E1.is_isogenous(5)
            Traceback (most recent call last):
            ...
            ValueError: Second argument is not an Elliptic Curve.
            sage: E1.is_isogenous(E1)
            True

            sage: E2 = EllipticCurve(GF(7^3,'b'),[3,1]); E2
            Elliptic Curve defined by y^2 = x^3 + 3*x + 1 over Finite Field in b of size 7^3
            sage: E1.is_isogenous(E2)
            Traceback (most recent call last):
            ...
            ValueError: The base fields must have the same characteristic.

            sage: E3 = EllipticCurve(GF(11^2,'c'),[4,3]); E3
            Elliptic Curve defined by y^2 = x^3 + 4*x + 3 over Finite Field in c of size 11^2
            sage: E1.is_isogenous(E3)
            False

            sage: E4 = EllipticCurve(GF(11^6,'d'),[6,5]); E4
            Elliptic Curve defined by y^2 = x^3 + 6*x + 5 over Finite Field in d of size 11^6
            sage: E1.is_isogenous(E4)
            True

            sage: E5 = EllipticCurve(GF(11^7,'e'),[4,2]); E5
            Elliptic Curve defined by y^2 = x^3 + 4*x + 2 over Finite Field in e of size 11^7
            sage: E1.is_isogenous(E5)
            Traceback (most recent call last):
            ...
            ValueError: Curves have different base fields: use the field parameter.

        When the field is given:

            sage: E1 = EllipticCurve(GF(13^2,'a'),[2,7]); E1
            Elliptic Curve defined by y^2 = x^3 + 2*x + 7 over Finite Field in a of size 13^2
            sage: E1.is_isogenous(5,GF(13^6,'f'))
            Traceback (most recent call last):
            ...
            ValueError: Second argument is not an Elliptic Curve.
            sage: E6 = EllipticCurve(GF(11^3,'g'),[9,3]); E6
            Elliptic Curve defined by y^2 = x^3 + 9*x + 3 over Finite Field in g of size 11^3
            sage: E1.is_isogenous(E6,QQ)
            Traceback (most recent call last):
            ...
            ValueError: The base fields must have the same characteristic.
            sage: E7 = EllipticCurve(GF(13^5,'h'),[2,9]); E7
            Elliptic Curve defined by y^2 = x^3 + 2*x + 9 over Finite Field in h of size 13^5
            sage: E1.is_isogenous(E7,GF(13^4,'i'))
            Traceback (most recent call last):
            ...
            ValueError: Field must be an extension of the base fields of both curves
            sage: E1.is_isogenous(E7,GF(13^10,'j'))
            False
            sage: E1.is_isogenous(E7,GF(13^30,'j'))
            False
        """
        from .ell_generic import is_EllipticCurve
        if not is_EllipticCurve(other):
            raise ValueError("Second argument is not an Elliptic Curve.")
        if self.is_isomorphic(other):
            return True
        if self.base_field().characteristic() != other.base_field().characteristic():
            raise ValueError("The base fields must have the same characteristic.")
        if field is None:
            if self.base_field().degree() == other.base_field().degree():
                return self.cardinality() == other.cardinality()

            elif self.base_field().degree() == gcd(self.base_field().degree(),
                                                   other.base_field().degree()):
                return self.cardinality(extension_degree=other.base_field().degree()//self.base_field().degree()) == other.cardinality()

            elif other.base_field().degree() == gcd(self.base_field().degree(),
                                                    other.base_field().degree()):
                return other.cardinality(extension_degree=self.base_field().degree()//other.base_field().degree()) == self.cardinality()

            else:
                raise ValueError("Curves have different base fields: use the field parameter.")
        else:
            f_deg = field.degree()
            s_deg = self.base_field().degree()
            o_deg = other.base_field().degree()
            if not lcm(s_deg, o_deg).divides(f_deg):
                raise ValueError("Field must be an extension of the base fields of both curves")
            else:
                sc = self.cardinality(extension_degree=f_deg // s_deg)
                oc = other.cardinality(extension_degree=f_deg // o_deg)
                return sc == oc

    def is_supersingular(self, proof=True):
        r"""
        Return True if this elliptic curve is supersingular, else False.

        INPUT:

        - ``proof`` (boolean, default True) -- If True, returns a
          proved result.  If False, then a return value of False is
          certain but a return value of True may be based on a
          probabilistic test.  See the documentation of the function
          :meth:`is_j_supersingular` for more details.

        EXAMPLES::

            sage: F = GF(101)
            sage: EllipticCurve(j=F(0)).is_supersingular()
            True
            sage: EllipticCurve(j=F(1728)).is_supersingular()
            False
            sage: EllipticCurve(j=F(66)).is_supersingular()
            True
            sage: EllipticCurve(j=F(99)).is_supersingular()
            False

        TESTS::

            sage: from sage.schemes.elliptic_curves.ell_finite_field import supersingular_j_polynomial, is_j_supersingular
            sage: F = GF(103)
            sage: ssjlist = [F(1728)] + supersingular_j_polynomial(103).roots(multiplicities=False)
            sage: Set([j for j in F if is_j_supersingular(j)]) == Set(ssjlist)
            True

        """
        return is_j_supersingular(self.j_invariant(), proof=proof)

    def is_ordinary(self, proof=True):
        r"""
        Return True if this elliptic curve is ordinary, else False.

        INPUT:

        - ``proof`` (boolean, default True) -- If True, returns a
          proved result.  If False, then a return value of True is
          certain but a return value of False may be based on a
          probabilistic test.  See the documentation of the function
          :meth:`is_j_supersingular` for more details.

        EXAMPLES::

            sage: F = GF(101)
            sage: EllipticCurve(j=F(0)).is_ordinary()
            False
            sage: EllipticCurve(j=F(1728)).is_ordinary()
            True
            sage: EllipticCurve(j=F(66)).is_ordinary()
            False
            sage: EllipticCurve(j=F(99)).is_ordinary()
            True

        """
        return not is_j_supersingular(self.j_invariant(), proof=proof)

    def set_order(self, value, num_checks=8):
        r"""
        Set the value of self._order to value.

        Use this when you know a priori the order of the curve to
        avoid a potentially expensive order calculation.

        INPUT:

        - ``value`` - Integer in the Hasse-Weil range for this
          curve.

        - ``num_checks`` - Integer (default: 8) number of times to
          check whether value*(a random point on this curve) is
          equal to the identity.


        OUTPUT:

        None

        EXAMPLES:

        This example illustrates basic usage.

        ::

            sage: E = EllipticCurve(GF(7), [0, 1]) # This curve has order 6
            sage: E.set_order(6)
            sage: E.order()
            6
            sage: E.order() * E.random_point()
            (0 : 1 : 0)

        We now give a more interesting case, the NIST-P521 curve. Its
        order is too big to calculate with Sage, and takes a long time
        using other packages, so it is very useful here.

        ::

            sage: p = 2^521 - 1
            sage: prev_proof_state = proof.arithmetic()
            sage: proof.arithmetic(False) # turn off primality checking
            sage: F = GF(p)
            sage: A = p - 3
            sage: B = 1093849038073734274511112390766805569936207598951683748994586394495953116150735016013708737573759623248592132296706313309438452531591012912142327488478985984
            sage: q = 6864797660130609714981900799081393217269435300143305409394463459185543183397655394245057746333217197532963996371363321113864768612440380340372808892707005449
            sage: E = EllipticCurve([F(A), F(B)])
            sage: E.set_order(q)
            sage: G = E.random_point()
            sage: G.order() * G  # This takes practically no time.
            (0 : 1 : 0)
            sage: proof.arithmetic(prev_proof_state) # restore state

        It is an error to pass a value which is not an integer in the
        Hasse-Weil range::

            sage: E = EllipticCurve(GF(7), [0, 1]) # This curve has order 6
            sage: E.set_order("hi")
            Traceback (most recent call last):
            ...
            TypeError: unable to convert 'hi' to an integer
            sage: E.set_order(0)
            Traceback (most recent call last):
            ...
            ValueError: Value 0 illegal (not an integer in the Hasse range)
            sage: E.set_order(1000)
            Traceback (most recent call last):
            ...
            ValueError: Value 1000 illegal (not an integer in the Hasse range)

        It is also very likely an error to pass a value which is not
        the actual order of this curve. How unlikely is determined by
        num_checks, the factorization of the actual order, and the
        actual group structure::

            sage: E = EllipticCurve(GF(7), [0, 1]) # This curve has order 6
            sage: E.set_order(11)
            Traceback (most recent call last):
            ...
            ValueError: Value 11 illegal (multiple of random point not the identity)

        However, set_order can be fooled, though it's not likely in
        "real cases of interest". For instance, the order can be set
        to a multiple of the actual order::

            sage: E = EllipticCurve(GF(7), [0, 1]) # This curve has order 6
            sage: E.set_order(12)  # 12 just fits in the Hasse range
            sage: E.order()
            12

        Or, the order can be set incorrectly along with num_checks set
        too small::

            sage: E = EllipticCurve(GF(7), [0, 1]) # This curve has order 6
            sage: E.set_order(4, num_checks=0)
            sage: E.order()
            4

        The value of num_checks must be an integer. Negative values
        are interpreted as zero, which means don't do any checking::

            sage: E = EllipticCurve(GF(7), [0, 1]) # This curve has order 6
            sage: E.set_order(4, num_checks=-12)
            sage: E.order()
            4

        AUTHORS:

         - Mariah Lenox (2011-02-16)
        """
        value = Integer(value)

        # Is value in the Hasse range?
        q = self.base_field().order()
        a,b = Hasse_bounds(q,1)
        if not a <= value <= b:
            raise ValueError('Value %s illegal (not an integer in the Hasse range)' % value)
        # Is value*random == identity?
        for i in range(num_checks):
            G = self.random_point()
            if value * G != self(0):
                raise ValueError('Value %s illegal (multiple of random point not the identity)' % value)
        self._order = value

# dict to hold precomputed coefficient vectors of supersingular j values (excluding 0, 1728):

supersingular_j_polynomials = {}

def fill_ss_j_dict():
    r"""
    Fill the global cache of supersingular j-_polynomials.

    This function does nothing except the first time it is called,
    when it fills ``supersingular_j_polynomials`` with precomputed
    values for `p<300`.  Setting the values this way avoids start-up
    costs.

    """
    global supersingular_j_polynomials
    if not supersingular_j_polynomials:
        supersingular_j_polynomials[13] = [8, 1]
        supersingular_j_polynomials[17] = [9, 1]
        supersingular_j_polynomials[19] = [12, 1]
        supersingular_j_polynomials[23] = [4, 1]
        supersingular_j_polynomials[29] = [21, 2, 1]
        supersingular_j_polynomials[31] = [8, 25, 1]
        supersingular_j_polynomials[37] = [11, 5, 23, 1]
        supersingular_j_polynomials[41] = [18, 10, 19, 1]
        supersingular_j_polynomials[43] = [32, 11, 21, 1]
        supersingular_j_polynomials[47] = [35, 33, 31, 1]
        supersingular_j_polynomials[53] = [24, 9, 30, 7, 1]
        supersingular_j_polynomials[59] = [39, 31, 35, 39, 1]
        supersingular_j_polynomials[61] = [60, 21, 27, 8, 60, 1]
        supersingular_j_polynomials[67] = [8, 36, 47, 4, 53, 1]
        supersingular_j_polynomials[71] = [18, 54, 28, 33, 1, 1]
        supersingular_j_polynomials[73] = [7, 39, 38, 9, 68, 60, 1]
        supersingular_j_polynomials[79] = [10, 25, 1, 63, 57, 55, 1]
        supersingular_j_polynomials[83] = [43, 72, 81, 81, 62, 11, 1]
        supersingular_j_polynomials[89] = [42, 79, 23, 22, 37, 86, 60, 1]
        supersingular_j_polynomials[97] = [19, 28, 3, 72, 2, 96, 10, 60, 1]
        supersingular_j_polynomials[101] = [9, 76, 45, 79, 1, 68, 87, 60, 1]
        supersingular_j_polynomials[103] = [64, 15, 24, 58, 70, 83, 84, 100, 1]
        supersingular_j_polynomials[107] = [6, 18, 72, 59, 43, 19, 17, 68, 1]
        supersingular_j_polynomials[109] = [107, 22, 39, 83, 30, 34, 108, 104, 60, 1]
        supersingular_j_polynomials[113] = [86, 71, 75, 6, 47, 97, 100, 4, 60, 1]
        supersingular_j_polynomials[127] = [32, 31, 5, 50, 115, 122, 114, 67, 38, 35, 1]
        supersingular_j_polynomials[131] = [65, 64, 10, 34, 129, 35, 94, 127, 7, 7, 1]
        supersingular_j_polynomials[137] = [104, 83, 3, 82, 112, 23, 77, 135, 18, 50, 60, 1]
        supersingular_j_polynomials[139] = [87, 79, 109, 21, 138, 9, 104, 130, 61, 118, 90, 1]
        supersingular_j_polynomials[149] = [135, 55, 80, 86, 87, 74, 32, 60, 130, 80, 146, 60, 1]
        supersingular_j_polynomials[151] = [94, 125, 8, 6, 93, 21, 114, 80, 107, 58, 42, 18, 1]
        supersingular_j_polynomials[157] = [14, 95, 22, 58, 110, 23, 71, 51, 47, 5, 147, 59, 60, 1]
        supersingular_j_polynomials[163] = [102, 26, 74, 95, 112, 151, 98, 107, 27, 37, 25, 111, 109, 1]
        supersingular_j_polynomials[167] = [14, 9, 27, 109, 97, 55, 51, 74, 145, 125, 36, 113, 89, 1]
        supersingular_j_polynomials[173] = [152, 73, 56, 12, 18, 96, 98, 49, 30, 43, 52, 79, 163, 60, 1]
        supersingular_j_polynomials[179] = [110, 51, 3, 94, 123, 90, 156, 90, 88, 119, 158, 27, 71, 29, 1]
        supersingular_j_polynomials[181] = [7, 65, 77, 29, 139, 34, 65, 84, 164, 73, 51, 136, 7, 141, 60, 1]
        supersingular_j_polynomials[191] = [173, 140, 144, 3, 135, 80, 182, 84, 93, 75, 83, 17, 22, 42, 160, 1]
        supersingular_j_polynomials[193] = [23, 48, 26, 15, 108, 141, 124, 44, 132, 49, 72, 173, 126, 101, 22, 60, 1]
        supersingular_j_polynomials[197] = [14, 111, 64, 170, 193, 32, 124, 91, 112, 163, 14, 112, 167, 191, 183, 60, 1]
        supersingular_j_polynomials[199] = [125, 72, 65, 30, 63, 45, 10, 177, 91, 102, 28, 27, 5, 150, 51, 128, 1]
        supersingular_j_polynomials[211] = [27, 137, 128, 90, 102, 141, 5, 77, 131, 144, 83, 108, 23, 105, 98, 13, 80, 1]
        supersingular_j_polynomials[223] = [56, 183, 46, 133, 191, 94, 20, 8, 92, 100, 57, 200, 166, 67, 59, 218, 28, 32, 1]
        supersingular_j_polynomials[227] = [79, 192, 142, 66, 11, 114, 100, 208, 57, 147, 32, 5, 144, 93, 185, 147, 92, 16, 1]
        supersingular_j_polynomials[229] = [22, 55, 182, 130, 228, 172, 63, 25, 108, 99, 100, 101, 220, 111, 205, 199, 91, 163, 60, 1]
        supersingular_j_polynomials[233] = [101, 148, 85, 113, 226, 68, 71, 103, 61, 44, 173, 175, 5, 225, 227, 99, 146, 170, 60, 1]
        supersingular_j_polynomials[239] = [225, 81, 47, 26, 133, 182, 238, 2, 144, 154, 234, 178, 165, 130, 35, 61, 144, 112, 207, 1]
        supersingular_j_polynomials[241] = [224, 51, 227, 139, 134, 186, 187, 152, 161, 175, 213, 59, 105, 88, 87, 124, 202, 40, 15, 60, 1]
        supersingular_j_polynomials[251] = [30, 183, 80, 127, 40, 56, 230, 168, 192, 48, 226, 61, 214, 54, 165, 147, 105, 88, 38, 171, 1]
        supersingular_j_polynomials[257] = [148, 201, 140, 146, 169, 147, 220, 4, 205, 224, 35, 42, 198, 97, 127, 7, 110, 229, 118, 202, 60, 1]
        supersingular_j_polynomials[263] = [245, 126, 72, 213, 14, 64, 152, 83, 169, 114, 9, 128, 138, 231, 103, 85, 114, 211, 173, 249, 135, 1]
        supersingular_j_polynomials[269] = [159, 32, 69, 95, 201, 266, 190, 176, 76, 151, 212, 21, 106, 49, 263, 105, 136, 194, 215, 181, 237, 60, 1]
        supersingular_j_polynomials[271] = [169, 87, 179, 109, 133, 101, 31, 167, 208, 99, 127, 120, 83, 62, 36, 23, 61, 50, 69, 263, 265, 111, 1]
        supersingular_j_polynomials[277] = [251, 254, 171, 72, 190, 237, 12, 231, 123, 217, 263, 151, 270, 183, 29, 228, 85, 4, 67, 101, 29, 169, 60, 1]
        supersingular_j_polynomials[281] = [230, 15, 146, 69, 41, 23, 142, 232, 18, 80, 58, 134, 270, 62, 272, 70, 247, 189, 118, 255, 274, 159, 60, 1]
        supersingular_j_polynomials[283] = [212, 4, 42, 155, 38, 1, 270, 175, 172, 256, 264, 232, 50, 82, 244, 127, 148, 46, 249, 72, 59, 124, 75, 1]
        supersingular_j_polynomials[293] = [264, 66, 165, 144, 243, 25, 163, 210, 18, 107, 160, 153, 70, 255, 91, 211, 22, 7, 256, 50, 150, 94, 225, 60, 1]

def supersingular_j_polynomial(p, use_cache=True):
    r"""
    Return a polynomial whose roots are the supersingular
    `j`-invariants in characteristic `p`, other than 0, 1728.

    INPUT:

    - `p` (integer) -- a prime number.

    - `use_cache` (boolean, default ``True``) -- use cached coefficients if they exist

    ALGORITHM:

    First compute H(X) whose roots are the Legendre
    `\lambda`-invariants of supersingular curves (Silverman V.4.1(b))
    in characteristic `p`.  Then, using a resultant computation with
    the polynomial relating `\lambda` and `j` (Silverman III.1.7(b)),
    we recover the polynomial (in variable ``j``) whose roots are the
    `j`-invariants.  Factors of `j` and `j-1728` are removed if
    present.

    .. note::

        The only point of the use_cache parameter is to allow checking
        the precomputed coefficients.

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.ell_finite_field import supersingular_j_polynomial
        sage: f = supersingular_j_polynomial(67); f
        j^5 + 53*j^4 + 4*j^3 + 47*j^2 + 36*j + 8
        sage: f.factor()
        (j + 1) * (j^2 + 8*j + 45) * (j^2 + 44*j + 24)

    ::

        sage: [supersingular_j_polynomial(p) for p in prime_range(30)]
        [1, 1, 1, 1, 1, j + 8, j + 9, j + 12, j + 4, j^2 + 2*j + 21]

    TESTS::

        sage: from sage.schemes.elliptic_curves.ell_finite_field import supersingular_j_polynomial
        sage: supersingular_j_polynomial(6)
        Traceback (most recent call last):
        ...
        ValueError: p (=6) should be a prime number

    Check the cached values are correct::

        sage: from sage.schemes.elliptic_curves.ell_finite_field import supersingular_j_polynomial as ssjpol
        sage: assert all(ssjpol(p,True) == ssjpol(p,False) for p in primes(300))

    """
    try:
        p = ZZ(p)
    except TypeError:
        raise ValueError("p (=%s) should be a prime number" % p)
    if not p.is_prime():
        raise ValueError("p (=%s) should be a prime number" % p)

    J = polygen(GF(p),'j')
    if p<13:
        return J.parent().one()
    if use_cache:
        fill_ss_j_dict()
        if p in supersingular_j_polynomials:
            return J.parent()(supersingular_j_polynomials[p])

    from sage.misc.misc_c import prod
    m=(p-1)//2
    X,T = PolynomialRing(GF(p),2,names=['X','T']).gens()
    H = sum(binomial(m, i) ** 2 * T ** i for i in range(m + 1))
    F = T**2 * (T-1)**2 * X - 256*(T**2-T+1)**3
    R = F.resultant(H, T)
    R = prod([fi for fi, e in R([J, 0]).factor()])
    if R(0) == 0:
        R = R // J
    if R(1728) == 0:
        R = R // (J - 1728)
    supersingular_j_polynomials[p] = R.coefficients(sparse=False)
    return R

def is_j_supersingular(j, proof=True):
    r"""
    Return True if `j` is a supersingular `j`-invariant.

    INPUT:

    - ``j`` (finite field element) -- an element of a finite field

    - ``proof`` (boolean, default True) -- If True, returns a proved
      result.  If False, then a return value of False is certain but a
      return value of True may be based on a probabilistic test.  See
      the ALGORITHM section below for more details.

    OUTPUT:

    (boolean) True if `j` is supersingular, else False.

    ALGORITHM:

    For small characteristics `p` we check whether the `j`-invariant
    is in a precomputed list of supersingular values.  Otherwise we
    next check the `j`-invariant.  If `j=0`, the curve is
    supersingular if and only if `p=2` or `p\equiv3\pmod{4}`; if
    `j=1728`, the curve is supersingular if and only if `p=3` or
    `p\equiv2\pmod{3}`.  Next, if the base field is the prime field
    `{\rm GF}(p)`, we check that `(p+1)P=0` for several random points
    `P`, returning False if any fail: supersingular curves over `{\rm
    GF}(p)` have cardinality `p+1`.  If Proof is false we now return
    True.  Otherwise we compute the cardinality and return True if and
    only if it is divisible by `p`.

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.ell_finite_field import is_j_supersingular, supersingular_j_polynomials
        sage: [(p,[j for j in GF(p) if is_j_supersingular(j)]) for p in prime_range(30)]
        [(2, [0]), (3, [0]), (5, [0]), (7, [6]), (11, [0, 1]), (13, [5]), (17, [0, 8]), (19, [7, 18]), (23, [0, 3, 19]), (29, [0, 2, 25])]

        sage: [j for j in GF(109) if is_j_supersingular(j)]
        [17, 41, 43]
        sage: PolynomialRing(GF(109),'j')(supersingular_j_polynomials[109]).roots()
        [(43, 1), (41, 1), (17, 1)]

        sage: [p for p in prime_range(100) if is_j_supersingular(GF(p)(0))]
        [2, 3, 5, 11, 17, 23, 29, 41, 47, 53, 59, 71, 83, 89]
        sage: [p for p in prime_range(100) if is_j_supersingular(GF(p)(1728))]
        [2, 3, 7, 11, 19, 23, 31, 43, 47, 59, 67, 71, 79, 83]
        sage: [p for p in prime_range(100) if is_j_supersingular(GF(p)(123456))]
        [2, 3, 59, 89]

    """
    if not is_FiniteFieldElement(j):
        raise ValueError("%s must be an element of a finite field" % j)

    F = j.parent()
    p = F.characteristic()
    d = F.degree()

    if j.is_zero():
        return p == 3 or p % 3 == 2

    if (j - 1728).is_zero():
        return p == 2 or p % 4 == 3

    # From now on we know that j != 0, 1728

    if p in (2, 3, 5, 7, 11):
        return False  # since j=0, 1728 are the only s.s. invariants

    # supersingular j-invariants have degree at most 2:

    jpol = j.minimal_polynomial()
    degj = jpol.degree()
    if degj > 2:
        return False

    # if p occurs in the precomputed list, use that:

    fill_ss_j_dict()
    if p in supersingular_j_polynomials:
        return supersingular_j_polynomial(p)(j).is_zero()

    # Over GF(p), supersingular elliptic curves have cardinality
    # exactly p+1, so we check some random points in order to detect
    # non-supersingularity.  Over GF(p^2) (for p at least 5) the
    # cardinality is either (p-1)^2 or (p+1)^2, and the group has
    # exponent p+1 or p-1, so we can do a similar random check: unless
    # (p+1)*P=0 for all the random points, or (p-1)*P=0 for all of
    # them, we can certainly return False.

    # First we replace j by an element of GF(p) or GF(p^2) (since F
    # might be a proper extension of these):

    if degj == 1:
        j = -jpol(0)  # = j, but in GF(p)
    elif d > 2:
        F = GF((p, 2), 'a')
        j = jpol.roots(F, multiplicities=False)[0]  # j, but in GF(p^2)

    E = EllipticCurve(j=j)
    if degj == 1:
        for i in range(10):
            P = E.random_element()
            if not ((p + 1) * P).is_zero():
                return False
    else:
        n = None  # will hold either p+1 or p-1 later
        for i in range(10):
            P = E.random_element()
            # avoid 2-torsion;  we know that a1=a3=0 and #E>4!
            while P[2].is_zero() or P[1].is_zero():
                P = E.random_element()

            if n is None:  # not yet decided between p+1 and p-1
                pP = p*P
                if pP[0] != P[0]:  # i.e. pP is neither P nor -P
                    return False
                if pP[1] == P[1]:  # then p*P == P != -P
                    n = p - 1
                else:           # then p*P == -P != P
                    n = p + 1
            else:
                if not (n*P).is_zero():
                    return False

    # when proof is False we return True for any curve which passes
    # the probabilistic test:

    if not proof:
        return True

    # otherwise we check the trace of Frobenius (which could be
    # expensive since it involves counting the number of points on E):

    return E.trace_of_frobenius() % p == 0
