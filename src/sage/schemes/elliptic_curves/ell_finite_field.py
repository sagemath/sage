"""
Elliptic curves over finite fields

AUTHORS:

- William Stein (2005): Initial version

- Robert Bradshaw et al....

- John Cremona (2008-02): Point counting and group structure for
  non-prime fields, Frobenius endomorphism and order, elliptic logs

- Mariah Lenox (2011-03): Added ``set_order`` method

- Lorenz Panny, John Cremona (2023-02): ``.twists()``

- Lorenz Panny (2023): ``special_supersingular_curve()``

- Martin Grenouilloux (2024): ``EllipticCurve_with_prime_order()``
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

import sage.groups.generic as generic

from sage.arith.functions import lcm
from sage.arith.misc import binomial, GCD as gcd
from sage.groups.additive_abelian.additive_abelian_wrapper import AdditiveAbelianGroupWrapper
from sage.misc.cachefunc import cached_method
from sage.rings.finite_rings.finite_field_base import FiniteField
from sage.rings.finite_rings.finite_field_constructor import FiniteField as GF
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.polynomial_ring import polygen
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.schemes.curves.projective_curve import Hasse_bounds
from sage.structure.element import Element

from . import ell_point
from .constructor import EllipticCurve
from .ell_field import EllipticCurve_field


class EllipticCurve_finite_field(EllipticCurve_field):
    r"""
    Elliptic curve over a finite field.

    EXAMPLES::

        sage: EllipticCurve(GF(101),[2,3])
        Elliptic Curve defined by y^2  = x^3 + 2*x + 3 over Finite Field of size 101

        sage: # needs sage.rings.finite_rings
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
        Category of abelian varieties over Ring of integers modulo 101

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

        - ``*args``, ``**kwds`` -- all other options are passed
          to the circle graphing primitive

        EXAMPLES::

            sage: E = EllipticCurve(FiniteField(17), [0,1])
            sage: P = plot(E, rgbcolor=(0,0,1))                                         # needs sage.plot
        """
        R = self.base_ring()
        if not R.is_prime_field():
            raise NotImplementedError

        from sage.plot.point import points

        return points([P[0:2] for P in self.points() if not P.is_zero()], *args, **kwds)

    def _points_via_group_structure(self):
        r"""
        Return a list of all the points on the curve using the group structure.

        For cyclic groups with generator `G` of order `n`, the set of points
        is computed from `[x]G` for `x \in [0, n - 1]` requiring `n - 2`
        additions on the curve.

        For non-cyclic groups, first two cyclic subgroups `H_i` are computed as
        above from `[x]G_i` for `x \in [0, n_i]` requiring `n_1 + n_2 - 4`
        additions. The set of all points is returned as the cartesian product
        of these two cyclic groups.

        When the group is trivial, only the point at infinity is returned.

        EXAMPLES::

            sage: S = EllipticCurve(GF(97),[2,3])._points_via_group_structure()
            sage: len(S)
            100

        See :issue:`4687`, where the following example did not work::

            sage: E = EllipticCurve(GF(2),[0, 0, 1, 1, 1])
            sage: E.points()
            [(0 : 1 : 0)]

        ::

            sage: E = EllipticCurve(GF(2),[0, 0, 1, 0, 1])
            sage: E.points()
            [(0 : 1 : 0), (1 : 0 : 1), (1 : 1 : 1)]

        ::

            sage: # needs sage.rings.finite_rings
            sage: E = EllipticCurve(GF(4,'a'),[0, 0, 1, 0, 1])
            sage: E.points()
            [(0 : 1 : 0), (0 : a : 1), (0 : a + 1 : 1), (1 : 0 : 1), (1 : 1 : 1), (a : 0 : 1), (a : 1 : 1), (a + 1 : 0 : 1), (a + 1 : 1 : 1)]
        """
        # TODO, eliminate when polynomial calling is fast
        # 11-03-2024 - G. Pope : it is not clear to me what the above TODO references

        # Compute the generators of the abelian group of the curve
        G = self.abelian_group()
        gens = [x.element() for x in G.gens()]

        # Zero element of the group
        zero = self(0)

        # Trivial group, return the identity element only
        if len(gens) == 0:
            return [zero]

        def __multiples(G):
            """
            Compute the list of points [i]G for i in [0, G.order())
            """
            H = [zero, G]
            P = G
            for _ in range(2, G.order()):
                P += G
                H.append(P)
            return H

        # Collect all multiples of the generator
        H1 = __multiples(gens[0])

        # Cyclic case, we now have all points
        if len(gens) == 1:
            return H1

        # Non-cyclic case we generate the second set of points and compute
        # the entire set of points from the Cartesian product
        H2 = __multiples(gens[1])
        return [P + Q for P in H1 for Q in H2]

    def points(self):
        r"""
        Return all rational points on this elliptic curve. The list of points is cached
        so subsequent calls are free.

        EXAMPLES::

            sage: p = 5
            sage: F = GF(p)
            sage: E = EllipticCurve(F, [1, 3])
            sage: len(E.points())
            4
            sage: E.order()
            4
            sage: E.points()
            [(0 : 1 : 0), (1 : 0 : 1), (4 : 1 : 1), (4 : 4 : 1)]

        ::

            sage: K = GF((p, 2), 'a')
            sage: E = E.change_ring(K)
            sage: len(E.points())
            32
            sage: E.order()
            32
            sage: w = E.points(); w
            [(0 : 1 : 0), (0 : 2*a + 4 : 1), (0 : 3*a + 1 : 1), (1 : 0 : 1), (2 : 2*a + 4 : 1), (2 : 3*a + 1 : 1), (3 : 2*a + 4 : 1), (3 : 3*a + 1 : 1), (4 : 1 : 1), (4 : 4 : 1), (a : 1 : 1), (a : 4 : 1), (a + 2 : a + 1 : 1), (a + 2 : 4*a + 4 : 1), (a + 3 : a : 1), (a + 3 : 4*a : 1), (a + 4 : 0 : 1), (2*a : 2*a : 1), (2*a : 3*a : 1), (2*a + 4 : a + 1 : 1), (2*a + 4 : 4*a + 4 : 1), (3*a + 1 : a + 3 : 1), (3*a + 1 : 4*a + 2 : 1), (3*a + 2 : 2*a + 3 : 1), (3*a + 2 : 3*a + 2 : 1), (4*a : 0 : 1), (4*a + 1 : 1 : 1), (4*a + 1 : 4 : 1), (4*a + 3 : a + 3 : 1), (4*a + 3 : 4*a + 2 : 1), (4*a + 4 : a + 4 : 1), (4*a + 4 : 4*a + 1 : 1)]

        Note that the returned list is an immutable sorted Sequence::

            sage: w[0] = 9
            Traceback (most recent call last):
            ...
            ValueError: object is immutable; please change a copy instead.
        """
        if hasattr(self, "__points"):
            return self.__points

        from sage.structure.sequence import Sequence
        v = self._points_via_group_structure()
        v.sort()
        self.__points = Sequence(v, immutable=True)
        return self.__points

    rational_points = points

    def count_points(self, n=1):
        """
        Return the cardinality of this elliptic curve over the base field or extensions.

        INPUT:

        - ``n`` -- positive integer

        OUTPUT: if `n=1`, returns the cardinality of the curve over its base field

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

            sage: # needs sage.rings.finite_rings
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
        and compute the possible y-coordinates. Return the i-th
        possible y-coordinate, where i is randomly chosen to be 0 or 1.
        If the i-th y-coordinate does not exist (either there is no
        point with the given x-coordinate or we hit a 2-torsion point
        with i == 1), try again.

        This gives a uniform distribution because you can imagine
        `2q + 1` buckets, one for the point at infinity and 2 for each
        element of the field (representing the x-coordinates). This
        gives a 1-to-1 map of elliptic curve points into buckets. At
        every iteration, we simply choose a random bucket until we find
        a bucket containing a point.

        AUTHORS:

        - Jeroen Demeyer (2014-09-09): choose points uniformly random,
          see :issue:`16951`.

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

            sage: # needs sage.rings.finite_rings
            sage: k.<a> = GF(7^5)
            sage: E = EllipticCurve(k,[2,4])
            sage: P = E.random_element(); P  # random
            (5*a^4 + 3*a^3 + 2*a^2 + a + 4 : 2*a^4 + 3*a^3 + 4*a^2 + a + 5 : 1)
            sage: type(P)
            <class 'sage.schemes.elliptic_curves.ell_point.EllipticCurvePoint_finite_field'>
            sage: P in E
            True

        ::

            sage: # needs sage.rings.finite_rings
            sage: k.<a> = GF(2^5)
            sage: E = EllipticCurve(k,[a^2,a,1,a+1,1])
            sage: P = E.random_element(); P  # random
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

        See :issue:`8311`::

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

            sage: # needs sage.rings.finite_rings
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
            sage: E = EllipticCurve(GF(11^5,'a'),[2,5])                                 # needs sage.rings.finite_rings
            sage: E.trace_of_frobenius()                                                # needs sage.rings.finite_rings
            802

        The following shows that the issue from :issue:`2849` is fixed::

            sage: E = EllipticCurve(GF(3^5,'a'),[-1,-1])                                # needs sage.rings.finite_rings
            sage: E.trace_of_frobenius()                                                # needs sage.rings.finite_rings
            -27
        """
        return 1 + self.base_field().order() - self.cardinality()

    def cardinality(self, algorithm=None, extension_degree=1):
        r"""
        Return the number of points on this elliptic curve.

        INPUT:

        - ``algorithm`` -- (optional) string:

          - ``'pari'`` -- use the PARI C-library function ``ellcard``

          - ``'bsgs'`` -- use the baby-step giant-step method as
             implemented in Sage, with the Cremona-Sutherland version
             of Mestre's trick

          - ``'exhaustive'`` -- naive point counting

          - ``'subfield'`` -- reduce to a smaller field, provided that
            the j-invariant lies in a subfield

          - ``'all'`` -- compute cardinality with both ``'pari'`` and
            ``'bsgs'``; return result if they agree or raise a
            :exc:`AssertionError` if they do not

        - ``extension_degree`` -- integer `d` (default: 1); if the
          base field is `\GF{q}`, return the cardinality of ``self``
          over the extension `\GF{q^d}` of degree `d`

        OUTPUT:

        The order of the group of rational points of ``self`` over its
        base field, or over an extension field of degree `d` as above.
        The result is cached.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: EllipticCurve(GF(4, 'a'), [1,2,3,4,5]).cardinality()
            8
            sage: k.<a> = GF(3^3)
            sage: l = [a^2 + 1, 2*a^2 + 2*a + 1, a^2 + a + 1, 2, 2*a]
            sage: EllipticCurve(k,l).cardinality()
            29

        ::

            sage: # needs sage.rings.finite_rings
            sage: l = [1, 1, 0, 2, 0]
            sage: EllipticCurve(k, l).cardinality()
            38

        An even bigger extension (which we check against Magma)::

            sage: # needs sage.rings.finite_rings
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

            sage: # needs sage.rings.finite_rings
            sage: E = EllipticCurve(GF(3^100, 'a'), [1,2,3,4,5])
            sage: E.cardinality() is E.cardinality()
            True

        The following is very fast since the curve is actually defined
        over the prime field::

            sage: # needs sage.rings.finite_rings
            sage: k.<a> = GF(11^100)
            sage: E1 = EllipticCurve(k, [3,3])
            sage: N1 = E1.cardinality(algorithm='subfield'); N1
            137806123398222701841183371720896367762643312000384671846835266941791510341065565176497846502742959856128
            sage: E1.cardinality_pari() == N1
            True
            sage: E2 = E1.quadratic_twist()
            sage: N2 = E2.cardinality(algorithm='subfield'); N2
            137806123398222701841183371720896367762643312000384656816094284101308193849980588362304472492174093035876
            sage: E2.cardinality_pari() == N2
            True
            sage: N1 + N2 == 2*(k.cardinality() + 1)
            True

        We can count points over curves defined as a reduction::

            sage: # needs sage.rings.number_field
            sage: x = polygen(QQ)
            sage: K.<w> = NumberField(x^2 + x + 1)
            sage: EK = EllipticCurve(K, [0, 0, w, 2, 1])
            sage: E = EK.base_extend(K.residue_field(2))
            sage: E
            Elliptic Curve defined by y^2 + wbar*y = x^3 + 1
             over Residue field in wbar of Fractional ideal (2)
            sage: E.cardinality()
            7
            sage: E = EK.base_extend(K.residue_field(w - 1))
            sage: E.abelian_group()
            Trivial group embedded in Abelian group of points on Elliptic Curve defined
             by y^2 + y = x^3 + 2*x + 1 over Residue field of Fractional ideal (w - 1)

        ::

            sage: R.<x> = GF(17)[]
            sage: pol = R.irreducible_element(5)
            sage: k.<a> = R.residue_field(pol)
            sage: E = EllipticCurve(R, [1, x]).base_extend(k)
            sage: E
            Elliptic Curve defined by y^2 = x^3 + x + a
             over Residue field in a of Principal ideal (x^5 + x + 14)
              of Univariate Polynomial Ring in x over Finite Field of size 17
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

        Check that a bug noted at :issue:`15667` is fixed::

            sage: # needs sage.rings.finite_rings
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

    @cached_method
    def multiplication_by_p_isogeny(self):
        r"""
        Return the multiplication-by-`p` isogeny.

        EXAMPLES::

            sage: p = 23
            sage: K.<a> = GF(p^3)
            sage: E = EllipticCurve(j=K.random_element())
            sage: phi = E.multiplication_by_p_isogeny()
            sage: assert phi.degree() == p**2
            sage: P = E.random_element()
            sage: assert phi(P) == P * p
        """
        frob = self.frobenius_isogeny()
        return frob.dual() * frob

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

            sage: # needs sage.rings.finite_rings
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
            Order of conductor 2 generated by phi
             in Number Field in phi with defining polynomial x^2 - 4*x + 11

        For some supersingular curves, Frobenius is in Z and the Frobenius
        order is Z::

            sage: # needs sage.rings.finite_rings
            sage: E = EllipticCurve(GF(25,'a'),[0,0,0,0,1])
            sage: R = E.frobenius_order()
            sage: R
            Order generated by []
             in Number Field in phi with defining polynomial x + 5
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

            sage: # needs sage.rings.finite_rings
            sage: E = EllipticCurve(GF(25,'a'),[0,0,0,0,1])
            sage: E.frobenius()
            -5
        """
        R = self.frobenius_order()
        if R.degree() == 1:
            return self.frobenius_polynomial().roots(multiplicities=False)[0]
        else:
            return R.gen(1)

    def frobenius_endomorphism(self):
        r"""
        Return the `q`-power Frobenius endomorphism of this elliptic
        curve, where `q` is the cardinality of the (finite) base field.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: F.<t> = GF(11^4)
            sage: E = EllipticCurve([t,t])
            sage: E.frobenius_endomorphism()
            Frobenius endomorphism of degree 14641 = 11^4:
              From: Elliptic Curve defined by y^2 = x^3 + t*x + t over Finite Field in t of size 11^4
              To:   Elliptic Curve defined by y^2 = x^3 + t*x + t over Finite Field in t of size 11^4
            sage: E.frobenius_endomorphism() == E.frobenius_isogeny(4)
            True

        .. SEEALSO::

            :meth:`~sage.schemes.elliptic_curves.ell_generic.EllipticCurve_generic.frobenius_isogeny`
        """
        return self.frobenius_isogeny(self.base_field().degree())

    def frobenius_discriminant(self):
        r"""
        Return the discriminant of the ring `\ZZ[\pi_E]` where `\pi_E` is the Frobenius endomorphism.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: F.<t> = GF(11^4)
            sage: E = EllipticCurve([t,t])
            sage: E.frobenius_discriminant()
            -57339
        """
        return self.frobenius_polynomial().discriminant()

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

        Since :issue:`16931`, this now works over finite fields which
        are not prime fields::

            sage: # needs sage.rings.finite_rings
            sage: k.<a> = GF(7^3)
            sage: E = EllipticCurve_from_j(a)
            sage: E.cardinality_pari()
            318
            sage: K.<a> = GF(3^20)
            sage: E = EllipticCurve(K,[1,0,0,1,a])
            sage: E.cardinality_pari()
            3486794310

        TESTS::

            sage: E.cardinality_pari().parent()                                         # needs sage.rings.finite_rings
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

        OUTPUT: a tuple of points on the curve

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
            Additive abelian group isomorphic to Z/22 + Z/2
             embedded in Abelian group of points on Elliptic Curve
             defined by y^2 = x^3 + 2*x + 5 over Finite Field of size 41
            sage: ab_gens = E.abelian_group().gens()
            sage: ab_gens == E.gens()
            True
            sage: E.gens()[0].order()
            22
            sage: E.gens()[1].order()
            2

        Examples with 1 and 0 generators::

            sage: # needs sage.rings.finite_rings
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

            sage: # needs sage.rings.finite_rings
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
        yield from self.points()

    def __getitem__(self, n):
        """
        Return the n-th point in ``self``'s ``__points`` list.

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

        .. SEEALSO::

            :meth:`AdditiveAbelianGroupWrapper.from_generators()<sage.groups.additive_abelian.additive_abelian_wrapper.AdditiveAbelianGroupWrapper.from_generators>`

        EXAMPLES::

            sage: E = EllipticCurve(GF(11),[2,5])
            sage: E.abelian_group()
            Additive abelian group isomorphic to Z/10 embedded in
             Abelian group of points on Elliptic Curve defined by y^2 = x^3 + 2*x + 5
              over Finite Field of size 11

        ::

            sage: E = EllipticCurve(GF(41),[2,5])
            sage: E.abelian_group()
            Additive abelian group isomorphic to Z/22 + Z/2 ...

        ::

            sage: # needs sage.rings.finite_rings
            sage: F.<a> = GF(3^6,'a')
            sage: E = EllipticCurve([a^4 + a^3 + 2*a^2 + 2*a, 2*a^5 + 2*a^3 + 2*a^2 + 1])
            sage: E.abelian_group()
            Additive abelian group isomorphic to Z/26 + Z/26 ...

        ::

            sage: # needs sage.rings.finite_rings
            sage: F.<a> = GF(101^3,'a')
            sage: E = EllipticCurve([2*a^2 + 48*a + 27, 89*a^2 + 76*a + 24])
            sage: E.abelian_group()
            Additive abelian group isomorphic to Z/1031352 ...

        The group can be trivial::

            sage: E = EllipticCurve(GF(2), [0,0,1,1,1])
            sage: E.abelian_group()
            Trivial group embedded in Abelian group of points on
             Elliptic Curve defined by y^2 + y = x^3 + x + 1 over Finite Field of size 2

        Of course, there are plenty of points if we extend the field::

            sage: E.cardinality(extension_degree=100)
            1267650600228231653296516890625

        This tests the patch for :issue:`3111`, using 10 primes randomly
        selected::

            sage: E = EllipticCurve('389a')
            sage: for p in [5927, 2297, 1571, 1709, 3851, 127, 3253, 5783, 3499, 4817]:
            ....:     G = E.change_ring(GF(p)).abelian_group()
            sage: for p in prime_range(10000):  # long time (19s on sage.math, 2011)
            ....:     if p != 389:
            ....:         G = E.change_ring(GF(p)).abelian_group()

        This tests that the bug reported in :issue:`3926` has been fixed::

            sage: # needs sage.rings.number_field
            sage: K.<i> = QuadraticField(-1)
            sage: OK = K.ring_of_integers()
            sage: P = K.factor(10007)[0][0]
            sage: OKmodP = OK.residue_field(P)
            sage: E = EllipticCurve([0, 0, 0, i, i + 3])
            sage: Emod = E.change_ring(OKmodP); Emod
            Elliptic Curve defined by y^2 = x^3 + ibar*x + (ibar+3)
             over Residue field in ibar of Fractional ideal (10007)
            sage: Emod.abelian_group() #random generators
            (Multiplicative Abelian group isomorphic to C50067594 x C2,
             ((3152*ibar + 7679 : 7330*ibar + 7913 : 1), (8466*ibar + 1770 : 0 : 1)))
        """

        gens = self.gens()
        assert len(gens) <= 2

        if len(gens) == 2:

            P, Q = gens
            n = self.cardinality()              # cached
            n1 = P.order()                      # cached
            n2 = n//n1
            assert not n1 * Q                   # PARI should guarantee this

            k = n1.prime_to_m_part(n2)
            Q *= k                              # don't need; kill that part
            nQ = n2 * generic.order_from_multiple(n2*Q, n1//k//n2)

            S = n//nQ * P
            T = n2 * Q
            S.set_order(nQ//n2, check=False)    # for .log()
            x = T.log(S)
            Q -= x * n1//nQ * P

            assert not n2 * Q                   # by construction
            Q.set_order(n2, check=False)

            gens = P, Q

        orders = [T.order() for T in gens]      # cached

        self.gens.set_cache(gens)
        return AdditiveAbelianGroupWrapper(self.point_homset(), gens, orders)

    def torsion_basis(self, n):
        r"""
        Return a basis of the `n`-torsion subgroup of this elliptic curve,
        assuming it is fully rational.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: E = EllipticCurve(GF(62207^2), [1,0])
            sage: E.abelian_group()
            Additive abelian group isomorphic to Z/62208 + Z/62208 embedded in
             Abelian group of points on Elliptic Curve defined by y^2 = x^3 + x
              over Finite Field in z2 of size 62207^2
            sage: PA,QA = E.torsion_basis(2^8)
            sage: PA.weil_pairing(QA, 2^8).multiplicative_order()
            256
            sage: PB,QB = E.torsion_basis(3^5)
            sage: PB.weil_pairing(QB, 3^5).multiplicative_order()
            243

        ::

            sage: E = EllipticCurve(GF(101), [4,4])
            sage: E.torsion_basis(23)
            Traceback (most recent call last):
            ...
            ValueError: curve does not have full rational 23-torsion
            sage: F = E.division_field(23); F
            Finite Field in t of size 101^11
            sage: EE = E.change_ring(F)
            sage: P, Q = EE.torsion_basis(23)
            sage: P  # random
            (89*z11^10 + 51*z11^9 + 96*z11^8 + 8*z11^7 + 67*z11^6
             + 31*z11^5 + 55*z11^4 + 59*z11^3 + 28*z11^2 + 8*z11 + 88
             : 40*z11^10 + 33*z11^9 + 80*z11^8 + 87*z11^7 + 97*z11^6
             + 69*z11^5 + 56*z11^4 + 17*z11^3 + 26*z11^2 + 69*z11 + 11
             : 1)
            sage: Q  # random
            (25*z11^10 + 61*z11^9 + 49*z11^8 + 17*z11^7 + 80*z11^6
             + 20*z11^5 + 49*z11^4 + 52*z11^3 + 61*z11^2 + 27*z11 + 61
             : 60*z11^10 + 91*z11^9 + 89*z11^8 + 7*z11^7 + 63*z11^6
             + 55*z11^5 + 23*z11^4 + 17*z11^3 + 90*z11^2 + 91*z11 + 68
             : 1)

        .. SEEALSO::

            Use :meth:`~sage.schemes.elliptic_curves.ell_field.EllipticCurve_field.division_field`
            to determine a field extension containing the full `\ell`-torsion subgroup.

        ALGORITHM:

        This method currently uses :meth:`abelian_group` and
        :meth:`AdditiveAbelianGroupWrapper.torsion_subgroup`.
        """
        # TODO: In many cases this is not the fastest algorithm.
        # Alternatives include factoring division polynomials and
        # random sampling (like PARI's ellgroup, but with a milder
        # termination condition). We should implement these too
        # and figure out when to use which.
        T = self.abelian_group().torsion_subgroup(n)
        if T.invariants() != (n, n):
            raise ValueError(f'curve does not have full rational {n}-torsion')
        return tuple(P.element() for P in T.gens())

    def is_isogenous(self, other, field=None, proof=True):
        """
        Return whether or not ``self`` is isogenous to ``other``.

        INPUT:

        - ``other`` -- another elliptic curve

        - ``field`` -- (default: ``None``) a field containing the base
          fields of the two elliptic curves into which the two curves
          may be extended to test if they are isogenous over this
          field. By default is_isogenous will not try to find this
          field unless one of the curves can be extended into the base
          field of the ``other``, in which case it will test over the
          larger base field.

        - ``proof`` -- boolean (default: ``True``); this parameter is here only
          to be consistent with versions for other types of elliptic curves

        OUTPUT:

        boolean; ``True`` if there is an isogeny from curve ``self`` to
        curve ``other`` defined over ``field``

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: E1 = EllipticCurve(GF(11^2,'a'),[2,7]); E1
            Elliptic Curve defined by y^2 = x^3 + 2*x + 7 over Finite Field in a of size 11^2
            sage: E1.is_isogenous(5)
            Traceback (most recent call last):
            ...
            ValueError: Second argument is not an Elliptic Curve.
            sage: E1.is_isogenous(E1)
            True

            sage: # needs sage.rings.finite_rings
            sage: E2 = EllipticCurve(GF(7^3,'b'),[3,1]); E2
            Elliptic Curve defined by y^2 = x^3 + 3*x + 1 over Finite Field in b of size 7^3
            sage: E1.is_isogenous(E2)
            Traceback (most recent call last):
            ...
            ValueError: The base fields must have the same characteristic.

            sage: # needs sage.rings.finite_rings
            sage: E3 = EllipticCurve(GF(11^2,'c'),[4,3]); E3
            Elliptic Curve defined by y^2 = x^3 + 4*x + 3 over Finite Field in c of size 11^2
            sage: E1.is_isogenous(E3)
            False

            sage: # needs sage.rings.finite_rings
            sage: E4 = EllipticCurve(GF(11^6,'d'),[6,5]); E4
            Elliptic Curve defined by y^2 = x^3 + 6*x + 5 over Finite Field in d of size 11^6
            sage: E1.is_isogenous(E4)
            True

            sage: # needs sage.rings.finite_rings
            sage: E5 = EllipticCurve(GF(11^7,'e'),[4,2]); E5
            Elliptic Curve defined by y^2 = x^3 + 4*x + 2 over Finite Field in e of size 11^7
            sage: E1.is_isogenous(E5)
            Traceback (most recent call last):
            ...
            ValueError: Curves have different base fields: use the field parameter.

        When the field is given::

            sage: # needs sage.rings.finite_rings
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
        from .ell_generic import EllipticCurve_generic
        if not isinstance(other, EllipticCurve_generic):
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
        Return ``True`` if this elliptic curve is supersingular, else ``False``.

        INPUT:

        - ``proof``-- boolean (default: ``True``); if ``True``, returns a
          proved result.  If ``False``, then a return value of ``False`` is
          certain but a return value of ``True`` may be based on a
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
        Return ``True`` if this elliptic curve is ordinary, else ``False``.

        INPUT:

        - ``proof``-- boolean (default: ``True``); if ``True``, returns a
          proved result.  If ``False``, then a return value of ``True`` is
          certain but a return value of ``False`` may be based on a
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

    def set_order(self, value, *, check=True, num_checks=8):
        r"""
        Set the value of ``self._order`` to ``value``.

        Use this when you know a priori the order of the curve to
        avoid a potentially expensive order calculation.

        INPUT:

        - ``value`` -- integer in the Hasse-Weil range for this curve

        - ``check``-- boolean (default: ``True``); whether or
          not to run sanity checks on the input

        - ``num_checks``-- integer (default: 8); if ``check`` is
          ``True``, the number of times to check whether ``value``
          times a random point on this curve equals the identity

        OUTPUT: none

        EXAMPLES:

        This example illustrates basic usage::

            sage: E = EllipticCurve(GF(7), [0, 1]) # This curve has order 12
            sage: E.set_order(12)
            sage: E.order()
            12
            sage: E.order() * E.random_point()
            (0 : 1 : 0)

        We now give a more interesting case, the NIST-P521 curve. Its
        order is too big to calculate with Sage, and takes a long time
        using other packages, so it is very useful here::

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

            sage: E = EllipticCurve(GF(7), [0, 1]) # This curve has order 12
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
        ``num_checks``, the factorization of the actual order, and the
        actual group structure::

            sage: E = EllipticCurve(GF(1009), [0, 1]) # This curve has order 948
            sage: E.set_order(947)
            Traceback (most recent call last):
            ...
            ValueError: Value 947 illegal (multiple of random point not the identity)

        For curves over small finite fields, the order is cheap to compute, so it is computed
        directly and compared::

            sage: E = EllipticCurve(GF(7), [0, 1]) # This curve has order 12
            sage: E.set_order(11)
            Traceback (most recent call last):
            ...
            ValueError: Value 11 illegal (correct order is 12)

        TESTS:

        The previous version's random tests are not strong enough. In particular, the following used
        to work::

            sage: E = EllipticCurve(GF(2), [0, 0, 1, 1, 1]) # This curve has order 1
            sage: E.set_order(3)
            Traceback (most recent call last):
            ...
            ValueError: Value 3 illegal (correct order is 1)

        ::

            sage: E = EllipticCurve(GF(7), [0, 1]) # This curve has order 12
            sage: E.set_order(4, num_checks=0)
            Traceback (most recent call last):
            ...
            ValueError: Value 4 illegal (correct order is 12)
            sage: E.order()
            12

        .. TODO:: Add provable correctness check by computing the abelian group structure and
            comparing.

        AUTHORS:

         - Mariah Lenox (2011-02-16): Initial implementation

         - Gareth Ma (2024-01-21): Fix bug for small curves
        """
        value = Integer(value)

        if check:
            # Is value in the Hasse range?
            q = self.base_field().order()
            a,b = Hasse_bounds(q,1)
            if not a <= value <= b:
                raise ValueError(f"Value {value} illegal (not an integer in the Hasse range)")

            # For really small values, the random tests are too weak to detect wrong orders
            # So we go with computing directly instead.
            if q <= 100:
                if self.order() != value:
                    raise ValueError(f"Value {value} illegal (correct order is {self.order()})")

            # Is value*random == identity?
            for _ in range(num_checks):
                G = self.random_point()
                if value * G != self(0):
                    raise ValueError(f"Value {value} illegal (multiple of random point not the identity)")

        # TODO: It might help some of PARI's algorithms if we
        # could copy this over to the .pari_curve() as well.
        # At the time of writing, this appears to be tricky to
        # do in a non-hacky way because cypari2 doesn't expose
        # "member functions" of PARI objects.

        self._order = value

    def _fetch_cached_order(self, other):
        r"""
        This method copies the ``_order`` member from ``other`` to
        ``self``. Both curves must have the same finite base field.

        This is used in
        :class:`~sage.schemes.elliptic_curves.hom.EllipticCurveHom`
        to keep track of an already computed curve order: According
        to Tate's theorem [Tate1966b]_, isogenous elliptic curves
        over a finite field have the same number of rational points.

        EXAMPLES::

            sage: E1 = EllipticCurve(GF(2^127-1), [1,2,3,4,5])
            sage: E1.set_order(170141183460469231746191640949390434666)
            sage: E2 = EllipticCurve(GF(2^127-1), [115649500210559831225094148253060920818, 36348294106991415644658737184600079491])
            sage: E2._fetch_cached_order(E1)
            sage: E2._order
            170141183460469231746191640949390434666

        TESTS::

            sage: E3 = EllipticCurve(GF(17), [1,2,3,4,5])
            sage: hasattr(E3, '_order')
            False
            sage: E3._fetch_cached_order(E1)
            Traceback (most recent call last):
            ...
            ValueError: curves have distinct base fields
        """
        if hasattr(self, '_order') or not hasattr(other, '_order'):
            return
        F = self.base_field()
        if F != other.base_field():
            raise ValueError('curves have distinct base fields')
        n = getattr(other, '_order', None)
        if n is not None:
            self._order = n

    def height_above_floor(self, ell, e):
        r"""
        Return the height of the `j`-invariant of this ordinary elliptic curve on its `\ell`-volcano.

        INPUT:

        - ``ell`` -- a prime number
        - ``e`` -- nonnegative integer, the `\ell`-adic valuation of
          the conductor the Frobenius order


        .. NOTE::

            For an ordinary `E/\GF{q}`, and a prime `\ell`, the height
            `e` of the `\ell`-volcano containing `j(E)` is the `\ell`-adic
            valuation of the conductor of the order generated by the
            Frobenius `\pi_E`; the height of `j(E)` on its
            ell-volcano is the `\ell`-adic valuation of the conductor
            of the order `\text{End}(E)`.

        ALGORITHM:

            See [RouSuthZur2022]_.

        EXAMPLES::

            sage: F = GF(312401)
            sage: E = EllipticCurve(F,(0, 0, 0, 309381, 93465))
            sage: D = E.frobenius_discriminant(); D
            -687104
            sage: D.factor()
            -1 * 2^10 * 11 * 61
            sage: E.height_above_floor(2,8)
            5
        """
        if self.is_supersingular():
            raise ValueError("{} is not ordinary".format(self))
        if e == 0:
            return 0
        j = self.j_invariant()
        if j in [0, 1728]:
            return e
        F = j.parent()
        x = polygen(F)
        from sage.rings.polynomial.polynomial_ring import polygens
        from sage.schemes.elliptic_curves.mod_poly import classical_modular_polynomial
        X, Y = polygens(F, "X, Y", 2)
        phi = classical_modular_polynomial(ell)(X, Y)
        j1 = phi([x,j]).roots(multiplicities=False)
        nj1 = len(j1)
        on_floor = self.two_torsion_rank() < 2 if ell == 2 else nj1 <= ell
        if on_floor:
            return 0
        if e == 1 or nj1 != ell+1:  # double roots can only happen at the surface
            return e
        if nj1 < 3:
            return 0
        j0 = [j,j,j]
        h = 1
        while True:
            for i in range(3):
                r = (phi([x,j1[i]])//(x-j0[i])).roots(multiplicities=False)
                if not r:
                    return h
                j0[i] = j1[i]
                j1[i] = r[0]
            h += 1

    def endomorphism_discriminant_from_class_number(self, h):
        r"""
        Return the endomorphism order discriminant of this ordinary elliptic curve, given its class number ``h``.

        INPUT:

        - ``h`` -- positive integer

        OUTPUT:

        integer; the discriminant of the endomorphism ring `\text{End}(E)`, if
        this has class number ``h``.  If `\text{End}(E)` does not have class
        number ``h``, a :exc:`ValueError` is raised.

        ALGORITHM:

        Compute the trace of Frobenius and hence the discriminant
        `D_0` and class number `h_0` of the maximal order containing
        the endomorphism order.  From the given value of `h`, which
        must be a multiple of `h_0`, compute the possible conductors,
        using :meth:`height_above_floor` for each prime `\ell`
        dividing the quotient `h/h_0`.  If exactly one conductor `f`
        remains, return `f^2D_0`, otherwise raise a :exc:`ValueError`;
        this can onlyhappen when the input value of `h` was incorrect.

        .. NOTE::

            Adapted from [RouSuthZur2022]_.  The application for which
            one knows the class number in advance is in the
            recognition of Hilbert Class Polynomials: see
            :func:`sage.schemes.elliptic_curves.cm.is_HCP`.

        EXAMPLES::

            sage: F = GF(312401)
            sage: E = EllipticCurve(F,(0, 0, 0, 309381, 93465))
            sage: E.endomorphism_discriminant_from_class_number(30)
            -671

        We check that this is the correct discriminant, and the input value of `h` was correct::

            sage: H = hilbert_class_polynomial(-671)
            sage: H(E.j_invariant()) == 0 and H.degree()==30
            True
        """
        F = self.base_field()
        if not F.is_finite():
            raise ValueError("Base field {} must be finite".format(F))
        if self.is_supersingular():
            raise ValueError("Elliptic curve ({}) must be ordinary".format(self))
        D1 = self.frobenius_discriminant()
        D0 = D1.squarefree_part()
        if D0 % 4 != 1:
            D0 *= 4
        v = ZZ(D1//D0).isqrt()
        h0 = D0.class_number()
        if h % h0:
            raise ValueError("Incorrect class number {}".format(h))
        from sage.schemes.elliptic_curves.cm import OrderClassNumber
        cs = [v//f for f in v.divisors() if OrderClassNumber(D0,h0,f) == h] # cofactors c=v/f compatible with h(f**2D0)=h
        if not cs:
            raise ValueError("Incorrect class number {}".format(h))
        if len(cs) == 1:
            return (v//cs[0])**2 * D0

        L = sorted(set(sum([c.prime_factors() for c in cs], [])))
        for ell in L:
            e = self.height_above_floor(ell,v.valuation(ell))
            cs = [c for c in cs if c.valuation(ell) == e]
            if not cs:
                raise ValueError("Incorrect class number {}".format(h))
            if len(cs) == 1:
                return (v//cs[0])**2 * D0
        raise ValueError("Incorrect class number {}".format(h))

    def twists(self):
        r"""
        Return a list of `k`-isomorphism representatives of all
        twists of this elliptic curve, where `k` is the base field.

        The input curve appears as the first entry of the result.

        .. NOTE::

            A *twist* of `E/k` is an elliptic curve `E'` defined over
            `k` that is isomorphic to `E` over the algebraic closure
            `\bar k`.

            Most elliptic curves over a finite field only admit a
            single nontrivial twist (the quadratic twist); the only
            exceptions are curves with `j`-invariant `0` or `1728`.

            In all cases the sum over all the twists `E'` of `1/|Aut(E')|` is 1.

        .. SEEALSO::

            - :meth:`~sage.schemes.elliptic_curves.ell_field.EllipticCurve_field.quadratic_twist`
            - :meth:`~sage.schemes.elliptic_curves.ell_field.EllipticCurve_field.quartic_twist`
            - :meth:`~sage.schemes.elliptic_curves.ell_field.EllipticCurve_field.sextic_twist`

        EXAMPLES::

            sage: E = EllipticCurve(GF(97), [1,1])
            sage: E.j_invariant()
            54
            sage: E.twists()
            [Elliptic Curve defined by y^2 = x^3 + x + 1 over Finite Field of size 97,
             Elliptic Curve defined by y^2 = x^3 + ... over Finite Field of size 97]

        ::

            sage: E = EllipticCurve(GF(97), [1,0])
            sage: E.j_invariant()
            79
            sage: E.twists()
            [Elliptic Curve defined by y^2 = x^3 + x over Finite Field of size 97,
             Elliptic Curve defined by y^2 = x^3 + ... over Finite Field of size 97,
             Elliptic Curve defined by y^2 = x^3 + ... over Finite Field of size 97,
             Elliptic Curve defined by y^2 = x^3 + ... over Finite Field of size 97]

        ::

            sage: E = EllipticCurve(GF(97), [0,1])
            sage: E.j_invariant()
            0
            sage: E.twists()
            [Elliptic Curve defined by y^2 = x^3 + 1 over Finite Field of size 97,
             Elliptic Curve defined by y^2 = x^3 + ... over Finite Field of size 97,
             Elliptic Curve defined by y^2 = x^3 + ... over Finite Field of size 97,
             Elliptic Curve defined by y^2 = x^3 + ... over Finite Field of size 97,
             Elliptic Curve defined by y^2 = x^3 + ... over Finite Field of size 97,
             Elliptic Curve defined by y^2 = x^3 + ... over Finite Field of size 97]

        This can be useful to quickly compute a list of all elliptic curves
        over a finite field `k` up to `k`-isomorphism::

            sage: Es = [E for j in GF(13) for E in EllipticCurve(j=j).twists()]
            sage: len(Es)
            32
            sage: Es
            [Elliptic Curve defined by y^2 = x^3 + 1 over Finite Field of size 13,
             ...
             Elliptic Curve defined by y^2 = x^3 + ... over Finite Field of size 13]

        In characteristic 3, the number of twists is 2 except for
        `j=0=1728`, when there are either 4 or 6 depending on whether the
        field has odd or even degree over `\GF{3}`::

            sage: # needs sage.rings.finite_rings
            sage: K = GF(3**5)
            sage: [E.ainvs() for E in EllipticCurve(j=K(1)).twists()]
            [(0, 1, 0, 0, 2), (0, z5, 0, 0, 2*z5^3)]

            sage: # needs sage.rings.finite_rings
            sage: K = GF(3**5)
            sage: [E.ainvs() for E in EllipticCurve(j=K(0)).twists()] # random
            [(0, 0, 0, 1, 0),
             (0, 0, 0, 2, 0),
             (0, 0, 0, 2, z5^4 + z5^3 + z5^2),
             (0, 0, 0, 2, 2*z5^4 + 2*z5^3 + 2*z5^2)]

            sage: # needs sage.rings.finite_rings
            sage: K = GF(3**4)
            sage: [E.ainvs() for E in EllipticCurve(j=K(1)).twists()]
            [(0, 1, 0, 0, 2), (0, z4, 0, 0, 2*z4^3)]

            sage: # needs sage.rings.finite_rings
            sage: K = GF(3**4)
            sage: [E.ainvs() for E in EllipticCurve(j=K(0)).twists()] # random
            [(0, 0, 0, 1, 0),
             (0, 0, 0, 2, 2*z4^3 + 2*z4^2 + 2*z4 + 2),
             (0, 0, 0, 1, 0),
             (0, 0, 0, 1, 2*z4^3 + 2*z4^2 + 2*z4 + 2),
             (0, 0, 0, z4, 0),
             (0, 0, 0, z4^3, 0)]

        In characteristic 2, the number of twists is 2 except for
        `j=0=1728`, when there are either 3 or 7 depending on whether the
        field has odd or even degree over `\GF{2}`::

            sage: # needs sage.rings.finite_rings
            sage: K = GF(2**7)
            sage: [E.ainvs() for E in EllipticCurve(j=K(1)).twists()]
            [(1, 0, 0, 0, 1), (1, 1, 0, 0, 1)]

            sage: # needs sage.rings.finite_rings
            sage: K = GF(2**7)
            sage: [E.ainvs() for E in EllipticCurve(j=K(0)).twists()]
            [(0, 0, 1, 0, 0), (0, 0, 1, 1, 0), (0, 0, 1, 1, 1)]

            sage: # needs sage.rings.finite_rings
            sage: K = GF(2**8)
            sage: [E.ainvs() for E in EllipticCurve(j=K(1)).twists()] # random
            [(1, 0, 0, 0, 1), (1, z8^7 + z8^6 + z8^5 + z8^4 + z8^2 + z8, 0, 0, 1)]

            sage: # needs sage.rings.finite_rings
            sage: K = GF(2**8)
            sage: [E.ainvs() for E in EllipticCurve(j=K(0)).twists()] # random
            [(0, 0, 1, 0, 0),
             (0, 0, 1, 0, z8^5 + z8^4 + z8^3),
             (0, 0, 1, z8^6 + z8^5 + z8^2 + 1, 0),
             (0, 0, z8^4 + z8^3 + z8^2 + 1, 0, 0),
             (0, 0, z8^4 + z8^3 + z8^2 + 1, 0, z8^3 + z8^2 + 1),
             (0, 0, z8^6 + z8^3 + z8^2, 0, 0),
             (0, 0, z8^6 + z8^3 + z8^2, 0, z8^3 + z8^2)]

        TESTS:

        Randomized check that we find all twists and there are no duplicates::

            sage: # needs sage.rings.finite_rings
            sage: p = next_prime(randrange(2,100))
            sage: e = randrange(1,10)
            sage: F.<t> = GF((p,e))
            sage: E = EllipticCurve(j=F.random_element())
            sage: twists1 = E.twists()
            sage: {sum(E1.is_isomorphic(E2) for E2 in twists1) == 1 for E1 in twists1}
            {True}
            sage: j = E.j_invariant()
            sage: A,B = polygens(F, 'A,B')
            sage: eq = 1728*4*A**3 - j * (4*A**3 + 27*B**2)
            sage: twists2 = []
            sage: for _ in range(10):
            ....:     I = Ideal([eq, A + B - F.random_element()])
            ....:     try:
            ....:         V = I.variety()
            ....:     except ValueError:
            ....:         if I.dimension() == 0:
            ....:              raise
            ....:     if not V:
            ....:         continue
            ....:     sol = choice(V)
            ....:     a, b = sol[A], sol[B]
            ....:     try:
            ....:         twists2.append(EllipticCurve([a, b]))
            ....:     except ArithmeticError:
            ....:         pass
            sage: all(any(E2.is_isomorphic(E1) for E1 in twists1) for E2 in twists2)
            True
        """
        K = self.base_field()
        j = self.j_invariant()
        twists = None
        if not j:
            twists = curves_with_j_0(K)
        elif j == 1728:
            twists = curves_with_j_1728(K)

        if twists:  # i.e. if j=0 or 1728
            # replace the one isomorphic to self with self and move to front
            for i, t in enumerate(twists):
                if self.is_isomorphic(t):
                    twists[i] = twists[0]
                    twists[0] = self
                break
            return twists

        # Now j is not 0 or 1728, and we only have a quadratic twist

        if K.characteristic() == 2:
            # find D with trace 1 for the additive twist
            D = K.one()
            while not D or D.trace() == 0:
                D = K.random_element()
        else:
            # find a nonsquare D.
            D = K.gen()
            q2 = (K.cardinality() - 1) // 2
            while not D or D**q2 == 1:
                D = K.random_element()
            # assert D and D**q2 != 1
            # assert not D.is_square()

        return [self, self.quadratic_twist(D)]


def curves_with_j_0(K):
    r"""
    Return a complete list of pairwise nonisomorphic elliptic curves with `j`-invariant 0 over the finite field `K`.

        .. NOTE::

            In characteristics 2 and 3 this function simply calls ``curves_with_j_0_char2`` or
            ``curves_with_j_0_char3``.  Otherwise there are either 2 or 6 curves, parametrised by
            `K^*/(K^*)^6`.

    Examples:

    For `K=\GF{q}` where `q\equiv1\mod{6}` there are six curves, the sextic twists of `y^2=x^3+1`::

        sage: # needs sage.rings.finite_rings
        sage: from sage.schemes.elliptic_curves.ell_finite_field import curves_with_j_0
        sage: sorted(curves_with_j_0(GF(7)), key = lambda E: E.a_invariants())
        [Elliptic Curve defined by y^2 = x^3 + 1 over Finite Field of size 7,
         Elliptic Curve defined by y^2 = x^3 + 2 over Finite Field of size 7,
         Elliptic Curve defined by y^2 = x^3 + 3 over Finite Field of size 7,
         Elliptic Curve defined by y^2 = x^3 + 4 over Finite Field of size 7,
         Elliptic Curve defined by y^2 = x^3 + 5 over Finite Field of size 7,
         Elliptic Curve defined by y^2 = x^3 + 6 over Finite Field of size 7]
        sage: curves = curves_with_j_0(GF(25)); len(curves)
        6
        sage: all(not curves[i].is_isomorphic(curves[j]) for i in range(6) for j in range(i + 1, 6))
        True
        sage: set(E.j_invariant() for E in curves)
        {0}

    For `K=\GF{q}` where `q\equiv5\mod{6}` there are two curves,
    quadratic twists of each other by `-3`: `y^2=x^3+1` and
    `y^2=x^3-27`::

        sage: from sage.schemes.elliptic_curves.ell_finite_field import curves_with_j_0
        sage: curves_with_j_0(GF(5))
        [Elliptic Curve defined by y^2 = x^3 + 1 over Finite Field of size 5,
         Elliptic Curve defined by y^2 = x^3 + 3 over Finite Field of size 5]
        sage: curves_with_j_0(GF(11))
        [Elliptic Curve defined by y^2 = x^3 + 1 over Finite Field of size 11,
         Elliptic Curve defined by y^2 = x^3 + 6 over Finite Field of size 11]
    """
    if not K.is_finite():
        raise ValueError("field must be finite")
    p = K.characteristic()
    if p == 2:
        return curves_with_j_0_char2(K)
    if p == 3:
        return curves_with_j_0_char3(K)
    q = K.cardinality()
    if q % 3 == 2:
        # Then we only have two quadratic twists (and -3 is non-square)
        return [EllipticCurve(K, [0, a]) for a in [1, -27]]
    # Now we have genuine sextic twists, find D generating K* mod 6th powers
    q2 = (q - 1) // 2
    q3 = (q - 1) // 3
    D = K.gen()
    while not D or D**q2 == 1 or D**q3 == 1:
        D = K.random_element()

    curves = [EllipticCurve(K, [0, D**i]) for i in range(6)]
    # TODO: issue 37110, Precompute orders of sextic twists + docs
    # The idea should be to evaluate the character (D / q) or something
    # Probably reference [RS2010]_ and [Connell1999]_
    # Also a necessary change is `curves_with_j_0` should take in an optional "starting curve"
    # (passed from the original .twists call), because if you start twisting from that curve,
    # then you can also compute the orders!
    return curves


def curves_with_j_1728(K):
    r"""
    Return a complete list of pairwise nonisomorphic elliptic curves with `j`-invariant 1728 over the finite field `K`.

        .. NOTE::

            In characteristics 2 and 3 (so 0=1728) this function simply calls ``curves_with_j_0_char2`` or
            ``curves_with_j_0_char3``.  Otherwise there are either 2 or 4 curves, parametrised by
            `K^*/(K^*)^4`.

    EXAMPLES:

    For `K=\GF{q}` where `q\equiv1\mod{4}`, there are four curves, the quartic twists of `y^2=x^3+x`::

        sage: from sage.schemes.elliptic_curves.ell_finite_field import curves_with_j_1728
        sage: sorted(curves_with_j_1728(GF(5)), key = lambda E: E.a_invariants())
        [Elliptic Curve defined by y^2 = x^3 + x over Finite Field of size 5,
         Elliptic Curve defined by y^2 = x^3 + 2*x over Finite Field of size 5,
         Elliptic Curve defined by y^2 = x^3 + 3*x over Finite Field of size 5,
         Elliptic Curve defined by y^2 = x^3 + 4*x over Finite Field of size 5]
        sage: curves_with_j_1728(GF(49))  # random                                      # needs sage.rings.finite_rings
        [Elliptic Curve defined by y^2 = x^3 + x over Finite Field in z2 of size 7^2,
         Elliptic Curve defined by y^2 = x^3 + z2*x over Finite Field in z2 of size 7^2,
         Elliptic Curve defined by y^2 = x^3 + (z2+4)*x over Finite Field in z2 of size 7^2,
         Elliptic Curve defined by y^2 = x^3 + (5*z2+4)*x over Finite Field in z2 of size 7^2]

    For `K=\GF{q}` where `q\equiv3\mod{4}`, there are two curves,
    quadratic twists of each other by `-1`: `y^2=x^3+x` and
    `y^2=x^3-x`::

        sage: from sage.schemes.elliptic_curves.ell_finite_field import curves_with_j_1728
        sage: curves_with_j_1728(GF(7))
        [Elliptic Curve defined by y^2 = x^3 + x over Finite Field of size 7,
         Elliptic Curve defined by y^2 = x^3 + 6*x over Finite Field of size 7]
        sage: curves_with_j_1728(GF(11))
        [Elliptic Curve defined by y^2 = x^3 + x over Finite Field of size 11,
         Elliptic Curve defined by y^2 = x^3 + 10*x over Finite Field of size 11]
    """
    if not K.is_finite():
        raise ValueError("field must be finite")
    p = K.characteristic()
    if p == 2:
        return curves_with_j_0_char2(K)
    if p == 3:
        return curves_with_j_0_char3(K)
    q = K.cardinality()
    if q % 4 == 3:
        return [EllipticCurve(K, [a,0]) for a in [1,-1]]
    # Now we have genuine quartic twists, find D generating K* mod 4th powers
    q2 = (q - 1) // 2
    D = K.gen()
    while not D or D**q2 == 1:
        D = K.random_element()
    curves = [EllipticCurve(K, [D**i, 0]) for i in range(4)]
    return curves


def curves_with_j_0_char2(K):
    r"""
    Return a complete list of pairwise nonisomorphic elliptic curves with `j`-invariant 0 over the finite field `K` of characteristic 2.

        .. NOTE::

            The number of twists is either 3 or 7 depending on whether
            the field has odd or even degree over `\GF{2}`.  See
            [Connell1999]_, pages 429-431.

    Examples:

    In odd degree, there are three isomorphism classes all with representatives defined over `\GF{2}`::

        sage: from sage.schemes.elliptic_curves.ell_finite_field import curves_with_j_0_char2
        sage: # needs sage.rings.finite_rings
        sage: K = GF(2**7)
        sage: curves = curves_with_j_0_char2(K)
        sage: len(curves)
        3
        sage: [E.ainvs() for E in curves]
        [(0, 0, 1, 0, 0), (0, 0, 1, 1, 0), (0, 0, 1, 1, 1)]

    Check that the curves are mutually non-isomorphic::

        sage: all((e1 == e2 or not e1.is_isomorphic(e2))                                # needs sage.rings.finite_rings
        ....:     for e1 in curves for e2 in curves)
        True

    Check that the weight formula holds::

        sage: sum(1/len(E.automorphisms()) for E in curves) == 1                        # needs sage.rings.finite_rings
        True

    In even degree there are seven isomorphism classes::

        sage: from sage.schemes.elliptic_curves.ell_finite_field import curves_with_j_0_char2
        sage: # needs sage.rings.finite_rings
        sage: K = GF(2**8)
        sage: curves = EllipticCurve(j=K(0)).twists()
        sage: len(curves)
        7
        sage: [E.ainvs() for E in curves] # random
        [(0, 0, 1, 0, 0),
         (0, 0, 1, 0, z8^5 + z8^4 + z8^3),
         (0, 0, 1, z8^6 + z8^5 + z8^2 + 1, 0),
         (0, 0, z8^4 + z8^3 + z8^2 + 1, 0, 0),
         (0, 0, z8^4 + z8^3 + z8^2 + 1, 0, z8^3 + z8^2 + 1),
         (0, 0, z8^6 + z8^3 + z8^2, 0, 0),
         (0, 0, z8^6 + z8^3 + z8^2, 0, z8^3 + z8^2)]

    Check that the twists are mutually non-isomorphic::

        sage: all((e1 == e2 or not e1.is_isomorphic(e2))                                # needs sage.rings.finite_rings
        ....:     for e1 in curves for e2 in curves)
        True

    Check that the weight formula holds::

        sage: sum(1/len(E.automorphisms()) for E in curves) == 1                        # needs sage.rings.finite_rings
        True
    """
    if not K.is_finite() or K.characteristic() != 2:
        raise ValueError("field must be finite of characteristic 2")
    if K.degree() % 2:
        return [EllipticCurve(K, [0, 0, 1, 0, 0]),
                EllipticCurve(K, [0, 0, 1, 1, 0]),
                EllipticCurve(K, [0, 0, 1, 1, 1])]
    # find a,b,c,d,e such that
    # a is not a cube, i.e. a**((q-1)//3)!=1
    # Tr(b)=1
    # X^4+X+c irreducible
    # X^2+a*X+d irreducible
    # X^2+a^2*X+e irreducible
    a = b = c = d = e = None
    x = polygen(K)
    q3 = (K.cardinality()-1)//3
    while not a or a**q3 == 1:
        a = K.random_element()
    asq = a*a
    while not b or not b.trace():
        b = K.random_element()
    c = K.one() # OK if degree is 2 mod 4
    if K.degree() % 4 == 0:
        while (x**4+x+c).roots():
            c = K.random_element()
    while not d or (x**2+a*x+d).roots():
        d = K.random_element()
    while not e or (x**2+asq*x+e).roots():
        e = K.random_element()
    return [EllipticCurve(K, ai) for ai in
            [[0,0,1,0,0], [0,0,1,0,b], [0,0,1,c,0], [0,0,a,0,0], [0,0,a,0,d], [0,0,asq,0,0], [0,0,asq,0,e]]]


def curves_with_j_0_char3(K):
    r"""
    Return a complete list of pairwise nonisomorphic elliptic curves with `j`-invariant 0 over the finite field `K` of characteristic 3.

        .. NOTE::

            The number of twists is either 4 or 6 depending on whether
            the field has odd or even degree over `\GF{3}`.  See
            [Connell1999]_, pages 429-431.

    Examples:

    In odd degree, there are four isomorphism classes::

        sage: from sage.schemes.elliptic_curves.ell_finite_field import curves_with_j_0_char3
        sage: # needs sage.rings.finite_rings
        sage: K = GF(3**5)
        sage: curves = curves_with_j_0_char3(K)
        sage: len(curves)
        4
        sage: [E.ainvs() for E in curves] # random
        [(0, 0, 0, 1, 0),
         (0, 0, 0, 2, 0),
         (0, 0, 0, 2, z5^4 + z5^3 + z5^2),
         (0, 0, 0, 2, 2*z5^4 + 2*z5^3 + 2*z5^2)]

    Check that the twists are mutually non-isomorphic::

        sage: all((e1 == e2 or not e1.is_isomorphic(e2))                                # needs sage.rings.finite_rings
        ....:     for e1 in curves for e2 in curves)
        True

    Check that the weight formula holds::

        sage: sum(1/len(E.automorphisms()) for E in curves) == 1                        # needs sage.rings.finite_rings
        True

    In even degree, there are six isomorphism classes::

        sage: from sage.schemes.elliptic_curves.ell_finite_field import curves_with_j_0_char3
        sage: # needs sage.rings.finite_rings
        sage: K = GF(3**4)
        sage: curves = EllipticCurve(j=K(0)).twists()
        sage: len(curves)
        6
        sage: [E.ainvs() for E in curves] # random
        [(0, 0, 0, 1, 0),
         (0, 0, 0, 2, 2*z4^3 + 2*z4^2 + 2*z4 + 2),
         (0, 0, 0, 1, 0),
         (0, 0, 0, 1, 2*z4^3 + 2*z4^2 + 2*z4 + 2),
         (0, 0, 0, z4, 0),
         (0, 0, 0, z4^3, 0)]

    Check that the twists are mutually non-isomorphic::

        sage: all((e1 == e2 or not e1.is_isomorphic(e2))                                # needs sage.rings.finite_rings
        ....:     for e1 in curves for e2 in curves)
        True

    Check that the weight formula holds::

        sage: sum(1/len(E.automorphisms()) for E in curves) == 1                        # needs sage.rings.finite_rings
        True
    """
    if not K.is_finite() or K.characteristic() != 3:
        raise ValueError("field must be finite of characteristic 3")
    # find b with nonzero trace
    b = None
    while not b or not b.trace():
        b = K.random_element()

    if K.degree() % 2:
        return [EllipticCurve(K, a4a6) for a4a6 in
            [[1,0], [-1,0], [-1,b], [-1,-b]]]

    # find a, i, c where:
    # a generates K* mod 4th powers, i.e. non-square,
    # i^2=-1
    # c with x^3+a^2*x+c irreducible
    a = K.gen()
    q2 = (K.cardinality()-1)//2
    while not a or a**q2 == 1:
        a = K.random_element()
    x = polygen(K)
    i = (x**2+1).roots()[0][0]
    c = None
    while not c or (x**3 + a**2*x + c).roots():
        c = K.random_element()
    return [EllipticCurve(K, a4a6) for a4a6 in
            [[1,0], [1,i*b], [a,0], [a**2,0], [a**2,c], [a**3,0]]]

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

    - ``p`` -- integer; a prime number

    - ``use_cache`` -- boolean (default: ``True``); use cached coefficients if they exist

    ALGORITHM:

    First compute H(X) whose roots are the Legendre
    `\lambda`-invariants of supersingular curves (Silverman V.4.1(b))
    in characteristic `p`.  Then, using a resultant computation with
    the polynomial relating `\lambda` and `j` (Silverman III.1.7(b)),
    we recover the polynomial (in variable ``j``) whose roots are the
    `j`-invariants.  Factors of `j` and `j-1728` are removed if
    present.

    .. NOTE::

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
    if p < 13:
        return J.parent().one()
    if use_cache:
        fill_ss_j_dict()
        if p in supersingular_j_polynomials:
            return J.parent()(supersingular_j_polynomials[p])

    from sage.misc.misc_c import prod
    m = (p-1)//2
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
    Return ``True`` if `j` is a supersingular `j`-invariant.

    INPUT:

    - ``j`` -- finite field element

    - ``proof``-- boolean (default: ``True``); if ``True``, returns a proved
      result.  If ``False``, then a return value of ``False`` is certain but a
      return value of ``True`` may be based on a probabilistic test.  See
      the ALGORITHM section below for more details.

    OUTPUT: boolean; ``True`` if `j` is supersingular, else ``False``

    ALGORITHM:

    For small characteristics `p` we check whether the `j`-invariant
    is in a precomputed list of supersingular values.  Otherwise we
    next check the `j`-invariant.  If `j=0`, the curve is
    supersingular if and only if `p=2` or `p\equiv3\pmod{4}`; if
    `j=1728`, the curve is supersingular if and only if `p=3` or
    `p\equiv2\pmod{3}`.  Next, if the base field is the prime field
    `{\rm GF}(p)`, we check that `(p+1)P=0` for several random points
    `P`, returning ``False`` if any fail: supersingular curves over `{\rm
    GF}(p)` have cardinality `p+1`.  If Proof is false we now return
    ``True``.  Otherwise we compute the cardinality and return ``True`` if and
    only if it is divisible by `p`.

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.ell_finite_field import is_j_supersingular, supersingular_j_polynomials
        sage: [(p,[j for j in GF(p) if is_j_supersingular(j)]) for p in prime_range(30)]
        [(2, [0]), (3, [0]), (5, [0]), (7, [6]), (11, [0, 1]), (13, [5]),
         (17, [0, 8]), (19, [7, 18]), (23, [0, 3, 19]), (29, [0, 2, 25])]

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
    if not (isinstance(j, Element) and isinstance(j.parent(), FiniteField)):
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


def special_supersingular_curve(F, *, endomorphism=False):
    r"""
    Given a finite field ``F``, construct a "special" supersingular
    elliptic curve `E` defined over ``F``.

    Such a curve

    - has coefficients in `\mathbb F_p`;

    - has group structure `E(\mathbb F_p) \cong \ZZ/(p+1)` and
      `E(\mathbb F_{p^2}) \cong \ZZ/(p+1) \times \ZZ/(p+1)`;

    - has an endomorphism `\vartheta` of small degree `q` that
      anticommutes with the `\mathbb F_p`-Frobenius on `E`.

    (The significance of `\vartheta` is that any such endomorphism,
    together with the `\mathbb F_p`-Frobenius, generates the endomorphism
    algebra `\mathrm{End}(E) \otimes \QQ`.)

    INPUT:

    - ``F`` -- finite field `\mathbb F_{p^r}`;

    - ``endomorphism`` -- boolean (default: ``False``); when set to ``True``,
      it is required that `2 \mid r`, and the function then additionally
      returns `\vartheta`

    EXAMPLES::

        sage: special_supersingular_curve(GF(1013^2), endomorphism=True)
        (Elliptic Curve defined by y^2 = x^3 + 1 over Finite Field in z2 of size 1013^2,
         Isogeny of degree 3 from Elliptic Curve defined by y^2 = x^3 + 1 over Finite Field in z2 of size 1013^2 to Elliptic Curve defined by y^2 = x^3 + 1 over Finite Field in z2 of size 1013^2)

        sage: special_supersingular_curve(GF(1019^2), endomorphism=True)
        (Elliptic Curve defined by y^2 = x^3 + x over Finite Field in z2 of size 1019^2,
         Elliptic-curve endomorphism of Elliptic Curve defined by y^2 = x^3 + x over Finite Field in z2 of size 1019^2
           Via:  (u,r,s,t) = (389*z2 + 241, 0, 0, 0))

        sage: special_supersingular_curve(GF(1021^2), endomorphism=True)
        (Elliptic Curve defined by y^2 = x^3 + 785*x + 794 over Finite Field in z2 of size 1021^2,
         Isogeny of degree 2 from Elliptic Curve defined by y^2 = x^3 + 785*x + 794 over Finite Field in z2 of size 1021^2 to Elliptic Curve defined by y^2 = x^3 + 785*x + 794 over Finite Field in z2 of size 1021^2)

        sage: special_supersingular_curve(GF(1031^2), endomorphism=True)
        (Elliptic Curve defined by y^2 = x^3 + x over Finite Field in z2 of size 1031^2,
         Elliptic-curve endomorphism of Elliptic Curve defined by y^2 = x^3 + x over Finite Field in z2 of size 1031^2
           Via:  (u,r,s,t) = (747*z2 + 284, 0, 0, 0))

        sage: special_supersingular_curve(GF(1033^2), endomorphism=True)
        (Elliptic Curve defined by y^2 = x^3 + 53*x + 980 over Finite Field in z2 of size 1033^2,
         Isogeny of degree 11 from Elliptic Curve defined by y^2 = x^3 + 53*x + 980 over Finite Field in z2 of size 1033^2 to Elliptic Curve defined by y^2 = x^3 + 53*x + 980 over Finite Field in z2 of size 1033^2)

        sage: special_supersingular_curve(GF(1039^2), endomorphism=True)
        (Elliptic Curve defined by y^2 = x^3 + x over Finite Field in z2 of size 1039^2,
         Elliptic-curve endomorphism of Elliptic Curve defined by y^2 = x^3 + x over Finite Field in z2 of size 1039^2
           Via:  (u,r,s,t) = (626*z2 + 200, 0, 0, 0))

        sage: special_supersingular_curve(GF(1049^2), endomorphism=True)
        (Elliptic Curve defined by y^2 = x^3 + 1 over Finite Field in z2 of size 1049^2,
         Isogeny of degree 3 from Elliptic Curve defined by y^2 = x^3 + 1 over Finite Field in z2 of size 1049^2 to Elliptic Curve defined by y^2 = x^3 + 1 over Finite Field in z2 of size 1049^2)

        sage: special_supersingular_curve(GF(1051^2), endomorphism=True)
        (Elliptic Curve defined by y^2 = x^3 + x over Finite Field in z2 of size 1051^2,
         Elliptic-curve endomorphism of Elliptic Curve defined by y^2 = x^3 + x over Finite Field in z2 of size 1051^2
           Via:  (u,r,s,t) = (922*z2 + 129, 0, 0, 0))

    TESTS::

        sage: p = random_prime(1000)
        sage: E = special_supersingular_curve(GF(p))
        sage: E.is_supersingular()
        True
        sage: E.order() == p + 1
        True
        sage: F.<t> = GF((p,2))
        sage: E, endo = special_supersingular_curve(F, endomorphism=True)
        sage: E.is_supersingular()
        True
        sage: E.j_invariant() in GF(p)
        True
        sage: E.abelian_group().invariants() == (p+1, p+1)
        True
        sage: endo.domain() is endo.codomain() is E
        True
        sage: endo.trace()
        0
        sage: pi = E.frobenius_isogeny()
        sage: pi.codomain() is pi.domain() is E
        True
        sage: pi * endo == -endo * pi
        True

    Also try it for larger-degree fields::

        sage: k = ZZ(randrange(3, 10, 2))
        sage: E = special_supersingular_curve(GF((p, k)))
        sage: E.is_supersingular()
        True
        sage: F.<t> = GF((p, 2*k))
        sage: E, endo = special_supersingular_curve(F, endomorphism=True)
        sage: E.is_supersingular()
        True
        sage: E.j_invariant() in GF(p)
        True
        sage: endo.domain() is endo.codomain() is E
        True
        sage: endo.trace()
        0
        sage: pi = E.frobenius_isogeny()
        sage: pi.codomain() is pi.domain() is E
        True
        sage: pi * endo == -endo * pi
        True

    .. NOTE::

        This function makes no guarantees about the distribution of
        the output. The current implementation is deterministic in
        many cases.

    ALGORITHM: [Bro2009]_, Algorithm 2.4
    """
    if not isinstance(F, FiniteField):
        raise TypeError('input must be a finite field')
    p = F.characteristic()
    deg = F.degree()

    if endomorphism and deg % 2:
        raise ValueError('endomorphism was requested but is not defined over given field')

    E = None

    # first find the degree q of our special endomorphism
    if p == 2:
        q = 3
        E = EllipticCurve(F, [0,0,1,0,0])

    elif p % 4 == 3:
        q = 1
        E = EllipticCurve(F, [1,0])

    elif p % 3 == 2:
        q = 3
        E = EllipticCurve(F, [0,1])

    elif p % 8 == 5:
        q = 2
        E = EllipticCurve(F, [-4320, 96768])

    else:
        from sage.arith.misc import legendre_symbol
        for q in map(ZZ, range(3,p,4)):
            if not q.is_prime():
                continue
            if legendre_symbol(-q, p) == -1:
                break
        else:
            assert False  # should never happen

    if E is None:
        from sage.arith.misc import fundamental_discriminant
        from sage.schemes.elliptic_curves.cm import hilbert_class_polynomial
        H = hilbert_class_polynomial(fundamental_discriminant(-q))
        j = H.change_ring(GF(p)).any_root()
        a = 27 * j / (4 * (1728-j))
        E = EllipticCurve(F, [a,-a])

    if ZZ(2).divides(deg):
        k = deg//2
        E.set_order((p**k - (-1)**k)**2)
    else:
        E.set_order(p**deg - (-1)**deg)

    if not endomorphism:
        return E

    if q == 1 or p <= 13:
        if q == 1:
            endos = E.automorphisms()
        else:
            endos = (iso*phi for phi in E.isogenies_prime_degree(q)
                             for iso in phi.codomain().isomorphisms(E))
        endo = next(endo for endo in endos if endo.trace().is_zero())

    else:
        from sage.schemes.elliptic_curves.weierstrass_morphism import WeierstrassIsomorphism
        iso = WeierstrassIsomorphism(None, (F(-q).sqrt(),0,0,0), E)
        if q == 3 and E.a_invariants() == (0,0,0,0,1):
            # workaround for #21883
            endo = E.isogeny(E(0,1))
        else:
            endo = E.isogeny(None, iso.domain(), degree=q)
        endo = iso * endo

    endo._degree = ZZ(q)
    endo.trace.set_cache(ZZ.zero())
    return E, endo


def EllipticCurve_with_order(m, *, D=None):
    r"""
    Return an iterator for elliptic curves over finite fields with the given order. The curves are
    computed using the Complex Multiplication (CM) method.

    A `:sage:`~sage.structure.factorization.Factorization` can be passed for ``m``, in which case
    the algorithm is more efficient.

    If ``D`` is specified, it is used as the discriminant.

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.ell_finite_field import EllipticCurve_with_order
        sage: E = next(EllipticCurve_with_order(1234)); E  # random
        Elliptic Curve defined by y^2 = x^3 + 1142*x + 1209 over Finite Field of size 1237
        sage: E.order() == 1234
        True

    When ``iter`` is set, the function returns an iterator of all elliptic curves with the given
    order::

        sage: from sage.schemes.elliptic_curves.ell_finite_field import EllipticCurve_with_order
        sage: it = EllipticCurve_with_order(21); it
        <generator object EllipticCurve_with_order at 0x...>
        sage: E = next(it); E  # random
        Elliptic Curve defined by y^2 = x^3 + 6*x + 14 over Finite Field of size 23
        sage: E.order() == 21
        True
        sage: Es = [E] + list(it); Es  # random
        [Elliptic Curve defined by y^2 = x^3 + 6*x + 14 over Finite Field of size 23,
         Elliptic Curve defined by y^2 = x^3 + 12*x + 4 over Finite Field of size 23,
         Elliptic Curve defined by y^2 = x^3 + 5*x + 2 over Finite Field of size 23,
         Elliptic Curve defined by y^2 = x^3 + (z2+3) over Finite Field in z2 of size 5^2,
         Elliptic Curve defined by y^2 = x^3 + (2*z2+2) over Finite Field in z2 of size 5^2,
         Elliptic Curve defined by y^2 = x^3 + 7*x + 1 over Finite Field of size 19,
         Elliptic Curve defined by y^2 = x^3 + 17*x + 10 over Finite Field of size 19,
         Elliptic Curve defined by y^2 = x^3 + 5*x + 12 over Finite Field of size 17,
         Elliptic Curve defined by y^2 = x^3 + 9*x + 1 over Finite Field of size 17,
         Elliptic Curve defined by y^2 = x^3 + 7*x + 6 over Finite Field of size 17,
         Elliptic Curve defined by y^2 = x^3 + z3^2*x^2 + (2*z3^2+z3) over Finite Field in z3 of size 3^3,
         Elliptic Curve defined by y^2 = x^3 + (z3^2+2*z3+1)*x^2 + (2*z3^2+2*z3) over Finite Field in z3 of size 3^3,
         Elliptic Curve defined by y^2 = x^3 + (z3^2+z3+1)*x^2 + (2*z3^2+1) over Finite Field in z3 of size 3^3,
         Elliptic Curve defined by y^2 + (z4^2+z4+1)*y = x^3 over Finite Field in z4 of size 2^4,
         Elliptic Curve defined by y^2 + (z4^2+z4)*y = x^3 over Finite Field in z4 of size 2^4,
         Elliptic Curve defined by y^2 = x^3 + 18*x + 26 over Finite Field of size 29,
         Elliptic Curve defined by y^2 = x^3 + 11*x + 19 over Finite Field of size 29,
         Elliptic Curve defined by y^2 = x^3 + 4 over Finite Field of size 19,
         Elliptic Curve defined by y^2 = x^3 + 19 over Finite Field of size 31,
         Elliptic Curve defined by y^2 = x^3 + 4 over Finite Field of size 13]
        sage: all(E.order() == 21 for E in Es)
        True

    Indeed, we can verify that this is correct. Hasse's bounds tell us that
    `p \leq 50` (approximately), and the rest can be checked via bruteforce::

        sage: for p in prime_range(50):
        ....:     for j in range(p):
        ....:         E0 = EllipticCurve(GF(p), j=j)
        ....:         for Et in E0.twists():
        ....:             if Et.order() == 21:
        ....:                 assert any(Et.is_isomorphic(E) for E in Es)

    .. NOTE::

        The output curves are not deterministic, as :func:`EllipticCurve_finite_field.twists` is not
        deterministic. However, the order of the j-invariants and base fields is fixed.

    AUTHORS:

     - Gareth Ma and Giacomo Pope (Sage Days 123): initial version
    """
    from sage.arith.misc import is_prime_power, factor
    from sage.quadratic_forms.binary_qf import BinaryQF
    from sage.structure.factorization import Factorization
    from sage.schemes.elliptic_curves.cm import hilbert_class_polynomial

    def find_q(m, m4_fac, D):
        for t, _ in BinaryQF(1, 0, -D).solve_integer(m4_fac, _flag=3):
            yield m + 1 - t
            yield m + 1 + t

    if isinstance(m, Factorization):
        m4_fac = m * factor(4)
        m_val = m.value()
    else:
        m4_fac = factor(m * 4)
        m_val = m

    if D is None:
        Ds = (D for D in range(-4 * m_val, 0) if D % 4 in [0, 1])
    else:
        assert D < 0 and D % 4 in [0, 1]
        Ds = [D]

    seen = set()
    for D in Ds:
        for q in find_q(m_val, m4_fac, D):
            if not is_prime_power(q):
                continue

            H = hilbert_class_polynomial(D)
            K = GF(q)
            roots = H.roots(ring=K)
            for j0, _ in roots:
                E = EllipticCurve(j=j0)
                for Et in E.twists():
                    if any(Et.is_isomorphic(E) for E in seen):
                        continue
                    try:
                        # This tests whether the curve has given order
                        Et.set_order(m_val)
                        seen.add(Et)
                        yield Et
                    except ValueError:
                        pass

def EllipticCurve_with_prime_order(N):
    r"""
    Given a prime number ``N``, find another prime number `p` and construct an
    elliptic curve `E` defined over `\mathbb F_p` such that
    `\#E(\mathbb F_p) = N`.

    INPUT:

    - ``N`` -- integer; the order for which we seek an elliptic curve. Must be a
      prime number.

    OUTPUT: an elliptic curve `E/\mathbb F_p` of order ``N``

    EXAMPLES::

        sage: N = next_prime(int(b'sagemath'.hex(), 16))
        sage: E = EllipticCurve_with_prime_order(N)
        sage: E
        Elliptic Curve defined by y^2 = x^3 + 4757897140353078952*x +
        1841350074072114366 over Finite Field of size 8314040074357871443
        sage: E.order() == N
        True

    The execution time largely depends on the input::

        sage: N = 125577861263605878504082476745517446213
        sage: E = EllipticCurve_with_prime_order(N) # Takes ~1 second.
        sage: E.order() == N
        True

    TESTS::

        sage: for N in prime_range(3, 100):
        ....:     E = EllipticCurve_with_prime_order(N)
        ....:     assert E.order() == N

        sage: N = 15175980689839334471
        sage: E = EllipticCurve_with_prime_order(N)
        sage: E.order() == N
        True

        sage: N = next_prime(123456789)
        sage: E = EllipticCurve_with_prime_order(N)
        sage: E.order() == N
        True

        sage: N = 123456789
        sage: E = EllipticCurve_with_prime_order(N)
        Traceback (most recent call last):
        ...
        ValueError: input order is not prime

        sage: E = EllipticCurve_with_prime_order(0)
        Traceback (most recent call last):
        ...
        ValueError: input order is not prime

        sage: E = EllipticCurve_with_prime_order(-7)
        Traceback (most recent call last):
        ...
        ValueError: input order is not prime

    .. NOTE::

        Depending on the input, this function may run for a *very* long time.
        This algorithm consists of multiple "search rounds" for a suitable
        discriminant `D`. We expect this algorithm to terminate after a number
        of rounds that is polynomial in `loglog N`. In practice (cf. Section 5),
        this number is usually 1.

    ALGORITHM: [BS2007]_, Algorithm 2.2
    """
    from sage.arith.misc import is_prime
    from sage.combinat.subset import powerset
    from sage.functions.other import ceil
    from sage.misc.functional import symbolic_prod as product
    from sage.quadratic_forms.binary_qf import BinaryQF
    from sage.rings.fast_arith import prime_range
    from sage.schemes.elliptic_curves.cm import hilbert_class_polynomial

    if not is_prime(N):
        raise ValueError("input order is not prime")

    # The algorithm consists of multiple ‘search rounds’ for a suitable
    # discriminant `D`, `r` defines the number of rounds. We expect this
    # algorithm to terminate after a number of rounds that is polynomial in
    # loglog N.
    r = 0

    while True:
        # Iterating over the odd primes by chunks of size log(`N`).
        S = prime_range(3, ceil((r + 1) * N.log()))

        # Every possible products of distinct elements of `S`.
        # There probably is a more optimal way to compute all possible products
        # of elements of S than using a powerset. Here many multiplications are
        # done multiple times.
        for e in powerset(S):
            D = product(e)
            if -D % 8 != 5:
                continue

            Q = BinaryQF([1, 0, D])
            sol = Q.solve_integer(4 * N)
            if sol is None:
                continue

            x, _ = sol
            p1 = N + 1 - x
            p2 = N + 1 + x
            for p_i in [p1, p2]:
                if is_prime(p_i):
                    H = hilbert_class_polynomial(-D)
                    for j0, _ in H.roots(ring=GF(p_i)):
                        E = EllipticCurve(j=j0)
                        if E.order() == N:
                            return E
                        else:
                            for Et in E.twists():
                                if Et.order() == N:
                                    return Et

        # At this point, no discriminant has been found, moving to next round
        # and extending the prime list.
        r += 1
