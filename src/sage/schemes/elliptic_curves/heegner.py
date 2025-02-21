# sage.doctest: needs sage.rings.number_field
r"""
Heegner points on elliptic curves over the rational numbers

AUTHORS:

- William Stein (August 2009)-- most of the initial version

- Robert Bradshaw (July 2009) -- an early version of some specific code

EXAMPLES::

    sage: E = EllipticCurve('433a')
    sage: P = E.heegner_point(-8,3)
    sage: z = P.point_exact(201); z
    (-4/3 : 1/9*a : 1)
    sage: parent(z)
    Abelian group of points on Elliptic Curve defined by y^2 + x*y = x^3 + 1
     over Number Field in a with defining polynomial x^2 - 12*x + 111
    sage: parent(z[0]).discriminant()
    -3
    sage: E.quadratic_twist(-3).rank()
    1
    sage: K.<a> = QuadraticField(-8)
    sage: K.factor(3)
    (Fractional ideal (1/2*a + 1)) * (Fractional ideal (-1/2*a + 1))

Next try an inert prime::

    sage: K.factor(5)
    Fractional ideal (5)
    sage: P = E.heegner_point(-8,5)
    sage: z = P.point_exact(300)
    sage: z[0].charpoly().factor()
    (x^6 + x^5 - 1/4*x^4 + 19/10*x^3 + 31/20*x^2 - 7/10*x + 49/100)^2
    sage: z[1].charpoly().factor()
    x^12 - x^11 + 6/5*x^10 - 33/40*x^9 - 89/320*x^8 + 3287/800*x^7
     - 5273/1600*x^6 + 993/4000*x^5 + 823/320*x^4 - 2424/625*x^3
     + 12059/12500*x^2 + 3329/25000*x + 123251/250000
    sage: f = P.x_poly_exact(300); f
    x^6 + x^5 - 1/4*x^4 + 19/10*x^3 + 31/20*x^2 - 7/10*x + 49/100
    sage: f.discriminant().factor()
    -1 * 2^-9 * 5^-9 * 7^2 * 281^2 * 1021^2

We find some Mordell-Weil generators in the rank 1 case using Heegner points::

    sage: E = EllipticCurve('43a'); P = E.heegner_point(-7)
    sage: P.x_poly_exact()
    x
    sage: z = P.point_exact(); z == E(0,0,1) or -z == E(0,0,1)
    True

    sage: E = EllipticCurve('997a')
    sage: E.rank()
    1
    sage: E.heegner_discriminants_list(10)
    [-19, -23, -31, -35, -39, -40, -52, -55, -56, -59]
    sage: P = E.heegner_point(-19)
    sage: P.x_poly_exact()
    x - 141/49
    sage: z = P.point_exact(); z == E(141/49, -162/343, 1)  or -z == E(141/49, -162/343, 1)
    True

Here we find that the Heegner point generates a subgroup of index 3::

    sage: E = EllipticCurve('92b1')
    sage: E.heegner_discriminants_list(1)
    [-7]
    sage: P = E.heegner_point(-7)
    sage: z = P.point_exact(); z == E(0, 1, 1)  or -z == E(0, 1, 1)
    True
    sage: E.regulator()
    0.0498083972980648
    sage: z.height()
    0.448275575682583
    sage: P = E(1,1); P # a generator
    (1 : 1 : 1)
    sage: -3*P
    (0 : 1 : 1)
    sage: E.tamagawa_product()
    3

The above is consistent with the following analytic computation::

    sage: E.heegner_index(-7)
    3.0000?
"""

# ****************************************************************************
#       Copyright (C) 2005-2009 William Stein <wstein@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from itertools import product

import sage.rings.abc

from sage.arith.functions import lcm
from sage.arith.misc import (binomial, factorial, prime_divisors,
                             GCD as gcd, XGCD as xgcd)
from sage.matrix.constructor import matrix
from sage.matrix.matrix_space import MatrixSpace
from sage.misc.cachefunc import cached_method
from sage.misc.misc_c import prod
from sage.misc.lazy_import import lazy_import
from sage.misc.verbose import verbose
from sage.modular.modsym.p1list import P1List
from sage.quadratic_forms.binary_qf import BinaryQF, BinaryQF_reduced_representatives
from sage.rings.factorint import factor_trial_division
from sage.rings.fast_arith import prime_range
from sage.rings.finite_rings.finite_field_constructor import FiniteField as GF
from sage.rings.finite_rings.integer_mod_ring import IntegerModRing as Integers
from sage.rings.infinity import Infinity as infinity
from sage.rings.integer_ring import ZZ
from sage.rings.number_field.number_field_element_base import NumberFieldElement_base
from sage.rings.rational_field import QQ
from sage.structure.sage_object import SageObject
from sage.structure.richcmp import (richcmp_method, richcmp,
                                    richcmp_not_equal, rich_to_bool)

lazy_import('sage.rings.complex_double', 'CDF')
lazy_import('sage.rings.complex_mpfr', 'ComplexField')
lazy_import('sage.rings.number_field.number_field', ['NumberField', 'QuadraticField'])
lazy_import('sage.rings.real_mpfi', 'RealIntervalField')
lazy_import('sage.rings.real_mpfr', 'RealField')

###############################################################################
#
# The exported functions, which are in most cases enough to get the
# user going working with Heegner points:
#
#    heegner_points -- all of them with given level, discriminant, conductor
#    heegner_point -- a specific one
#
###############################################################################


def heegner_points(N, D=None, c=None):
    """
    Return all Heegner points of given level `N`.  Can also restrict
    to Heegner points with specified discriminant `D` and optionally
    conductor `c`.

    INPUT:

    - ``N`` -- level (positive integer)

    - ``D`` -- discriminant (negative integer)

    - ``c`` -- conductor (positive integer)

    EXAMPLES::

        sage: heegner_points(389, -7)
        Set of all Heegner points on X_0(389) associated to QQ[sqrt(-7)]
        sage: heegner_points(389, -7, 1)
        All Heegner points of conductor 1 on X_0(389) associated to QQ[sqrt(-7)]
        sage: heegner_points(389, -7, 5)
        All Heegner points of conductor 5 on X_0(389) associated to QQ[sqrt(-7)]
    """
    if D is None and c is None:
        return HeegnerPoints_level(N)
    if D is not None and c is None:
        return HeegnerPoints_level_disc(N, D)
    if D is not None and c is not None:
        return HeegnerPoints_level_disc_cond(N, D, c)
    raise TypeError


def heegner_point(N, D=None, c=1):
    """
    Return a specific Heegner point of level `N` with given
    discriminant and conductor.  If `D` is not specified, then the
    first valid Heegner discriminant is used.  If `c` is not given,
    then `c=1` is used.

    INPUT:

    - ``N`` -- level (positive integer)

    - ``D`` -- discriminant (optional: default first valid `D`)

    - ``c`` -- conductor (positive integer, default: 1)

    EXAMPLES::

        sage: heegner_point(389)
        Heegner point 1/778*sqrt(-7) - 185/778 of discriminant -7 on X_0(389)
        sage: heegner_point(389,-7)
        Heegner point 1/778*sqrt(-7) - 185/778 of discriminant -7 on X_0(389)
        sage: heegner_point(389,-7,5)
        Heegner point 5/778*sqrt(-7) - 147/778 of discriminant -7 and conductor 5 on X_0(389)
        sage: heegner_point(389,-20)
        Heegner point 1/778*sqrt(-20) - 165/389 of discriminant -20 on X_0(389)
    """
    if D is not None:
        return heegner_points(N, D, c)[0]
    H = heegner_points(N)
    D = H.discriminants(1)[0]
    return heegner_points(N, D, c)[0]


###############################################################################
#
# Ring class fields, represented as abstract objects.  These do not
# derive from number fields, since we do not need to work with their
# elements, and explicitly representing them as number fields would be
# far too difficult.
#
###############################################################################

class RingClassField(SageObject):
    """
    A Ring class field of a quadratic imaginary field of given conductor.

    .. NOTE::

        This is a *ring* class field, not a ray class field. In
        general, the ring class field of given conductor is a subfield
        of the ray class field of the same conductor.

    EXAMPLES::

        sage: heegner_point(37,-7).ring_class_field()
        Hilbert class field of QQ[sqrt(-7)]
        sage: heegner_point(37,-7,5).ring_class_field()
        Ring class field extension of QQ[sqrt(-7)] of conductor 5
        sage: heegner_point(37,-7,55).ring_class_field()
        Ring class field extension of QQ[sqrt(-7)] of conductor 55

    TESTS::

        sage: K_c = heegner_point(37,-7).ring_class_field()
        sage: type(K_c)
        <class 'sage.schemes.elliptic_curves.heegner.RingClassField'>
        sage: loads(dumps(K_c)) == K_c
        True
    """
    def __init__(self, D, c, check=True):
        """
        INPUT:

        - ``D`` -- discriminant of quadratic imaginary field

        - ``c`` -- conductor (positive integer coprime to `D`)

        - ``check`` -- boolean (default: ``True``); whether to check
          validity of input

        EXAMPLES::

            sage: sage.schemes.elliptic_curves.heegner.RingClassField(-7, 5, False)
            Ring class field extension of QQ[sqrt(-7)] of conductor 5
        """
        if check:
            D = ZZ(D)
            c = ZZ(c)
        self.__D = D
        self.__c = c

    def __eq__(self, other):
        """
        Used for equality testing.

        EXAMPLES::

            sage: E = EllipticCurve('389a')
            sage: K5 = E.heegner_point(-7,5).ring_class_field()
            sage: K11 = E.heegner_point(-7,11).ring_class_field()
            sage: K5 == K11
            False
            sage: K5 == K5
            True
            sage: K11 == 11
            False
        """
        return isinstance(other, RingClassField) and self.__D == other.__D and self.__c == other.__c

    def __ne__(self, other):
        """
        Check whether ``self`` is not equal to ``other``.

        EXAMPLES::

            sage: E = EllipticCurve('389a')
            sage: K5 = E.heegner_point(-7,5).ring_class_field()
            sage: K11 = E.heegner_point(-7,11).ring_class_field()
            sage: K5 != K11
            True
            sage: K5 != K5
            False
            sage: K11 != 11
            True
        """
        return not (self == other)

    def __hash__(self):
        """
        Used for computing hash of ``self``.

        .. NOTE::

            The hash is equal to the hash of the pair
            ``(discriminant, conductor)``.

        EXAMPLES::

            sage: E = EllipticCurve('389a'); K5 = E.heegner_point(-7,5).ring_class_field()
            sage: hash(K5) == hash((-7,5))
            True
        """
        return hash((self.__D, self.__c))

    def conductor(self):
        """
        Return the conductor of this ring class field.

        EXAMPLES::

            sage: E = EllipticCurve('389a'); K5 = E.heegner_point(-7,5).ring_class_field()
            sage: K5.conductor()
            5
        """
        return self.__c

    def discriminant_of_K(self):
        """
        Return the discriminant of the quadratic imaginary field `K` contained in ``self``.

        EXAMPLES::

            sage: E = EllipticCurve('389a'); K5 = E.heegner_point(-7,5).ring_class_field()
            sage: K5.discriminant_of_K()
            -7
        """
        return self.__D

    @cached_method
    def ramified_primes(self):
        r"""
        Return the primes of `\ZZ` that ramify in this ring class field.

        EXAMPLES::

            sage: E = EllipticCurve('389a'); K55 = E.heegner_point(-7,55).ring_class_field()
            sage: K55.ramified_primes()
            [5, 7, 11]
            sage: E.heegner_point(-7).ring_class_field().ramified_primes()
            [7]
        """
        return prime_divisors(self.__D * self.__c)

    def _repr_(self):
        """
        EXAMPLES::

            sage: heegner_point(37,-7,55).ring_class_field()._repr_()
            'Ring class field extension of QQ[sqrt(-7)] of conductor 55'
            sage: heegner_point(37,-7).ring_class_field()._repr_()
            'Hilbert class field of QQ[sqrt(-7)]'
        """
        c = self.__c
        if c == 1:
            return "Hilbert class field of QQ[sqrt(%s)]" % self.__D
        return "Ring class field extension of QQ[sqrt(%s)] of conductor %s" % (self.__D, self.__c)

    @cached_method
    def degree_over_K(self):
        """
        Return the relative degree of this ring class field over the
        quadratic imaginary field `K`.

        EXAMPLES::

            sage: E = EllipticCurve('389a'); P = E.heegner_point(-7,5)
            sage: K5 = P.ring_class_field(); K5
            Ring class field extension of QQ[sqrt(-7)] of conductor 5
            sage: K5.degree_over_K()
            6
            sage: type(K5.degree_over_K())
            <... 'sage.rings.integer.Integer'>

            sage: E = EllipticCurve('389a'); E.heegner_point(-20).ring_class_field().degree_over_K()
            2
            sage: E.heegner_point(-20,3).ring_class_field().degree_over_K()
            4
            sage: kronecker(-20,11)
            -1
            sage: E.heegner_point(-20,11).ring_class_field().degree_over_K()
            24
        """
        K = self.quadratic_field()

        # Multiply class number by relative degree of the Hilbert class field H over K.
        return K.class_number() * self.degree_over_H()

    @cached_method
    def degree_over_H(self):
        """
        Return the degree of this field over the Hilbert class field `H` of `K`.

        EXAMPLES::

            sage: E = EllipticCurve('389a')
            sage: E.heegner_point(-59).ring_class_field().degree_over_H()
            1
            sage: E.heegner_point(-59).ring_class_field().degree_over_K()
            3
            sage: QuadraticField(-59,'a').class_number()
            3

        Some examples in which prime dividing c is inert::

            sage: heegner_point(37,-7,3).ring_class_field().degree_over_H()
            4
            sage: heegner_point(37,-7,3^2).ring_class_field().degree_over_H()
            12
            sage: heegner_point(37,-7,3^3).ring_class_field().degree_over_H()
            36

        The prime dividing c is split.  For example, in the first case
        `O_K/cO_K` is isomorphic to a direct sum of two copies of
        ``GF(2)``, so the units are trivial::

            sage: heegner_point(37,-7,2).ring_class_field().degree_over_H()
            1
            sage: heegner_point(37,-7,4).ring_class_field().degree_over_H()
            2
            sage: heegner_point(37,-7,8).ring_class_field().degree_over_H()
            4

        Now c is ramified::

            sage: heegner_point(37,-7,7).ring_class_field().degree_over_H()
            7
            sage: heegner_point(37,-7,7^2).ring_class_field().degree_over_H()
            49

        Check that :issue:`15218` is solved::

            sage: E = EllipticCurve("19a");
            sage: s = E.heegner_point(-3,2).ring_class_field().galois_group().complex_conjugation()
            sage: H = s.domain(); H.absolute_degree()
            2
        """
        c = self.__c
        if c == 1:
            return ZZ.one()

        # Let K_c be the ring class field.  We have by class field theory that
        #           Gal(K_c / H) = (O_K / c O_K)^* / ((Z/cZ)^* M),
        # where M is the image of the roots of unity of K in (O_K / c O_K)^*.
        #
        # To compute the cardinality of the above Galois group, we
        # first reduce to the case that c = p^e is a prime power
        # (since the expression is multiplicative in c).
        # Of course, note also that #(Z/cZ)^* = phi(c)
        #
        # Case 1: p splits in O_K.  Then
        #         #(O_K/p^e*O_K)^* = (#(Z/p^eZ)^*)^2 = phi(p^e)^2, so
        #           #(O_K/p^e*O_K)^*/(Z/p^eZ)^* = phi(p^e) = p^e - p^(e-1)
        #
        # Case 2: p is inert in O_K.  Then
        #         #(O_K/p^e O_K)^* = p^(2*e)-p^(2*(e-1))
        #         so #(O_K/p^e*O_K)^*/(Z/p^eZ)^*
        #              = (p^(2*e)-p^(2*(e-1)))/(p^e-p^(e-1)) = p^e + p^(e-1).
        #
        # Case 3: p ramified in O_K. Then
        #         #(O_K/p^e O_K)^* = p^(2*e) - p^(2*e-1),
        #         so #(O_K/p^e O_K)^*/#(Z/p^eZ)^* = p^e.
        #
        # Section 4.2 of Cohen's "Advanced Computational Algebraic
        # Number Theory" GTM is also relevant, though Cohen is working
        # with *ray* class fields and here we want the cardinality
        # of the *ring* class field, which is a subfield.

        K = self.quadratic_field()

        n = ZZ(1)
        for p, e in c.factor():
            F = K.factor(p)
            if len(F) == 2:
                # split case
                n *= p**e - p**(e-1)
            else:
                if F[0][1] > 1:
                    # ramified case
                    n *= p**e
                else:
                    # inert case
                    n *= p**e + p**(e-1)
        return (n * ZZ(2)) // K.number_of_roots_of_unity()

    @cached_method
    def absolute_degree(self):
        r"""
        Return the absolute degree of this field over `\QQ`.

        EXAMPLES::

            sage: E = EllipticCurve('389a'); K = E.heegner_point(-7,5).ring_class_field()
            sage: K.absolute_degree()
            12
            sage: K.degree_over_K()
            6
        """
        return 2*self.degree_over_K()

    degree_over_Q = absolute_degree

    @cached_method
    def quadratic_field(self):
        r"""
        Return the quadratic imaginary field `K = \QQ(\sqrt{D})`.

        EXAMPLES::

            sage: E = EllipticCurve('389a'); K = E.heegner_point(-7,5).ring_class_field()
            sage: K.quadratic_field()
            Number Field in sqrt_minus_7 with defining polynomial x^2 + 7
             with sqrt_minus_7 = 2.645751311064591?*I
        """
        D = self.__D
        var = 'sqrt_minus_%s' % (-D)
        return QuadraticField(D,var)

    @cached_method
    def galois_group(self, base=QQ):
        r"""
        Return the Galois group of ``self`` over base.

        INPUT:

        - ``base`` -- (default: `\QQ`) a subfield of ``self`` or `\QQ`

        EXAMPLES::

            sage: E = EllipticCurve('389a')
            sage: A = E.heegner_point(-7,5).ring_class_field()
            sage: A.galois_group()
            Galois group of Ring class field extension of QQ[sqrt(-7)] of conductor 5
            sage: B = E.heegner_point(-7).ring_class_field()
            sage: C = E.heegner_point(-7,15).ring_class_field()
            sage: A.galois_group()
            Galois group of Ring class field extension of QQ[sqrt(-7)] of conductor 5
            sage: A.galois_group(B)
            Galois group of Ring class field extension of QQ[sqrt(-7)] of conductor 5
             over Hilbert class field of QQ[sqrt(-7)]
            sage: A.galois_group().cardinality()
            12
            sage: A.galois_group(B).cardinality()
            6
            sage: C.galois_group(A)
            Galois group of Ring class field extension of QQ[sqrt(-7)] of conductor 15
             over Ring class field extension of QQ[sqrt(-7)] of conductor 5
            sage: C.galois_group(A).cardinality()
            4
        """
        return GaloisGroup(self, base)

    def is_subfield(self, M):
        """
        Return ``True`` if this ring class field is a subfield of the ring class field `M`.
        If `M` is not a ring class field, then a :exc:`TypeError` is raised.

        EXAMPLES::

            sage: E = EllipticCurve('389a')
            sage: A = E.heegner_point(-7,5).ring_class_field()
            sage: B = E.heegner_point(-7).ring_class_field()
            sage: C = E.heegner_point(-20).ring_class_field()
            sage: D = E.heegner_point(-7,15).ring_class_field()
            sage: B.is_subfield(A)
            True
            sage: B.is_subfield(B)
            True
            sage: B.is_subfield(D)
            True
            sage: B.is_subfield(C)
            False
            sage: A.is_subfield(B)
            False
            sage: A.is_subfield(D)
            True
        """
        if not isinstance(M, RingClassField):
            raise TypeError("M must be a ring class field")
        return self.quadratic_field() == M.quadratic_field() and \
               M.conductor() % self.conductor() == 0

##################################################################################
#
# Galois groups of ring class fields
#
##################################################################################


class GaloisGroup(SageObject):
    """
    A Galois group of a ring class field.

    EXAMPLES::

        sage: E = EllipticCurve('389a')
        sage: G = E.heegner_point(-7,5).ring_class_field().galois_group(); G
        Galois group of Ring class field extension of QQ[sqrt(-7)] of conductor 5
        sage: G.field()
        Ring class field extension of QQ[sqrt(-7)] of conductor 5
        sage: G.cardinality()
        12
        sage: G.complex_conjugation()
        Complex conjugation automorphism of Ring class field extension of QQ[sqrt(-7)]
         of conductor 5

    TESTS::

        sage: G = heegner_point(37,-7).ring_class_field().galois_group()
        sage: loads(dumps(G)) == G
        True
        sage: type(G)
        <class 'sage.schemes.elliptic_curves.heegner.GaloisGroup'>
    """
    def __init__(self, field, base=QQ):
        r"""
        INPUT:

        - ``field`` -- a ring class field

        - ``base`` -- subfield of field (default: `\QQ`)

        EXAMPLES::

            sage: K5 = heegner_points(389,-7,5).ring_class_field()
            sage: K1 = heegner_points(389,-7,1).ring_class_field()
            sage: sage.schemes.elliptic_curves.heegner.GaloisGroup(K5,K1)
            Galois group of Ring class field extension of QQ[sqrt(-7)] of conductor 5
             over Hilbert class field of QQ[sqrt(-7)]
            sage: K5.galois_group(K1)
            Galois group of Ring class field extension of QQ[sqrt(-7)] of conductor 5
             over Hilbert class field of QQ[sqrt(-7)]
        """
        if not isinstance(field, RingClassField):
            raise TypeError("field must be of type RingClassField")
        if base != QQ and base != field.quadratic_field():
            if not isinstance(base, RingClassField):
                raise TypeError("base must be of type RingClassField or QQ or quadratic field")
            if not base.is_subfield(field):
                raise TypeError("base must be a subfield of field")
        self.__field = field
        self.__base = base

    def __eq__(self, G):
        """
        EXAMPLES::

            sage: G = EllipticCurve('389a').heegner_point(-7,5).ring_class_field().galois_group()
            sage: G == G
            True
            sage: G == 0
            False
            sage: H = EllipticCurve('389a').heegner_point(-7,11).ring_class_field().galois_group()
            sage: G == H
            False
        """
        return isinstance(G, GaloisGroup) and (G.__field,G.__base) == (self.__field,self.__base)

    def __ne__(self, other):
        """
        EXAMPLES::

            sage: G = EllipticCurve('389a').heegner_point(-7,5).ring_class_field().galois_group()
            sage: G != G
            False
            sage: G != 0
            True
            sage: H = EllipticCurve('389a').heegner_point(-7,11).ring_class_field().galois_group()
            sage: G != H
            True
        """
        return not (self == other)

    def __hash__(self):
        """
        Return hash of this Galois group, which is the same as the
        hash of the pair, the field and its base.

        EXAMPLES::

            sage: G = EllipticCurve('389a').heegner_point(-7,5).ring_class_field().galois_group()
            sage: hash(G) == hash((G.field(), G.base_field()))
            True
        """
        return hash((self.__field, self.__base))

    def __call__(self, x):
        """
        Coerce `x` into ``self``, where `x` is a Galois group element, or
        in case ``self`` has base field the Hilbert class field, `x` can
        also be an element of the ring of integers.

        INPUT:

        - ``x`` -- automorphism or quadratic field element

        OUTPUT: an automorphism (or :exc:`TypeError`)

        EXAMPLES::

            sage: K5 = heegner_points(389,-52,5).ring_class_field()
            sage: K1 = heegner_points(389,-52,1).ring_class_field()
            sage: G = K5.galois_group(K1)
            sage: G(1)
            Class field automorphism defined by x^2 + 325*y^2
            sage: G(G[0])
            Class field automorphism defined by x^2 + 325*y^2
            sage: alpha = 2 + K1.quadratic_field().gen(); alpha
            sqrt_minus_52 + 2
            sage: G(alpha)
            Class field automorphism defined by 14*x^2 - 10*x*y + 25*y^2

        A :exc:`TypeError` is raised when the coercion is not possible::

            sage: G(0)
            Traceback (most recent call last):
            ...
            TypeError: x does not define element of (O_K/c*O_K)^*
        """
        if isinstance(x, GaloisAutomorphism) and x.parent() == self:
            return x
        try:
            return self._alpha_to_automorphism(x)
        except (ZeroDivisionError, TypeError):
            raise TypeError("x does not define element of (O_K/c*O_K)^*")

    def _repr_(self):
        """
        Return string representation of this Galois group.

        EXAMPLES::

            sage: E = EllipticCurve('389a')
            sage: G = E.heegner_point(-7,5).ring_class_field().galois_group()
            sage: G._repr_()
            'Galois group of Ring class field extension of QQ[sqrt(-7)] of conductor 5'
        """
        if self.base_field() != QQ:
            s = " over %s" % self.base_field()
        else:
            s = ''
        return "Galois group of %s%s" % (self.field(), s)

    def field(self):
        """
        Return the ring class field that this Galois group acts on.

        EXAMPLES::

            sage: G = heegner_point(389,-7,5).ring_class_field().galois_group()
            sage: G.field()
            Ring class field extension of QQ[sqrt(-7)] of conductor 5
        """
        return self.__field

    def base_field(self):
        """
        Return the base field, which the field fixed by all the
        automorphisms in this Galois group.

        EXAMPLES::

            sage: x = heegner_point(37,-7,5)
            sage: Kc = x.ring_class_field(); Kc
            Ring class field extension of QQ[sqrt(-7)] of conductor 5
            sage: K = x.quadratic_field()
            sage: G = Kc.galois_group(); G
            Galois group of Ring class field extension of QQ[sqrt(-7)] of conductor 5
            sage: G.base_field()
            Rational Field
            sage: G.cardinality()
            12
            sage: Kc.absolute_degree()
            12
            sage: G = Kc.galois_group(K); G
            Galois group of Ring class field extension of QQ[sqrt(-7)] of conductor 5
             over Number Field in sqrt_minus_7 with defining polynomial x^2 + 7
              with sqrt_minus_7 = 2.645751311064591?*I
            sage: G.cardinality()
            6
            sage: G.base_field()
            Number Field in sqrt_minus_7 with defining polynomial x^2 + 7
             with sqrt_minus_7 = 2.645751311064591?*I
            sage: G = Kc.galois_group(Kc); G
            Galois group of Ring class field extension of QQ[sqrt(-7)] of conductor 5
             over Ring class field extension of QQ[sqrt(-7)] of conductor 5
            sage: G.cardinality()
            1
            sage: G.base_field()
            Ring class field extension of QQ[sqrt(-7)] of conductor 5
        """
        return self.__base

    @cached_method
    def kolyvagin_generators(self):
        r"""
        Assuming this Galois group `G` is of the form
        `G=\textrm{Gal}(K_c/K_1)`, with `c=p_1\dots p_n` satisfying the
        Kolyvagin hypothesis, this function returns noncanonical
        choices of lifts of generators for each of the cyclic factors
        of `G` corresponding to the primes dividing `c`.  Thus the
        `i`-th returned valued is an element of `G` that maps to the
        identity element of `\textrm{Gal}(K_p/K_1)` for all `p \neq p_i` and
        to a choice of generator of `\textrm{Gal}(K_{p_i}/K_1)`.

        OUTPUT: list of elements of ``self``

        EXAMPLES::

            sage: K3 = heegner_points(389,-52,3).ring_class_field()
            sage: K1 = heegner_points(389,-52,1).ring_class_field()
            sage: G = K3.galois_group(K1)
            sage: G.kolyvagin_generators()
            (Class field automorphism defined by 9*x^2 - 6*x*y + 14*y^2,)

            sage: K5 = heegner_points(389,-52,5).ring_class_field()
            sage: K1 = heegner_points(389,-52,1).ring_class_field()
            sage: G = K5.galois_group(K1)
            sage: G.kolyvagin_generators()
            (Class field automorphism defined by 17*x^2 - 14*x*y + 22*y^2,)
        """
        M = self.field()
        c = M.conductor()
        if not (self._base_is_hilbert_class_field() and self.is_kolyvagin()):
            raise ValueError("field must be of the form Gal(K_c/K_1)")
        if not c.is_prime():
            raise NotImplementedError("only implemented when c is prime")

        # Since c satisfies Kolyvagin and is prime, the group is cyclic,
        # so we just find a generator.
        for sigma in self:
            if sigma.order() == self.cardinality():
                return (sigma,)

        raise NotImplementedError

    @cached_method
    def lift_of_hilbert_class_field_galois_group(self):
        r"""
        Assuming this Galois group `G` is of the form `G=\textrm{Gal}(K_c/K)`,
        this function returns noncanonical choices of lifts of the
        elements of the quotient group `\textrm{Gal}(K_1/K)`.

        OUTPUT: tuple of elements of ``self``

        EXAMPLES::

            sage: K5 = heegner_points(389,-52,5).ring_class_field()
            sage: G = K5.galois_group(K5.quadratic_field())
            sage: G.lift_of_hilbert_class_field_galois_group()
            (Class field automorphism defined by x^2 + 325*y^2,
             Class field automorphism defined by 2*x^2 + 2*x*y + 163*y^2)
            sage: G.cardinality()
            12
            sage: K5.quadratic_field().class_number()
            2
        """
        if not self._base_is_quad_imag_field():
            raise ValueError("Galois group must be of the form Gal(K_c/K)")
        K = self.base_field()
        C = K.class_group()
        v = []
        lifts = []
        for sigma in self:
            I = sigma.ideal()
            g = C(I)
            if g not in v:
                v.append(g)
                lifts.append(sigma)
        return tuple(lifts)

    @cached_method
    def _list(self):
        r"""
        Enumerate the elements of ``self``.

        EXAMPLES:

        Example with order 1 (a special case)::

            sage: E = EllipticCurve('389a'); F = E.heegner_point(-7,1).ring_class_field()
            sage: G = F.galois_group(F.quadratic_field())
            sage: G._list()
            (Class field automorphism defined by x^2 + x*y + 2*y^2,)

        Example over quadratic imaginary field::

            sage: E = EllipticCurve('389a'); F = E.heegner_point(-7,5).ring_class_field()
            sage: G = F.galois_group(F.quadratic_field())
            sage: G._list()
            (Class field automorphism defined by x^2 + x*y + 44*y^2,
             Class field automorphism defined by 2*x^2 - x*y + 22*y^2,
             Class field automorphism defined by 2*x^2 + x*y + 22*y^2,
             Class field automorphism defined by 4*x^2 - x*y + 11*y^2,
             Class field automorphism defined by 4*x^2 + x*y + 11*y^2,
             Class field automorphism defined by 7*x^2 + 7*x*y + 8*y^2)

        Example over `\QQ` (it is not implemented yet)::

            sage: K3 = heegner_points(389,-52,3).ring_class_field()
            sage: K3.galois_group()._list()
            Traceback (most recent call last):
            ...
            NotImplementedError: Galois group over QQ not yet implemented

        Example over Hilbert class field::

            sage: K3 = heegner_points(389,-52,3).ring_class_field()
            sage: K1 = heegner_points(389,-52,1).ring_class_field()
            sage: G = K3.galois_group(K1)
            sage: G._list()
            (Class field automorphism defined by x^2 + 117*y^2, Class field automorphism defined by 9*x^2 - 6*x*y + 14*y^2,
             Class field automorphism defined by 9*x^2 + 13*y^2, Class field automorphism defined by 9*x^2 + 6*x*y + 14*y^2)
        """
        if self._base_is_QQ():
            raise NotImplementedError("Galois group over QQ not yet implemented")
        elif self._base_is_quad_imag_field():
            # Over the quadratic imaginary field, so straightforward
            # enumeration of all reduced primitive binary quadratic
            # forms of discriminant D*c^2.
            D = self.base_field().discriminant()
            c = self.field().conductor()
            Q = [f for f in BinaryQF_reduced_representatives(D*c*c) if f.is_primitive()]
            v = [GaloisAutomorphismQuadraticForm(self, f) for f in Q]

        elif self._base_is_hilbert_class_field() and self.is_kolyvagin():
            # Take only the automorphisms in the quad imag case that map to
            # a principal ideal.
            M = self.field()
            K = M.quadratic_field()
            v = []
            self.__p1_to_automorphism = {}
            for sigma in M.galois_group(K)._list():
                I = sigma.ideal()
                if I.is_principal():
                    # sigma does define an element of our Galois subgroup.
                    alpha = sigma.ideal().gens_reduced()[0]
                    t = GaloisAutomorphismQuadraticForm(self, sigma.quadratic_form(), alpha=alpha)
                    self.__p1_to_automorphism[t.p1_element()] = t
                    v.append(t)
        else:
            raise NotImplementedError("general Galois group not yet implemented")

        v.sort()
        assert len(v) == self.cardinality(), "bug enumerating Galois group elements"
        return tuple(v)

    def _quadratic_form_to_alpha(self, f):
        """
        INPUT:

        - ``f`` -- a binary quadratic form with discriminant `c^2 D`

        OUTPUT: an element of the ring of integers of the quadratic imaginary field

        EXAMPLES::

            sage: K3 = heegner_points(389,-52,3).ring_class_field(); K1 = heegner_points(389,-52,1).ring_class_field()
            sage: G = K3.galois_group(K1)
            sage: [G._quadratic_form_to_alpha(s.quadratic_form()) for s in G]
            [3/2*sqrt_minus_52, 1/6*sqrt_minus_52 + 1/3, 1/6*sqrt_minus_52, 1/6*sqrt_minus_52 - 1/3]

        What happens when we input a quadratic form that has nothing
        to do with `G`::

            sage: G._quadratic_form_to_alpha(BinaryQF([1,2,3]))
            Traceback (most recent call last):
            ...
            ValueError: quadratic form has the wrong discriminant
        """
        A,B,C = f
        K = self.field().quadratic_field()
        if f.discriminant() != self.field().conductor()**2 * K.discriminant():
            raise ValueError("quadratic form has the wrong discriminant")

        R = K['X']
        v = R([C,B,A]).roots()[0][0]
        return v

    def _alpha_to_automorphism(self, alpha):
        r"""
        Assuming ``self`` has base field the Hilbert class field, make an
        automorphism from the element `\alpha` of the ring of integers
        into ``self``.

        INPUT:

        - ``alpha`` -- element of quadratic imaginary field coprime to conductor

        EXAMPLES::

            sage: K3 = heegner_points(389,-52,3).ring_class_field()
            sage: K1 = heegner_points(389,-52,1).ring_class_field()
            sage: G = K3.galois_group(K1)
            sage: G._alpha_to_automorphism(1)
            Class field automorphism defined by x^2 + 117*y^2
            sage: [G._alpha_to_automorphism(s.alpha()) for s in G] == list(G)
            True
        """
        if not self._base_is_hilbert_class_field() and self.is_kolyvagin():
            raise TypeError("base must be Hilbert class field with Kolyvagin condition on conductor")
        R = self.field().quadratic_field().maximal_order()
        uv = self._alpha_to_p1_element(R(alpha))
        try:
            d = self.__p1_to_automorphism
        except AttributeError:
            self._list()  # computes attribute as side-effect
            d = self.__p1_to_automorphism
        return d[uv]

    def _alpha_to_p1_element(self, alpha):
        r"""
        Given an element of the ring of integers that is nonzero
        modulo c, return canonical (after our fixed choice of basis)
        element of the project line corresponding to it.

        INPUT:

        - ``alpha`` -- element of the ring of integers of the quadratic
          imaginary field

        OUTPUT: 2-tuple of integers

        EXAMPLES::

            sage: K3 = heegner_points(389,-52,3).ring_class_field()
            sage: K1 = heegner_points(389,-52,1).ring_class_field()
            sage: G = K3.galois_group(K1)
            sage: G._alpha_to_p1_element(1)
            (1, 0)
            sage: sorted([G._alpha_to_p1_element(s.alpha()) for s in G])
            [(0, 1), (1, 0), (1, 1), (1, 2)]
        """
        try:
            A, P1 = self.__alpha_to_p1_element
        except AttributeError:
            # todo (optimize) -- this whole function can be massively optimized:
            M = self.field()
            A = M.quadratic_field().maximal_order().free_module()
            P1 = P1List(M.conductor())
            self.__alpha_to_p1_element = A, P1
        alpha = self.field().quadratic_field()(alpha)
        w = A.coordinate_vector(alpha.vector())
        w *= w.denominator()
        w = w.change_ring(ZZ)
        n = gcd(w)
        w /= n
        c = P1.N()
        w = P1.normalize(ZZ(w[0]) % c, ZZ(w[1]) % c)
        if w == (0,0):
            w = (1,0)
        return w

    def _p1_element_to_alpha(self, uv):
        """
        Convert a normalized pair ``uv=(u,v)`` of integers to the
        corresponding element of the ring of integers got by taking `u
        b_0 + v b_1` where `b_0, b_1` are the basis for the ring of
        integers.

        INPUT:

        - ``uv`` -- pair of integers

        OUTPUT: element of maximal order of quadratic field

        EXAMPLES::

            sage: K5 = heegner_points(389,-52,5).ring_class_field()
            sage: K1 = heegner_points(389,-52,1).ring_class_field()
            sage: G = K5.galois_group(K1)
            sage: v = [G._alpha_to_p1_element(s.alpha()) for s in G]
            sage: [G._p1_element_to_alpha(z) for z in v]
            [1, 1/2*sqrt_minus_52, 1/2*sqrt_minus_52 + 1, 2*sqrt_minus_52 + 1, sqrt_minus_52 + 1, 3/2*sqrt_minus_52 + 1]
            sage: [G(G._p1_element_to_alpha(z)) for z in v] == list(G)
            True
        """
        B = self.field().quadratic_field().maximal_order().basis()
        return uv[0]*B[0] + uv[1]*B[1]

    def _base_is_QQ(self):
        r"""
        Return ``True`` if the base field of this ring class field is `\QQ`.

        EXAMPLES::

            sage: H = heegner_points(389,-20,3); M = H.ring_class_field()
            sage: M.galois_group(H.quadratic_field())._base_is_QQ()
            False
            sage: M.galois_group(QQ)._base_is_QQ()
            True
            sage: M.galois_group(heegner_points(389,-20,1).ring_class_field())._base_is_QQ()
            False
        """
        return self.__base == QQ

    def _base_is_quad_imag_field(self):
        """
        Return ``True`` if the base field of this ring class field is the
        quadratic imaginary field `K`.

        EXAMPLES::

            sage: H = heegner_points(389,-20,3); M = H.ring_class_field()
            sage: M.galois_group(H.quadratic_field())._base_is_quad_imag_field()
            True
            sage: M.galois_group(QQ)._base_is_quad_imag_field()
            False
            sage: M.galois_group(heegner_points(389,-20,1).ring_class_field())._base_is_quad_imag_field()
            False
        """
        return isinstance(self.__base, sage.rings.abc.NumberField_quadratic)

    def is_kolyvagin(self):
        """
        Return ``True`` if conductor `c` is prime to the discriminant of the
        quadratic field, `c` is squarefree and each prime dividing `c`
        is inert.

        EXAMPLES::

            sage: K5 = heegner_points(389,-52,5).ring_class_field()
            sage: K1 = heegner_points(389,-52,1).ring_class_field()
            sage: K5.galois_group(K1).is_kolyvagin()
            True
            sage: K7 = heegner_points(389,-52,7).ring_class_field()
            sage: K7.galois_group(K1).is_kolyvagin()
            False
            sage: K25 = heegner_points(389,-52,25).ring_class_field()
            sage: K25.galois_group(K1).is_kolyvagin()
            False
        """
        M = self.field()
        c = M.conductor()
        D = M.quadratic_field().discriminant()
        if c.gcd(D) != 1:
            return False
        if not c.is_squarefree():
            return False
        for p in c.prime_divisors():
            if not is_inert(D,p):
                return False
        return True

    def _base_is_hilbert_class_field(self):
        """
        Return ``True`` if the base field of this ring class field is the
        Hilbert class field of `K` viewed as a ring class field (so
        not of data type QuadraticField).

        EXAMPLES::

            sage: H = heegner_points(389,-20,3); M = H.ring_class_field()
            sage: M.galois_group(H.quadratic_field())._base_is_hilbert_class_field()
            False
            sage: M.galois_group(QQ)._base_is_hilbert_class_field()
            False
            sage: M.galois_group(heegner_points(389,-20,1).ring_class_field())._base_is_hilbert_class_field()
            True
        """
        M = self.__base
        return isinstance(M, RingClassField) and M.conductor() == 1

    def __getitem__(self, i):
        """
        EXAMPLES::

            sage: E = EllipticCurve('389a'); F = E.heegner_point(-7,5).ring_class_field()
            sage: G = F.galois_group(F.quadratic_field())
            sage: G[0]
            Class field automorphism defined by x^2 + x*y + 44*y^2
        """
        return self._list()[i]

    def __len__(self):
        """
        EXAMPLES::

            sage: K5 = heegner_points(389,-52,5).ring_class_field()
            sage: K1 = heegner_points(389,-52,1).ring_class_field()
            sage: G = K5.galois_group(K1)
            sage: G.cardinality()
            6
            sage: len(G)
            6
        """
        return self.cardinality()

    @cached_method
    def cardinality(self):
        """
        Return the cardinality of this Galois group.

        EXAMPLES::

            sage: E = EllipticCurve('389a')
            sage: G = E.heegner_point(-7,5).ring_class_field().galois_group(); G
            Galois group of Ring class field extension of QQ[sqrt(-7)] of conductor 5
            sage: G.cardinality()
            12
            sage: G = E.heegner_point(-7).ring_class_field().galois_group()
            sage: G.cardinality()
            2
            sage: G = E.heegner_point(-7,55).ring_class_field().galois_group()
            sage: G.cardinality()
            120
        """
        return self.__field.absolute_degree() // self.__base.absolute_degree()

    @cached_method
    def complex_conjugation(self):
        """
        Return the automorphism of ``self`` determined by complex
        conjugation.  The base field must be the rational numbers.

        EXAMPLES::

            sage: E = EllipticCurve('389a')
            sage: G = E.heegner_point(-7,5).ring_class_field().galois_group()
            sage: G.complex_conjugation()
            Complex conjugation automorphism of Ring class field extension
             of QQ[sqrt(-7)] of conductor 5
        """
        if self.base_field() != QQ:
            raise ValueError("the base field must be fixed by complex conjugation")
        return GaloisAutomorphismComplexConjugation(self)


##################################################################################
#
# Elements of Galois groups
#
##################################################################################


class GaloisAutomorphism(SageObject):
    """
    An abstract automorphism of a ring class field.

    .. TODO::

        make :class:`GaloisAutomorphism` derive from GroupElement, so
        that one gets powers for free, etc.
    """
    def __init__(self, parent):
        """
        INPUT:

        - ``parent`` -- a group of automorphisms of a ring class field

        EXAMPLES::

            sage: G = heegner_points(389,-7,5).ring_class_field().galois_group(); G
            Galois group of Ring class field extension of QQ[sqrt(-7)] of conductor 5
            sage: sage.schemes.elliptic_curves.heegner.GaloisAutomorphism(G)
            <sage.schemes.elliptic_curves.heegner.GaloisAutomorphism object at ...>
        """
        self.__parent = parent

    def parent(self):
        """
        Return the parent of this automorphism, which is a Galois
        group of a ring class field.

        EXAMPLES::

            sage: E = EllipticCurve('389a')
            sage: G = E.heegner_point(-7,5).ring_class_field().galois_group()
            sage: s = G.complex_conjugation()
            sage: s.parent()
            Galois group of Ring class field extension of QQ[sqrt(-7)] of conductor 5
        """
        return self.__parent

    def domain(self):
        """
        Return the domain of this automorphism.

        EXAMPLES::

            sage: E = EllipticCurve('389a')
            sage: G = E.heegner_point(-7,5).ring_class_field().galois_group()
            sage: s = G.complex_conjugation()
            sage: s.domain()
            Ring class field extension of QQ[sqrt(-7)] of conductor 5
        """
        return self.parent().field()


class GaloisAutomorphismComplexConjugation(GaloisAutomorphism):
    """
    The complex conjugation automorphism of a ring class field.

    EXAMPLES::

        sage: G = heegner_point(37,-7,5).ring_class_field().galois_group()
        sage: conj = G.complex_conjugation()
        sage: conj
        Complex conjugation automorphism of Ring class field extension
         of QQ[sqrt(-7)] of conductor 5
        sage: conj.domain()
        Ring class field extension of QQ[sqrt(-7)] of conductor 5

    TESTS::

        sage: type(conj)
        <class 'sage.schemes.elliptic_curves.heegner.GaloisAutomorphismComplexConjugation'>
        sage: loads(dumps(conj)) == conj
        True
    """
    def __init__(self, parent):
        """
        INPUT:

            - ``parent`` -- a group of automorphisms of a ring class field

        EXAMPLES::

            sage: G = heegner_point(37,-7,5).ring_class_field().galois_group()
            sage: sage.schemes.elliptic_curves.heegner.GaloisAutomorphismComplexConjugation(G)
            Complex conjugation automorphism of Ring class field extension
             of QQ[sqrt(-7)] of conductor 5
        """
        GaloisAutomorphism.__init__(self, parent)

    def __hash__(self):
        """
        The hash value is the same as the hash value of the
        pair ``(self.parent(), 1)``.

        EXAMPLES::

            sage: G = EllipticCurve('389a').heegner_point(-7,5).ring_class_field().galois_group()
            sage: conj = G.complex_conjugation()
            sage: hash(conj) == hash((conj.parent(), 1))
            True
        """
        return hash((self.parent(), 1))

    def __eq__(self, right):
        """
        EXAMPLES::

            sage: G = EllipticCurve('389a').heegner_point(-7,5).ring_class_field().galois_group()
            sage: conj = G.complex_conjugation()
            sage: conj2 = sage.schemes.elliptic_curves.heegner.GaloisAutomorphismComplexConjugation(G)
            sage: conj is conj2
            False
            sage: conj == conj2
            True
        """
        return isinstance(right, GaloisAutomorphismComplexConjugation) and \
               self.parent() == right.parent()

    def __ne__(self, other):
        """
        EXAMPLES::

            sage: G = EllipticCurve('389a').heegner_point(-7,5).ring_class_field().galois_group()
            sage: conj = G.complex_conjugation()
            sage: conj2 = sage.schemes.elliptic_curves.heegner.GaloisAutomorphismComplexConjugation(G)
            sage: conj != conj2
            False
        """
        return not (self == other)

    def _repr_(self):
        """
        Return print representation of the complex conjugation automorphism.

        EXAMPLES::

            sage: conj = heegner_point(37,-7,5).ring_class_field().galois_group().complex_conjugation()
            sage: conj._repr_()
            'Complex conjugation automorphism of Ring class field extension of QQ[sqrt(-7)] of conductor 5'
        """
        return "Complex conjugation automorphism of %s" % self.domain()

#     def __mul__(self, right):
#         """
#         Return the composition of two automorphisms.

#         EXAMPLES::

#             sage: ?
#         """
#         if self.parent() != right.__parent():
#             raise TypeError("automorphisms must be of the same class field")
#         raise NotImplementedError

    def __invert__(self):
        """
        Return the inverse of ``self``, which is just ``self`` again.

        EXAMPLES::

            sage: conj = heegner_point(37,-7,5).ring_class_field().galois_group().complex_conjugation()
            sage: ~conj
            Complex conjugation automorphism of Ring class field extension of QQ[sqrt(-7)] of conductor 5
        """
        return self

    def order(self):
        """
        EXAMPLES::

            sage: G = heegner_point(37,-7,5).ring_class_field().galois_group()
            sage: conj = G.complex_conjugation()
            sage: conj.order()
            2
        """
        return ZZ(2)


@richcmp_method
class GaloisAutomorphismQuadraticForm(GaloisAutomorphism):
    """
    An automorphism of a ring class field defined by a quadratic form.

    EXAMPLES::

        sage: H = heegner_points(389,-20,3)
        sage: sigma = H.ring_class_field().galois_group(H.quadratic_field())[0]; sigma
        Class field automorphism defined by x^2 + 45*y^2
        sage: type(sigma)
        <class 'sage.schemes.elliptic_curves.heegner.GaloisAutomorphismQuadraticForm'>
        sage: loads(dumps(sigma)) == sigma
        True
    """
    def __init__(self, parent, quadratic_form, alpha=None):
        r"""
        INPUT:

        - ``parent`` -- a group of automorphisms of a ring class field

        - ``quadratic_form`` -- a binary quadratic form that
          defines an element of the Galois group of `K_c` over `K`

        - ``alpha`` -- (default: ``None``) optional data that specified
          element corresponding element of `(\mathcal{O}_K /
          c\mathcal{O}_K)^* / (\ZZ/c\ZZ)^*`, via class field theory.

        EXAMPLES::

            sage: H = heegner_points(389,-20,3)
            sage: G = H.ring_class_field().galois_group(H.quadratic_field())
            sage: f = BinaryQF_reduced_representatives(-20*9)[0]
            sage: sage.schemes.elliptic_curves.heegner.GaloisAutomorphismQuadraticForm(G, f)
            Class field automorphism defined by x^2 + 45*y^2
        """
        self.__quadratic_form = quadratic_form.reduced_form()
        self.__alpha = alpha
        GaloisAutomorphism.__init__(self, parent)

    @cached_method
    def order(self):
        """
        Return the multiplicative order of this Galois group automorphism.

        EXAMPLES::

            sage: K3 = heegner_points(389,-52,3).ring_class_field()
            sage: K1 = heegner_points(389,-52,1).ring_class_field()
            sage: G = K3.galois_group(K1)
            sage: sorted([g.order() for g in G])
            [1, 2, 4, 4]
            sage: K5 = heegner_points(389,-52,5).ring_class_field()
            sage: K1 = heegner_points(389,-52,1).ring_class_field()
            sage: G = K5.galois_group(K1)
            sage: sorted([g.order() for g in G])
            [1, 2, 3, 3, 6, 6]
        """
        alpha = self.__alpha
        if alpha is None:
            raise NotImplementedError("order only currently implemented when alpha given in construction")
        G = self.parent()
        one = G(1).p1_element()
        ans = ZZ(1)
        z = alpha
        for i in range(G.cardinality()):
            if G._alpha_to_p1_element(z) == one:
                return ans
            ans += 1
            z *= alpha
        assert False, "bug in order"

    def alpha(self):
        r"""
        Optional data that specified element corresponding element of
        `(\mathcal{O}_K / c\mathcal{O}_K)^* / (\ZZ/c\ZZ)^*`, via class
        field theory.

        This is a generator of the ideal corresponding to this
        automorphism.

        EXAMPLES::

            sage: K3 = heegner_points(389,-52,3).ring_class_field()
            sage: K1 = heegner_points(389,-52,1).ring_class_field()
            sage: G = K3.galois_group(K1)
            sage: orb = sorted([g.alpha() for g in G]); orb # random (the sign depends on the database being installed or not)
            [1, 1/2*sqrt_minus_52 + 1, -1/2*sqrt_minus_52, 1/2*sqrt_minus_52 - 1]
            sage: sorted([x^2 for x in orb]) # this is just for testing
            [-13, -sqrt_minus_52 - 12, sqrt_minus_52 - 12, 1]

            sage: K5 = heegner_points(389,-52,5).ring_class_field()
            sage: K1 = heegner_points(389,-52,1).ring_class_field()
            sage: G = K5.galois_group(K1)
            sage: orb = sorted([g.alpha() for g in G]); orb # random (the sign depends on the database being installed or not)
            [1, -1/2*sqrt_minus_52, 1/2*sqrt_minus_52 + 1, 1/2*sqrt_minus_52 - 1,
             1/2*sqrt_minus_52 - 2, -1/2*sqrt_minus_52 - 2]
            sage: sorted([x^2 for x in orb]) # just for testing
            [-13, -sqrt_minus_52 - 12, sqrt_minus_52 - 12,
             -2*sqrt_minus_52 - 9, 2*sqrt_minus_52 - 9, 1]
        """
        if self.__alpha is None:
            raise ValueError("alpha data not defined")
        return self.__alpha

    @cached_method
    def p1_element(self):
        r"""
        Return element of the projective line corresponding to this
        automorphism.

        This only makes sense if this automorphism is in the Galois
        group `\textrm{Gal}(K_c/K_1)`.

        EXAMPLES::

            sage: K3 = heegner_points(389,-52,3).ring_class_field()
            sage: K1 = heegner_points(389,-52,1).ring_class_field()
            sage: G = K3.galois_group(K1)
            sage: sorted([g.p1_element() for g in G])
            [(0, 1), (1, 0), (1, 1), (1, 2)]

            sage: K5 = heegner_points(389,-52,5).ring_class_field()
            sage: K1 = heegner_points(389,-52,1).ring_class_field()
            sage: G = K5.galois_group(K1)
            sage: sorted([g.p1_element() for g in G])
            [(0, 1), (1, 0), (1, 1), (1, 2), (1, 3), (1, 4)]
        """
        return self.parent()._alpha_to_p1_element(self.__alpha)

    def __hash__(self):
        """
        The hash value is the hash of the pair formed by the parent
        and the quadratic form read as tuple.

        EXAMPLES::

            sage: H = heegner_points(389,-20,3)
            sage: s = H.ring_class_field().galois_group(H.quadratic_field())[0]
            sage: hash(s) == hash((s.parent(), tuple(s.quadratic_form())))
            True
        """
        return hash((self.parent(), tuple(self.__quadratic_form)))

    def __richcmp__(self, right, op):
        """
        Comparison.

        EXAMPLES::

            sage: H = heegner_points(389,-7,5)
            sage: s = H.ring_class_field().galois_group(H.quadratic_field())[1]
            sage: s == s
            True
            sage: s == s*s
            False
            sage: s == s*s*s*s*s
            False
            sage: s == s*s*s*s*s*s*s
            True

            sage: H = heegner_points(389,-20,3)
            sage: s = H.ring_class_field().galois_group(H.quadratic_field())[0]
            sage: s == s
            True
            sage: s == 0
            False
        """
        if not isinstance(right, GaloisAutomorphismQuadraticForm):
            return NotImplemented
        lx = self.parent()
        rx = right.parent()
        if lx != rx:
            return richcmp_not_equal(lx, rx, op)
        if self.quadratic_form().is_equivalent(right.quadratic_form()):
            return rich_to_bool(op, 0)
        return richcmp(self.quadratic_form(), right.quadratic_form(), op)

    def _repr_(self):
        """
        Return string representation of this automorphism.

        EXAMPLES::

            sage: H = heegner_points(389,-20,3); s = H.ring_class_field().galois_group(H.quadratic_field())[0]
            sage: s._repr_()
            'Class field automorphism defined by x^2 + 45*y^2'
        """
        return "Class field automorphism defined by %s" % self.__quadratic_form

    def __mul__(self, right):
        """
        Return the composition of two automorphisms.

        EXAMPLES::

            sage: H = heegner_points(389,-20,3); s = H.ring_class_field().galois_group(H.quadratic_field())[0]
            sage: s * s
            Class field automorphism defined by x^2 + 45*y^2
            sage: G = s.parent(); list(G)
            [Class field automorphism defined by x^2 + 45*y^2, Class field automorphism defined by 2*x^2 + 2*x*y + 23*y^2,
             Class field automorphism defined by 5*x^2 + 9*y^2, Class field automorphism defined by 7*x^2 + 4*x*y + 7*y^2]
            sage: G[0]*G[0]
            Class field automorphism defined by x^2 + 45*y^2
            sage: G[1]*G[2] == G[3]
            True
        """
        if self.parent() != right.parent():
            raise TypeError("automorphisms must be of the same class field")
        if not isinstance(right, GaloisAutomorphismQuadraticForm):
            # TODO: special case when right is complex conjugation
            raise NotImplementedError
        Q = (self.__quadratic_form * right.__quadratic_form).reduced_form()
        if self.__alpha and right.__alpha:
            alpha = self.__alpha * right.__alpha
        else:
            alpha = None
        return GaloisAutomorphismQuadraticForm(self.parent(), Q, alpha=alpha)

    def quadratic_form(self):
        """
        Return reduced quadratic form corresponding to this Galois
        automorphism.

        EXAMPLES::

            sage: H = heegner_points(389,-20,3); s = H.ring_class_field().galois_group(H.quadratic_field())[0]
            sage: s.quadratic_form()
            x^2 + 45*y^2
        """
        return self.__quadratic_form

    @cached_method
    def ideal(self):
        r"""
        Return ideal of ring of integers of quadratic imaginary field
        corresponding to this quadratic form.  This is the ideal

         `I = \left(A, \frac{-B+ c\sqrt{D}}{2}\right) \mathcal{O}_K`.

        EXAMPLES::

            sage: E = EllipticCurve('389a'); F = E.heegner_point(-20,3).ring_class_field()
            sage: G = F.galois_group(F.quadratic_field())
            sage: G[1].ideal()
            Fractional ideal (2, 1/2*sqrt_minus_20 + 1)
            sage: [s.ideal().gens() for s in G]
            [(1, 3/2*sqrt_minus_20), (2, 3/2*sqrt_minus_20 - 1),
             (5, 3/2*sqrt_minus_20), (7, 3/2*sqrt_minus_20 - 2)]
        """
        M = self.parent().field()
        K = M.quadratic_field()
        f = self.quadratic_form()
        c = M.conductor()
        sqrtD = K.gen()
        A, B, C = f
        if A % c == 0:
            A, C = C, A
        return K.fractional_ideal([A, (-B+c*sqrtD)/2])

##     def __call__(self, z):
##         """
##         Return image of the Heegner point `z` under this automorphism.
##
##         INPUT:
##
##             - ``z`` -- a Heegner point on `X_0(N)` or an elliptic curve
##
##         OUTPUT:
##
##             - a Heegner point
##
##         EXAMPLES::
##
##             sage: x = heegner_point(389,-20,3); F = x.ring_class_field()
##             sage: sigma = F.galois_group(F.quadratic_field())[1]; sigma
##             Class field automorphism defined by 2*x^2 + 2*x*y + 23*y^2
##             sage: sigma(x)
##             Heegner point 3/1556*sqrt(-20) - 495/778 of discriminant -20 and conductor 3 on X_0(389)
##         """
##         if isinstance(z, HeegnerPointOnX0N):
##             if z.ring_class_field() != self.domain():
##                 raise NotImplementedError("class fields must be the same")
##             # TODO -- check more compatibilities?
##             # TODO -- this is surely backwards -- something must be inverted?
##             f = z.quadratic_form() * self.quadratic_form()
##             # TODO -- put f into the correct form with A divisible by N, etc.?
##             # That could be done by looking up reduced form of f in a canonical
##             # list of best reps.
##             N,D,c = z.level(),z.discriminant(),z.conductor()
##             return HeegnerPointOnX0N(N,D,c, f = f)
##         else:
##             raise NotImplementedError

##################################################################################
#
# Specific Heegner points
#
##################################################################################


@richcmp_method
class HeegnerPoint(SageObject):
    r"""
    A Heegner point of level `N`, discriminant `D` and conductor `c`
    is any point on a modular curve or elliptic curve that is
    concocted in some way from a quadratic imaginary `\tau` in the upper
    half plane with `\Delta(\tau) = D c = \Delta(N \tau)`.

    EXAMPLES::

        sage: x = sage.schemes.elliptic_curves.heegner.HeegnerPoint(389,-7,13); x
        Heegner point of level 389, discriminant -7, and conductor 13
        sage: type(x)
        <class 'sage.schemes.elliptic_curves.heegner.HeegnerPoint'>
        sage: loads(dumps(x)) == x
        True
    """
    def __init__(self, N, D, c):
        """
        INPUT:

        - ``N`` -- (positive integer) the level

        - ``D`` -- (negative integer) fundamental discriminant

        - ``c`` -- (positive integer) conductor

        Since this is an abstract base class, no type or compatibility
        checks are done, as those are all assumed to be done in the
        derived class.

        EXAMPLES::

            sage: H = sage.schemes.elliptic_curves.heegner.HeegnerPoint(389,-7,5)
            sage: type(H)
            <class 'sage.schemes.elliptic_curves.heegner.HeegnerPoint'>
        """
        self.__N = N
        self.__D = D
        self.__c = c

    def __richcmp__(self, x, op):
        """
        Compare two Heegner points.

        EXAMPLES::

            sage: H = sage.schemes.elliptic_curves.heegner.HeegnerPoint(389,-7,5)
            sage: H == H
            True

            sage: H = sage.schemes.elliptic_curves.heegner.HeegnerPoint(389,-7,5); type(H)
            <class 'sage.schemes.elliptic_curves.heegner.HeegnerPoint'>
            sage: J = sage.schemes.elliptic_curves.heegner.HeegnerPoint(389,-7,11)
            sage: H == H
            True
            sage: H == J
            False
            sage: J == H
            False
            sage: H == 0
            False
        """
        if not isinstance(x, HeegnerPoint):
            return NotImplemented
        return richcmp((self.__N, self.__D, self.__c),
                       (x.__N, x.__D, x.__c), op)

    def _repr_(self):
        """
        EXAMPLES::

            sage: H = sage.schemes.elliptic_curves.heegner.HeegnerPoint(389,-7,5)
            sage: type(H)
            <class 'sage.schemes.elliptic_curves.heegner.HeegnerPoint'>
            sage: H._repr_()
            'Heegner point of level 389, discriminant -7, and conductor 5'
        """
        return "Heegner point of level %s, discriminant %s, and conductor %s" % (
            self.__N, self.__D, self.__c)

    def __hash__(self):
        """
        The hash value is obtained from level, discriminant, and conductor.

        EXAMPLES::

            sage: H = sage.schemes.elliptic_curves.heegner.HeegnerPoint(389,-7,5); type(H)
            <class 'sage.schemes.elliptic_curves.heegner.HeegnerPoint'>
            sage: hash(H)  == hash((H.level(), H.discriminant(), H.conductor()))
            True
        """
        return hash((self.__N, self.__D, self.__c))

    def level(self):
        """
        Return the level of this Heegner point, which is the level of the
        modular curve `X_0(N)` on which this is a Heegner point.

        EXAMPLES::

            sage: heegner_point(389,-7,5).level()
            389
        """
        return self.__N

    def conductor(self):
        """
        Return the conductor of this Heegner point.

        EXAMPLES::

            sage: heegner_point(389,-7,5).conductor()
            5
            sage: E = EllipticCurve('37a1'); P = E.kolyvagin_point(-67,7); P
            Kolyvagin point of discriminant -67 and conductor 7
             on elliptic curve of conductor 37
            sage: P.conductor()
            7
            sage: E = EllipticCurve('389a'); P = E.heegner_point(-7, 5); P.conductor()
            5
        """
        return self.__c

    def discriminant(self):
        """
        Return the discriminant of the quadratic imaginary field
        associated to this Heegner point.

        EXAMPLES::

            sage: heegner_point(389,-7,5).discriminant()
            -7
            sage: E = EllipticCurve('37a1'); P = E.kolyvagin_point(-67,7); P
            Kolyvagin point of discriminant -67 and conductor 7
             on elliptic curve of conductor 37
            sage: P.discriminant()
            -67
            sage: E = EllipticCurve('389a'); P = E.heegner_point(-7, 5); P.discriminant()
            -7
        """
        return self.__D

    @cached_method
    def quadratic_field(self):
        """
        Return the quadratic number field of discriminant `D`.

        EXAMPLES::

            sage: x = heegner_point(37,-7,5)
            sage: x.quadratic_field()
            Number Field in sqrt_minus_7 with defining polynomial x^2 + 7
             with sqrt_minus_7 = 2.645751311064591?*I

            sage: E = EllipticCurve('37a'); P = E.heegner_point(-40)
            sage: P.quadratic_field()
            Number Field in sqrt_minus_40 with defining polynomial x^2 + 40
             with sqrt_minus_40 = 6.324555320336759?*I
            sage: P.quadratic_field() is P.quadratic_field()
            True
            sage: type(P.quadratic_field())
            <class 'sage.rings.number_field.number_field.NumberField_quadratic_with_category'>
        """
        return self.ring_class_field().quadratic_field()

    @cached_method
    def quadratic_order(self):
        """
        Return the order in the quadratic imaginary field of conductor
        `c`, where `c` is the conductor of this Heegner point.

        EXAMPLES::

            sage: heegner_point(389,-7,5).quadratic_order()
            Order of conductor 10 generated by 5*sqrt_minus_7
             in Number Field in sqrt_minus_7 with defining polynomial x^2 + 7
             with sqrt_minus_7 = 2.645751311064591?*I
            sage: heegner_point(389,-7,5).quadratic_order().basis()
            [1, 5*sqrt_minus_7]

            sage: E = EllipticCurve('37a'); P = E.heegner_point(-40,11)
            sage: P.quadratic_order()
            Order of conductor 22 generated by 11*sqrt_minus_40
             in Number Field in sqrt_minus_40 with defining polynomial x^2 + 40
             with sqrt_minus_40 = 6.324555320336759?*I
            sage: P.quadratic_order().basis()
            [1, 11*sqrt_minus_40]
        """
        K = self.quadratic_field()
        return K.order([1,self.conductor()*K.gen()])

    @cached_method
    def ring_class_field(self):
        """
        Return the ring class field associated to this Heegner point.
        This is an extension `K_c` over `K`, where `K` is the
        quadratic imaginary field and `c` is the conductor associated
        to this Heegner point.  This Heegner point is defined over
        `K_c` and the Galois group `Gal(K_c/K)` acts transitively on
        the Galois conjugates of this Heegner point.

        EXAMPLES::

            sage: E = EllipticCurve('389a'); K.<a> = QuadraticField(-5)
            sage: len(K.factor(5))
            1
            sage: len(K.factor(23))
            2
            sage: E.heegner_point(-7, 5).ring_class_field().degree_over_K()
            6
            sage: E.heegner_point(-7, 23).ring_class_field().degree_over_K()
            22
            sage: E.heegner_point(-7, 5*23).ring_class_field().degree_over_K()
            132
            sage: E.heegner_point(-7, 5^2).ring_class_field().degree_over_K()
            30
            sage: E.heegner_point(-7, 7).ring_class_field().degree_over_K()
            7
        """
        return RingClassField(self.discriminant(), self.conductor())


##################################################################################
#
# Sets of Heegner points
#
##################################################################################

class HeegnerPoints(SageObject):
    """
    The set of Heegner points with given parameters.

    EXAMPLES::

        sage: H = heegner_points(389); H
        Set of all Heegner points on X_0(389)
        sage: type(H)
        <class 'sage.schemes.elliptic_curves.heegner.HeegnerPoints_level'>
        sage: isinstance(H, sage.schemes.elliptic_curves.heegner.HeegnerPoints)
        True
    """
    def __init__(self, N):
        """
        INPUT:

        - ``N`` -- level, a positive integer

        EXAMPLES::

            sage: heegner_points(37)
            Set of all Heegner points on X_0(37)
            sage: heegner_points(0)
            Traceback (most recent call last):
            ...
            ValueError: N must a positive integer
        """
        self.__N = ZZ(N)
        if self.__N <= 0:
            raise ValueError("N must a positive integer")

    def level(self):
        """
        Return the level `N` of the modular curve `X_0(N)`.

        EXAMPLES::

            sage: heegner_points(389).level()
            389
        """
        return self.__N


class HeegnerPoints_level(HeegnerPoints):
    """
    Return the infinite set of all Heegner points on `X_0(N)` for all
    quadratic imaginary fields.

    EXAMPLES::

        sage: H = heegner_points(11); H
        Set of all Heegner points on X_0(11)
        sage: type(H)
        <class 'sage.schemes.elliptic_curves.heegner.HeegnerPoints_level'>
        sage: loads(dumps(H)) == H
        True
    """
    def __eq__(self, other):
        """
        EXAMPLES::

            sage: H = heegner_points(11)
            sage: H == heegner_points(13)
            False
            sage: H == heegner_points(11)
            True
            sage: H == 0
            False
        """
        return isinstance(other, HeegnerPoints_level) and self.level() == other.level()

    def __ne__(self, other):
        """
        EXAMPLES::

            sage: H = heegner_points(11)
            sage: H != heegner_points(13)
            True
            sage: H != heegner_points(11)
            False
            sage: H != 0
            True
        """
        return not (self == other)

    def _repr_(self):
        """
        Return string representation of the set of Heegner points.

        EXAMPLES::

            sage: heegner_points(389)._repr_()
            'Set of all Heegner points on X_0(389)'
        """
        return "Set of all Heegner points on X_0(%s)" % self.level()

    def reduce_mod(self, ell):
        r"""
        Return object that allows for computation with Heegner points
        of level `N` modulo the prime `\ell`, represented using
        quaternion algebras.

        INPUT:

        - `\ell` -- prime

        EXAMPLES::

            sage: heegner_points(389).reduce_mod(7).quaternion_algebra()
            Quaternion Algebra (-1, -7) with base ring Rational Field
        """
        return HeegnerQuatAlg(self.level(), ell)

    def discriminants(self, n=10, weak=False):
        r"""
        Return the first `n` quadratic imaginary discriminants that
        satisfy the Heegner hypothesis for `N`.

        INPUT:

        - ``n`` -- nonnegative integer

        - ``weak`` -- boolean (default: ``False``); if ``True`` only require
          weak Heegner hypothesis, which is the same as usual but
          without the condition that `\gcd(D,N)=1`.

        EXAMPLES::

            sage: X = heegner_points(37)
            sage: X.discriminants(5)
            [-7, -11, -40, -47, -67]

        The default is 10::

            sage: X.discriminants()
            [-7, -11, -40, -47, -67, -71, -83, -84, -95, -104]
            sage: X.discriminants(15)
            [-7, -11, -40, -47, -67, -71, -83, -84, -95, -104, -107, -115, -120, -123, -127]

        The discriminant -111 satisfies only the weak Heegner hypothesis, since it
        is divisible by 37::

            sage: X.discriminants(15, weak=True)
            [-7, -11, -40, -47, -67, -71, -83, -84, -95, -104, -107, -111, -115, -120, -123]
        """
        N = self.level()
        n = ZZ(n)
        v = []
        D = ZZ(-4)
        while len(v) < n:
            D -= 1
            if satisfies_weak_heegner_hypothesis(N,D):
                # if not weak, then also require gcd(D,N)=1
                if not weak and D.gcd(N) != 1:
                    continue
                v.append(D)
        return v


class HeegnerPoints_level_disc(HeegnerPoints):
    """
    Set of Heegner points of given level and all conductors associated
    to a quadratic imaginary field.

    EXAMPLES::

        sage: H = heegner_points(389,-7); H
        Set of all Heegner points on X_0(389) associated to QQ[sqrt(-7)]
        sage: type(H)
        <class 'sage.schemes.elliptic_curves.heegner.HeegnerPoints_level_disc'>
        sage: H._repr_()
        'Set of all Heegner points on X_0(389) associated to QQ[sqrt(-7)]'
        sage: H.discriminant()
        -7
        sage: H.quadratic_field()
        Number Field in sqrt_minus_7 with defining polynomial x^2 + 7
         with sqrt_minus_7 = 2.645751311064591?*I
        sage: H.kolyvagin_conductors()
        [1, 3, 5, 13, 15, 17, 19, 31, 39, 41]

        sage: loads(dumps(H)) == H
        True
    """
    def __init__(self, N, D):
        """
        INPUT:

        - ``N`` -- positive integer

        - ``D`` -- negative fundamental discriminant

        EXAMPLES::

            sage: sage.schemes.elliptic_curves.heegner.HeegnerPoints_level_disc(37,-7)
            Set of all Heegner points on X_0(37) associated to QQ[sqrt(-7)]
        """
        HeegnerPoints.__init__(self, N)
        D = ZZ(D)
        if not satisfies_weak_heegner_hypothesis(N,D):
            raise ValueError("D (=%s) must satisfy the weak Heegner hypothesis for N (=%s)" % (D,N))
        self.__D = D

    def __eq__(self, other):
        """
        EXAMPLES::

            sage: H = heegner_points(389,-7)
            sage: H == heegner_points(389,-7)
            True
            sage: H == 0
            False
            sage: H == heegner_points(389,-11)
            False
        """
        return isinstance(other, HeegnerPoints_level_disc) and \
               self.level() == other.level() and self.__D == other.__D

    def __ne__(self, other):
        """
        EXAMPLES::

            sage: H = heegner_points(389,-7)
            sage: H != heegner_points(389,-7)
            False
            sage: H != 0
            True
            sage: H != heegner_points(389,-11)
            True
        """
        return not (self == other)

    def _repr_(self):
        """
        Return string representation of the set of Heegner points for a given
        quadratic field.

        EXAMPLES::

            sage: heegner_points(389,-7)._repr_()
            'Set of all Heegner points on X_0(389) associated to QQ[sqrt(-7)]'
        """
        return "Set of all Heegner points on X_0(%s) associated to QQ[sqrt(%s)]" % (
            self.level(), self.discriminant())

    def discriminant(self):
        r"""
        Return the discriminant of the quadratic imaginary extension `K`.

        EXAMPLES::

            sage: heegner_points(389,-7).discriminant()
            -7
        """
        return self.__D

    @cached_method
    def quadratic_field(self):
        r"""
        Return the quadratic imaginary field `K = \QQ(\sqrt{D})`.

        EXAMPLES::

            sage: E = EllipticCurve('389a'); K = E.heegner_point(-7,5).ring_class_field()
            sage: K.quadratic_field()
            Number Field in sqrt_minus_7 with defining polynomial x^2 + 7
             with sqrt_minus_7 = 2.645751311064591?*I
        """
        D = self.__D
        var = 'sqrt_minus_%s' % (-D)
        return QuadraticField(D,var)

    def kolyvagin_conductors(self, r=None, n=10, E=None, m=None):
        r"""
        Return the first `n` conductors that are squarefree products
        of distinct primes inert in the quadratic imaginary field
        `K = \QQ(\sqrt{D})`.  If `r` is specified, return only
        conductors that are a product of `r` distinct primes all inert
        in `K`.  If `r = 0`, always return the list ``[1]``,
        no matter what.

        If the optional elliptic curve `E` and integer `m` are given,
        then only include conductors `c` such that for each prime
        divisor `p` of `c` we have `m \mid \gcd(a_p(E), p+1)`.

        INPUT:

        - ``r`` -- (default: ``None``) nonnegative integer or ``None``

        - ``n`` -- positive integer

        - ``E`` -- an elliptic curve

        - ``m`` -- positive integer

        EXAMPLES::

            sage: H = heegner_points(389, -7)
            sage: H.kolyvagin_conductors(0)
            [1]
            sage: H.kolyvagin_conductors(1)
            [3, 5, 13, 17, 19, 31, 41, 47, 59, 61]
            sage: H.kolyvagin_conductors(1,15)
            [3, 5, 13, 17, 19, 31, 41, 47, 59, 61, 73, 83, 89, 97, 101]
            sage: H.kolyvagin_conductors(1,5)
            [3, 5, 13, 17, 19]
            sage: H.kolyvagin_conductors(1, 5, EllipticCurve('389a'), 3)
            [5, 17, 41, 59, 83]
            sage: H.kolyvagin_conductors(2, 5, EllipticCurve('389a'), 3)
            [85, 205, 295, 415, 697]
        """
        D = self.__D
        if not satisfies_weak_heegner_hypothesis(self.level(),D):
            raise ValueError("D must satisfy the weak Heegner hypothesis")
        n = ZZ(n)
        if n <= 0:
            raise ValueError("n must be a positive integer")
        if r is not None:
            r = ZZ(r)
            if r < 0:
                raise ValueError("n must be a nonnegative integer")
        if r == 0:
            return [ZZ(1)]

        c = ZZ(1)
        v = []
        N = self.level()

        if E is not None:
            m = ZZ(m)

        while len(v) < n:
            if is_kolyvagin_conductor(N, E, D, r, m, c):
                v.append(c)
            c += 1

        return v


def is_kolyvagin_conductor(N, E, D, r, n, c):
    r"""
    Return ``True`` if `c` is a Kolyvagin conductor for level `N`,
    discriminant `D`, mod `n`, etc., i.e., `c` is divisible by exactly
    `r` prime factors, is coprime to `ND`, each prime dividing `c` is
    inert, and if `E` is not ``None`` then `n | \gcd(p+1, a_p(E))`
    for each prime `p` dividing `c`.

    INPUT:

    - ``N`` -- level (positive integer)

    - ``E`` -- elliptic curve or ``None``

    - ``D`` -- negative fundamental discriminant

    - ``r`` -- number of prime factors (nonnegative integer) or ``None``

    - ``n`` -- torsion order (i.e., do we get class in `(E(K_c)/n E(K_c))^{Gal(K_c/K)}`?)

    - ``c`` -- conductor (positive integer)

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.heegner import is_kolyvagin_conductor
        sage: is_kolyvagin_conductor(389, None, -7, 1, None, 5)
        True
        sage: is_kolyvagin_conductor(389, None, -7, 1, None, 7)
        False
        sage: is_kolyvagin_conductor(389, None, -7, 1, None, 11)
        False
        sage: is_kolyvagin_conductor(389, EllipticCurve('389a'), -7, 1, 3, 5)
        True
        sage: is_kolyvagin_conductor(389, EllipticCurve('389a'), -7, 1, 11, 5)
        False
    """
    ND = N*D
    if ND.gcd(c) != 1:
        return False
    if not c.is_squarefree():
        return False
    P = c.prime_factors()
    if r is not None and len(P) != r:
        return False
    # check that each prime in P is inert in K
    for p in P:
        if D.kronecker(p) != -1:
            return False
    if E is not None and n is not None:
        for p in P:
            if (p+1).gcd(E.ap(p)) % n != 0:
                return False
    return True


class HeegnerPoints_level_disc_cond(HeegnerPoints_level, HeegnerPoints_level_disc):
    """
    The set of Heegner points of given level, discriminant, and conductor.

    EXAMPLES::

        sage: H = heegner_points(389,-7,5); H
        All Heegner points of conductor 5 on X_0(389) associated to QQ[sqrt(-7)]
        sage: type(H)
        <class 'sage.schemes.elliptic_curves.heegner.HeegnerPoints_level_disc_cond'>
        sage: H.discriminant()
        -7
        sage: H.level()
        389

        sage: len(H.points())
        12
        sage: H.points()[0]
        Heegner point 5/778*sqrt(-7) - 147/778 of discriminant -7 and conductor 5 on X_0(389)
        sage: H.betas()
        (147, 631)

        sage: H.quadratic_field()
        Number Field in sqrt_minus_7 with defining polynomial x^2 + 7
         with sqrt_minus_7 = 2.645751311064591?*I
        sage: H.ring_class_field()
        Ring class field extension of QQ[sqrt(-7)] of conductor 5

        sage: H.kolyvagin_conductors()
        [1, 3, 5, 13, 15, 17, 19, 31, 39, 41]
        sage: H.satisfies_kolyvagin_hypothesis()
        True

        sage: H = heegner_points(389,-7,5)
        sage: loads(dumps(H)) == H
        True
    """
    def __init__(self, N, D, c=ZZ(1)):
        """
        Create set of Heegner points.

        INPUT:

        - ``N`` -- positive integer (the level)

        - ``D`` -- negative fundamental discriminant

        - ``c`` -- conductor (default: 1)

        EXAMPLES::

            sage: H = heegner_points(389,-7,5); H
            All Heegner points of conductor 5 on X_0(389) associated to QQ[sqrt(-7)]
            sage: type(H)
            <class 'sage.schemes.elliptic_curves.heegner.HeegnerPoints_level_disc_cond'>
        """
        HeegnerPoints_level.__init__(self, N)
        HeegnerPoints_level_disc.__init__(self, N, D)
        self.__c = ZZ(c)

    def __eq__(self, other):
        """
        EXAMPLES::

            sage: H = heegner_points(389,-7, 3)
            sage: H == heegner_points(389,-7, 3)
            True
            sage: H == heegner_points(389,-7, 1)
            False
            sage: H == 0
            False
        """
        return isinstance(other, HeegnerPoints_level_disc_cond) and \
               self.level() == other.level() and self.discriminant() == other.discriminant() \
               and self.conductor() == other.conductor()

    def __ne__(self, other):
        """
        EXAMPLES::

            sage: H = heegner_points(389,-7, 3)
            sage: H != heegner_points(389,-7, 3)
            False
            sage: H != heegner_points(389,-7, 1)
            True
            sage: H != 0
            True
        """
        return not (self == other)

    def _repr_(self):
        """
        Return string representation of this set of Heegner points.

        EXAMPLES::

            sage: H = heegner_points(37,-7,5); H._repr_()
            'All Heegner points of conductor 5 on X_0(37) associated to QQ[sqrt(-7)]'
        """
        return "All Heegner points of conductor %s on X_0(%s) associated to QQ[sqrt(%s)]" % (
            self.conductor(), self.level(), self.discriminant())

    def conductor(self):
        """
        Return the level of the conductor.

        EXAMPLES::

            sage: heegner_points(389,-7,5).conductor()
            5
        """
        return self.__c

    @cached_method
    def satisfies_kolyvagin_hypothesis(self):
        """
        Return ``True`` if ``self`` satisfies the Kolyvagin hypothesis, i.e.,
        that each prime dividing the conductor `c` of ``self`` is inert in
        `K` and coprime to `ND`.

        EXAMPLES:

        The prime 5 is inert, but the prime 11 is not::

            sage: heegner_points(389,-7,5).satisfies_kolyvagin_hypothesis()
            True
            sage: heegner_points(389,-7,11).satisfies_kolyvagin_hypothesis()
            False
        """
        return is_kolyvagin_conductor(N=self.level(), E=None, D=self.discriminant(),
                                      r=None, n=None, c=self.conductor())

    @cached_method
    def ring_class_field(self):
        """
        Return the ring class field associated to this set of Heegner
        points.  This is an extension `K_c` over `K`, where `K` is the
        quadratic imaginary field and `c` the conductor associated to
        this Heegner point.  This Heegner point is defined over `K_c`
        and the Galois group `Gal(K_c/K)` acts transitively on the
        Galois conjugates of this Heegner point.

        EXAMPLES::

            sage: heegner_points(389,-7,5).ring_class_field()
            Ring class field extension of QQ[sqrt(-7)] of conductor 5
        """
        return RingClassField(self.discriminant(), self.conductor())

    def __getitem__(self, i):
        """
        Return the `i`-th Heegner point.

        EXAMPLES::

            sage: H = heegner_points(389,-7,5)
            sage: len(H)
            12
            sage: H[0]
            Heegner point 5/778*sqrt(-7) - 147/778 of discriminant -7 and conductor 5
             on X_0(389)
            sage: H[-1]
            Heegner point 5/5446*sqrt(-7) - 757/778 of discriminant -7 and conductor 5
             on X_0(389)
        """
        return self.points()[i]

    def __len__(self):
        """
        Return the number of Heegner points.

        EXAMPLES::

            sage: len(heegner_points(389,-7,5))
            12

        When the conductor is 1 the length is a power of 2 (number of
        square roots of `D` mod `4N` reduced mod `2N`) times the class
        number::

            sage: len(heegner_points(389,-20,1))
            4
            sage: QQ[sqrt(-20)].class_number()
            2
        """
        return len(self.points())

    @cached_method
    def betas(self):
        """
        Return the square roots of `D c^2` modulo `4 N` all reduced
        mod `2 N`, without multiplicity.

        EXAMPLES::

            sage: X = heegner_points(45,-11,1); X
            All Heegner points of conductor 1 on X_0(45) associated to QQ[sqrt(-11)]
            sage: [x.quadratic_form() for x in X]
            [45*x^2 + 13*x*y + y^2,
             45*x^2 + 23*x*y + 3*y^2,
             45*x^2 + 67*x*y + 25*y^2,
             45*x^2 + 77*x*y + 33*y^2]
            sage: X.betas()
            (13, 23, 67, 77)
            sage: X.points(13)
            (Heegner point 1/90*sqrt(-11) - 13/90 of discriminant -11 on X_0(45),)
            sage: [x.quadratic_form() for x in X.points(13)]
            [45*x^2 + 13*x*y + y^2]
        """
        c = self.__c
        D = self.discriminant()*c*c
        N = self.level()
        R = Integers(4*N)
        m = 2*N
        return tuple(sorted( {a % m for a in R(D).sqrt(all=True)} ))

    @cached_method
    def points(self, beta=None):
        r"""
        Return the Heegner points in ``self``.  If `\beta` is given,
        return only those Heegner points with given `\beta`, i.e.,
        whose quadratic form has `B` congruent to `\beta` modulo `2 N`.

        Use ``self.beta()`` to get a list of betas.

        EXAMPLES::

            sage: H = heegner_points(389, -7, 5); H
            All Heegner points of conductor 5 on X_0(389) associated to QQ[sqrt(-7)]
            sage: H.points()
            (Heegner point 5/778*sqrt(-7) - 147/778 of discriminant -7 and conductor 5
              on X_0(389),
             ...,
             Heegner point 5/5446*sqrt(-7) - 757/778 of discriminant -7 and conductor 5
              on X_0(389))
            sage: H.betas()
            (147, 631)
            sage: [x.tau() for x in H.points(147)]
            [5/778*sqrt_minus_7 - 147/778, 5/1556*sqrt_minus_7 - 147/1556,
             5/1556*sqrt_minus_7 - 925/1556, 5/3112*sqrt_minus_7 - 1703/3112,
             5/3112*sqrt_minus_7 - 2481/3112, 5/5446*sqrt_minus_7 - 21/778]

            sage: [x.tau() for x in H.points(631)]
            [5/778*sqrt_minus_7 - 631/778, 5/1556*sqrt_minus_7 - 631/1556,
             5/1556*sqrt_minus_7 - 1409/1556, 5/3112*sqrt_minus_7 - 631/3112,
             5/3112*sqrt_minus_7 - 1409/3112, 5/5446*sqrt_minus_7 - 757/778]

        The result is cached and is a tuple (since it is immutable)::

            sage: H.points() is H.points()
            True
            sage: type(H.points())
            <... 'tuple'>
        """
        if beta is None:
            SDN = self.betas()
            return tuple(sorted(sum([list(self.points(b)) for b in SDN], [])))

        c = self.conductor()
        N = self.level()
        D = self.discriminant()
        b = ZZ(beta) % (2*N)

        disc = D*c*c

        U = []
        R = []
        h = self.ring_class_field().degree_over_K()
        a = 1
        while len(U) < h:
            if c.gcd(a) != 1:
                a += 1
                continue
            # todo (optimize) -- replace for over all s with for over solution set
            y = ZZ((b*b - disc)/(4*N))
            for s in Integers(a):
                if N*s*s + b*s + y == 0:
                    s = s.lift()
                    f = (a*N, b+2*N*s, ZZ( ((b + 2*N*s)**2 - disc)/(4*a*N)) )
                    g = BinaryQF(f).reduced_form()
                    assert g.discriminant() == disc
                    if g not in U:
                        U.append(g)
                        R.append(HeegnerPointOnX0N(N,D,c,f))
                        if len(U) >= h:
                            break
            a += 1
        return tuple(sorted(R))

    def plot(self, *args, **kwds):
        """
        Return plot of all the representatives in the upper half
        plane of the Heegner points in this set of Heegner points.

        The inputs to this function get passed onto the point command.

        EXAMPLES::

            sage: heegner_points(389,-7,5).plot(pointsize=50, rgbcolor='red')           # needs sage.plot
            Graphics object consisting of 12 graphics primitives
            sage: heegner_points(53,-7,15).plot(pointsize=50, rgbcolor='purple')        # needs sage.plot
            Graphics object consisting of 48 graphics primitives
        """
        return sum(z.plot(*args, **kwds) for z in self)


class HeegnerPointOnX0N(HeegnerPoint):
    r"""
    A Heegner point as a point on the modular curve `X_0(N)`, which we
    view as the upper half plane modulo the action of `\Gamma_0(N)`.

    EXAMPLES::

        sage: x = heegner_point(37, -7, 5); x
        Heegner point 5/74*sqrt(-7) - 11/74 of discriminant -7 and conductor 5 on X_0(37)
        sage: type(x)
        <class 'sage.schemes.elliptic_curves.heegner.HeegnerPointOnX0N'>
        sage: x.level()
        37
        sage: x.conductor()
        5
        sage: x.discriminant()
        -7
        sage: x.quadratic_field()
        Number Field in sqrt_minus_7 with defining polynomial x^2 + 7
         with sqrt_minus_7 = 2.645751311064591?*I
        sage: x.quadratic_form()
        37*x^2 + 11*x*y + 2*y^2
        sage: x.quadratic_order()
        Order of conductor 10 generated by 5*sqrt_minus_7
         in Number Field in sqrt_minus_7 with defining polynomial x^2 + 7
         with sqrt_minus_7 = 2.645751311064591?*I
        sage: x.tau()
        5/74*sqrt_minus_7 - 11/74
        sage: loads(dumps(x)) == x
        True
    """
    def __init__(self, N, D, c=ZZ(1), f=None, check=True):
        r"""
        INPUT:

        - ``N`` -- positive integer

        - ``D`` -- fundamental discriminant, a negative integer

        - ``c`` -- conductor, a positive integer coprime to `N`

        - ``f`` -- binary quadratic form, 3-tuple `(A,B,C)` of coefficients
          of `AX^2 + BXY + CY^2`, or element of quadratic imaginary
          field `\QQ(\sqrt{D})` in the upper half plane

        - ``check`` -- boolean (default: ``True``); should not be used
          except internally

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.heegner import HeegnerPointOnX0N
            sage: x = heegner_point(389, -7, 5); x
            Heegner point 5/778*sqrt(-7) - 147/778 of discriminant -7 and conductor 5
             on X_0(389)
            sage: type(x)
            <class 'sage.schemes.elliptic_curves.heegner.HeegnerPointOnX0N'>
            sage: HeegnerPointOnX0N(389, -7, 5, None, check=False)
            Heegner point 5/778*sqrt(-7) - 147/778 of discriminant -7 and conductor 5
             on X_0(389)
        """
        if check:
            N = ZZ(N)
            D = ZZ(D)
            c = ZZ(c)
            if c.gcd(N) != 1:
                raise ValueError("conductor c (=%s) must be coprime to N (=%s)" % (c, N))
            if not satisfies_weak_heegner_hypothesis(N, D):
                raise ValueError("N (=%s) and D (=%s) must satisfy the Heegner hypothesis" % (N, D))
            if f is not None:
                if isinstance(f, tuple):
                    if len(f) != 3:
                        raise ValueError("if f is a tuple, it must have length 3")
                    f = tuple(ZZ(a) for a in f)
                elif isinstance(f, BinaryQF):
                    # convert from BinaryQF
                    f = tuple(f)
                elif isinstance(f, NumberFieldElement_base):
                    # tau = number field element
                    g = f.minpoly()
                    if g.degree() != 2:
                        raise TypeError("number field element f must have degree 2")
                    g *= g.denominator()  # make integral
                    f = (ZZ(g[2]), ZZ(g[1]), ZZ(g[0]))
                else:
                    raise TypeError("f must be a 3-tuple, quadratic form, or element of the upper half plane")
                A, B, C = f
                if B*B - 4*A*C != D*c*c:
                    raise ValueError("f (=%s) must have discriminant %s" % (f, D*c*c))
        HeegnerPoint.__init__(self, N, D, c)
        if f is None:
            # We know that N|A, so A = N is optimal.
            A = N
            B = ZZ(Integers(4*N)(D*c*c).sqrt(extend=False) % (2*N))
            C = ZZ((B*B - D*c*c)/(4*A))
            f = (A,B,C)
        self.__f = f

    def __hash__(self):
        """
        The hash is obtained from the hash provided by :class:`HeegnerPoint`,
        together with the reduced quadratic form.

        EXAMPLES::

            sage: x = heegner_point(37, -7, 5)
            sage: from sage.schemes.elliptic_curves.heegner import HeegnerPoint
            sage: hash(x) == hash( (HeegnerPoint.__hash__(x), x.reduced_quadratic_form()) )
            True
        """
        return hash((HeegnerPoint.__hash__(self), self.reduced_quadratic_form()))

    def __richcmp__(self, x, op):
        """
        Compare two Heegner points with character.

        EXAMPLES::

            sage: x1 = EllipticCurve('389a').heegner_point(-7).heegner_point_on_X0N()
            sage: x5 = EllipticCurve('389a').heegner_point(-7, 5).heegner_point_on_X0N()
            sage: x1 == x1
            True
            sage: x1 < x5
            True
            sage: x5 > x1
            True
        """
        if not isinstance(x, HeegnerPointOnX0N):
            return NotImplemented
        return richcmp((self.level(), self.discriminant(),
                        self.conductor(), self.__f),
                       (x.level(), x.discriminant(),
                        x.conductor(), x.__f), op)

    def _repr_(self):
        """
        Return string representation of this Heegner point.

        EXAMPLES::

            sage: x = heegner_point(37, -7, 5); x._repr_()
            'Heegner point 5/74*sqrt(-7) - 11/74 of discriminant -7 and conductor 5 on X_0(37)'
        """
        c = self.conductor()
        s = " and conductor %s" % c if c != 1 else ""
        N = self.level()
        D = self.discriminant()
        tau = repr(self.tau()).replace('sqrt_minus_%s' % (-D),'sqrt(%s)' % D)
        return "Heegner point %s of discriminant %s%s on X_0(%s)" % (tau, D, s, N)

    def atkin_lehner_act(self, Q=None):
        r"""
        Given an integer Q dividing the level N such that `\gcd(Q, N/Q) = 1`, return the
        image of this Heegner point under the Atkin-Lehner operator `W_Q`.

        INPUT:

        - ``Q`` -- positive divisor of `N`; if not given, default to `N`

        EXAMPLES::

            sage: x = heegner_point(389, -7, 5)
            sage: x.atkin_lehner_act()
            Heegner point 5/199168*sqrt(-7) - 631/199168 of discriminant -7
             and conductor 5 on X_0(389)

            sage: x = heegner_point(45, D=-11, c=1); x
            Heegner point 1/90*sqrt(-11) - 13/90 of discriminant -11 on X_0(45)
            sage: x.atkin_lehner_act(5)
            Heegner point 1/90*sqrt(-11) + 23/90 of discriminant -11 on X_0(45)
            sage: y = x.atkin_lehner_act(9); y
            Heegner point 1/90*sqrt(-11) - 23/90 of discriminant -11 on X_0(45)
            sage: z = y.atkin_lehner_act(9); z
            Heegner point 1/90*sqrt(-11) - 13/90 of discriminant -11 on X_0(45)
            sage: z == x
            True
        """
        N = self.level()
        if Q is None:
            Q = N
        if Q == 1:
            return self  # trivial special case
        g, u, v = xgcd(Q * Q, -N)
        if g != Q:
            raise ValueError("Q must divide N and be coprime to N/Q")
        tau = self.tau()
        WQ_tau = ((u * Q * tau + v) / (N * tau + Q))
        return HeegnerPointOnX0N(N, self.discriminant(), self.conductor(),
                                 f=WQ_tau, check=True)

    @cached_method
    def quadratic_form(self):
        """
        Return the integral primitive positive-definite binary
        quadratic form associated to this Heegner point.

        EXAMPLES::

            sage: heegner_point(389, -7, 5).quadratic_form()
            389*x^2 + 147*x*y + 14*y^2
        """
        # It is good/important that this return a copy, since
        # BinaryQF's stupidly are mutable and cannot be made immutable.
        # In particular, they have a stupid reduce method that changes
        # them in place.
        return BinaryQF(self.__f)

    def reduced_quadratic_form(self):
        """
        Return reduced binary quadratic corresponding to this Heegner point.

        EXAMPLES::

            sage: x = heegner_point(389, -7, 5)
            sage: x.quadratic_form()
            389*x^2 + 147*x*y + 14*y^2
            sage: x.reduced_quadratic_form()
            4*x^2 - x*y + 11*y^2
        """
        return self.quadratic_form().reduced_form()

    @cached_method
    def tau(self):
        """
        Return an element ``tau`` in the upper half plane that corresponds
        to this particular Heegner point.

        Actually, ``tau`` is in the quadratic imaginary field K associated
        to this Heegner point.

        EXAMPLES::

            sage: x = heegner_point(37, -7, 5); tau = x.tau(); tau
            5/74*sqrt_minus_7 - 11/74
            sage: 37 * tau.minpoly()
            37*x^2 + 11*x + 2
            sage: x.quadratic_form()
            37*x^2 + 11*x*y + 2*y^2
        """
        K = self.quadratic_field()
        d = K.gen() * self.conductor()
        A, B, _ = self.__f
        return (-B + d) / (2 * A)

    def map_to_curve(self, E):
        """
        Return the image of this Heegner point on the elliptic curve
        `E`, which must also have conductor `N`, where `N` is the
        level of ``self``.

        EXAMPLES::

            sage: x = heegner_point(389, -7, 5); x
            Heegner point 5/778*sqrt(-7) - 147/778 of discriminant -7
             and conductor 5 on X_0(389)
            sage: y = x.map_to_curve(EllipticCurve('389a')); y
            Heegner point of discriminant -7 and conductor 5
             on elliptic curve of conductor 389
            sage: y.curve().cremona_label()
            '389a1'
            sage: y.heegner_point_on_X0N()
            Heegner point 5/778*sqrt(-7) - 147/778 of discriminant -7
             and conductor 5 on X_0(389)

        You can also directly apply the modular parametrization of the elliptic curve::

            sage: x = heegner_point(37,-7); x
            Heegner point 1/74*sqrt(-7) - 17/74 of discriminant -7 on X_0(37)
            sage: E = EllipticCurve('37a'); phi = E.modular_parametrization()
            sage: phi(x)
            Heegner point of discriminant -7 on elliptic curve of conductor 37
        """
        return HeegnerPointOnEllipticCurve(E, self)

    @cached_method
    def galois_orbit_over_K(self):
        r"""
        Return the `Gal(K_c/K)`-orbit of this Heegner point.

        EXAMPLES::

            sage: x = heegner_point(389, -7, 3); x
            Heegner point 3/778*sqrt(-7) - 223/778 of discriminant -7
             and conductor 3 on X_0(389)
            sage: x.galois_orbit_over_K()
            [Heegner point 3/778*sqrt(-7) - 223/778 of discriminant -7 and conductor 3 on X_0(389),
             Heegner point 3/1556*sqrt(-7) - 223/1556 of discriminant -7 and conductor 3 on X_0(389),
             Heegner point 3/1556*sqrt(-7) - 1001/1556 of discriminant -7 and conductor 3 on X_0(389),
             Heegner point 3/3112*sqrt(-7) - 223/3112 of discriminant -7 and conductor 3 on X_0(389)]
        """
        c = self.conductor()
        N = self.level()
        D = self.discriminant()
        b = self.__f[1] % (2*N)  # B

        disc = D*c*c

        U = []
        R = []
        h = self.ring_class_field().degree_over_K()
        a = 1
        while len(U) < h:
            if c.gcd(a) != 1:
                a += 1
                continue
            # todo (optimize) -- replace for over all s with for over solution set
            y = ZZ((b*b - disc)/(4*N))
            for s in Integers(a):
                if N*s*s + b*s + y == 0:
                    s = s.lift()
                    f = (a*N, b+2*N*s, ZZ( ((b + 2*N*s)**2 - disc)/(4*a*N)) )
                    g = BinaryQF(f).reduced_form()
                    assert g.discriminant() == disc
                    if g not in U:
                        U.append(g)
                        R.append(HeegnerPointOnX0N(N,D,c,f))
            a += 1
        return R

    def plot(self, **kwds):
        r"""
        Draw a point at `(x,y)` where this Heegner point is
        represented by the point `\tau = x + i y` in the upper half
        plane.

        The ``kwds`` get passed onto the point plotting command.

        EXAMPLES::

            sage: heegner_point(389,-7,1).plot(pointsize=50)
            Graphics object consisting of 1 graphics primitive
        """
        from sage.plot.all import point
        return point(CDF(self.tau()), **kwds)


class HeegnerPointOnEllipticCurve(HeegnerPoint):
    """
    A Heegner point on a curve associated to an order in a quadratic
    imaginary field.

    EXAMPLES::

        sage: E = EllipticCurve('37a'); P = E.heegner_point(-7,5); P
        Heegner point of discriminant -7 and conductor 5 on elliptic curve of conductor 37
        sage: type(P)
        <class 'sage.schemes.elliptic_curves.heegner.HeegnerPointOnEllipticCurve'>
    """
    def __init__(self, E, x, check=True):
        r"""
        INPUT:

        - ``E`` -- an elliptic curve over the rational numbers

        - ``x`` -- Heegner point on `X_0(N)`

        - ``check`` -- boolean (default: ``True``); if ``True``, ensure that `D`,
          `c` are of type Integer and define a Heegner point on `E`

        EXAMPLES::

            sage: x = heegner_point(389,-7,5)
            sage: E = EllipticCurve('389a')
            sage: sage.schemes.elliptic_curves.heegner.HeegnerPointOnEllipticCurve(E, x)
            Heegner point of discriminant -7 and conductor 5 on elliptic curve
             of conductor 389
        """
        if check:
            if E.conductor() != x.level():
                raise ValueError("conductor of curve must equal level of Heegner point")
        self.__E = E
        self.__x = x
        HeegnerPoint.__init__(self, x.level(), x.discriminant(), x.conductor())

    @cached_method
    def satisfies_kolyvagin_hypothesis(self, n=None):
        r"""
        Return ``True`` if this Heegner point and `n` satisfy the
        Kolyvagin hypothesis, i.e., that each prime dividing the
        conductor `c` of ``self`` is inert in K and coprime to `ND`.
        Moreover, if `n` is not ``None``, also check that for each prime
        `p` dividing `c` we have that `n | \gcd(a_p(E), p+1)`.

        INPUT:

        - ``n`` -- positive integer

        EXAMPLES::

            sage: EllipticCurve('389a').heegner_point(-7).satisfies_kolyvagin_hypothesis()
            True
            sage: EllipticCurve('389a').heegner_point(-7, 5).satisfies_kolyvagin_hypothesis()
            True
            sage: EllipticCurve('389a').heegner_point(-7, 11).satisfies_kolyvagin_hypothesis()
            False
        """
        if n is not None:
            n = ZZ(n)
            if n <= 0:
                raise ValueError("n must be a positive integer")
        return is_kolyvagin_conductor(N=self.level(), E=self.__E, D=self.discriminant(),
                                      r=None, n=n, c=self.conductor())

    def __hash__(self):
        """
        The hash value is obtained from the elliptic curve and the Heegner
        point on `X_0(N)`.

        EXAMPLES::

            sage: x = EllipticCurve('389a').heegner_point(-7, 5)
            sage: hash(x) == hash( (x.curve(), x.heegner_point_on_X0N()) )
            True
        """
        return hash((self.__E, self.__x))

    def __eq__(self, right):
        """
        EXAMPLES::

            sage: y1 = EllipticCurve('389a').heegner_point(-7)
            sage: y5 = EllipticCurve('389a').heegner_point(-7, 5)
            sage: y1 == y1
            True
            sage: y5 == y5
            True
            sage: y1 == y5
            False
            sage: y1 == 10
            False
        """
        return isinstance(right, HeegnerPointOnEllipticCurve) and \
               (self.__E, self.__x) == (right.__E, right.__x)

    def __ne__(self, other):
        """
        EXAMPLES::

            sage: y1 = EllipticCurve('389a').heegner_point(-7)
            sage: y5 = EllipticCurve('389a').heegner_point(-7, 5)
            sage: y1 != y1
            False
            sage: y5 != y5
            False
            sage: y1 != y5
            True
            sage: y1 != 10
            True
        """
        return not (self == other)

    def _repr_(self):
        """
        Return string representation of this Heegner point.

        EXAMPLES::

            sage: E = EllipticCurve('389a'); P = E.heegner_point(-7, 97)
            sage: P._repr_()
            'Heegner point of discriminant -7 and conductor 97 on elliptic curve of conductor 389'
        """
        s = " and conductor %s" % self.conductor() if self.conductor() != 1 else ""
        N = self.__E.conductor()
        return "Heegner point of discriminant %s%s on elliptic curve of conductor %s" % (self.discriminant(), s, N)

    def heegner_point_on_X0N(self):
        r"""
        Return Heegner point on `X_0(N)` that maps to this Heegner point on `E`.

        EXAMPLES::

            sage: E = EllipticCurve('37a'); P = E.heegner_point(-7,5); P
            Heegner point of discriminant -7 and conductor 5 on elliptic curve
             of conductor 37
            sage: P.heegner_point_on_X0N()
            Heegner point 5/74*sqrt(-7) - 11/74 of discriminant -7 and conductor 5
             on X_0(37)
        """
        return self.__x

    def map_to_complex_numbers(self, prec=53):
        """
        Return the point in the subfield `M` of the complex numbers
        (well defined only modulo the period lattice) corresponding to
        this Heegner point.

        EXAMPLES:

        We compute a nonzero Heegner point over a ring class field on
        a curve of rank 2::

            sage: E = EllipticCurve('389a'); y = E.heegner_point(-7,5)
            sage: y.map_to_complex_numbers()
            1.49979679635196 + 0.369156204821526*I
            sage: y.map_to_complex_numbers(100)
            1.4997967963519640592142411892 + 0.36915620482152626830089145962*I
            sage: y.map_to_complex_numbers(10)
            1.5 + 0.37*I

        Here we see that the Heegner point is 0 since it lies in the
        lattice::

            sage: E = EllipticCurve('389a'); y = E.heegner_point(-7)
            sage: y.map_to_complex_numbers(10)
            0.0034 - 3.9*I
            sage: y.map_to_complex_numbers()
            4.71844785465692e-15 - 3.94347540310330*I
            sage: E.period_lattice().basis()
            (2.49021256085505, 1.97173770155165*I)
            sage: 2*E.period_lattice().basis()[1]
            3.94347540310330*I

        You can also directly coerce to the complex field::

            sage: E = EllipticCurve('389a'); y = E.heegner_point(-7)
            sage: z = ComplexField(100)(y); z # real part approx. 0
            -... - 3.9434754031032964088448153963*I
            sage: E.period_lattice().elliptic_exponential(z)
            (0.00000000000000000000000000000 : 1.0000000000000000000000000000 : 0.00000000000000000000000000000)
        """
        phi = self.__E.modular_parametrization()
        tau = self.heegner_point_on_X0N().tau()
        return phi.map_to_complex_numbers(tau, prec)

    def _complex_mpfr_field_(self, C):
        """
        Used internally for coercing Heegner point to a complex field.

        EXAMPLES::

            sage: E = EllipticCurve('37a'); y = E.heegner_point(-7)
            sage: CC(y)                          # indirect doctest
            0.929592715285395 - 1.22569469099340*I
            sage: ComplexField(100)(y)
            0.92959271528539567440519934446 - 1.2256946909933950304271124159*I
        """
        phi = self.__E.modular_parametrization()
        tau = C(self.heegner_point_on_X0N().tau())
        return phi.map_to_complex_numbers(tau)

    @cached_method
    def kolyvagin_point(self):
        """
        Return the Kolyvagin point corresponding to this Heegner
        point.

        This is the point obtained by applying the Kolyvagin
        operator `J_c I_c` in the group ring of the Galois group to
        this Heegner point.   It is a point that defines an element
        of `H^1(K, E[n])`, under certain hypotheses on `n`.

        EXAMPLES::

            sage: E = EllipticCurve('37a1'); y = E.heegner_point(-7); y
            Heegner point of discriminant -7 on elliptic curve of conductor 37
            sage: P = y.kolyvagin_point(); P
            Kolyvagin point of discriminant -7 on elliptic curve of conductor 37
            sage: P.numerical_approx()  # abs tol 1e-15
            (-3.36910401903861e-16 - 2.22076195576076e-16*I
             : 3.33066907387547e-16 + 2.22076195576075e-16*I : 1.00000000000000)
        """
        return KolyvaginPoint(self)

    @cached_method
    def _trace_index(self, *args, **kwds):
        """
        Return index of the trace of this Heegner point down to `K` in
        the group of `K`-rational points.

        IMPORTANT: See the help for ``E=self.curve(); E.index?`` for
        the inputs to this function and more details about what is
        computed.  In particular, the returned index can be off at 2.

        OUTPUT: ``Integer`` -- returns an integer

        EXAMPLES::

            sage: E = EllipticCurve('77a1')
            sage: P = E.heegner_point(-19); y = P._trace_numerical_conductor_1()
            sage: [c.real() for c in y]
            [-1.2...e-16, -1.00000000000000, 1.00000000000000]
            sage: -2*E.gens()[0]
            (0 : -1 : 1)
            sage: P._trace_index()
            2

            sage: P = E.heegner_point(-68); P
            Heegner point of discriminant -68 on elliptic curve of conductor 77
            sage: N(P)
            (0.219223593595584 - 1.87443160153148*I
             : -1.34232921921325 - 1.52356748877889*I : 1.00000000000000)
            sage: P._trace_index()
            0
        """
        if self.conductor() != 1:
            raise ValueError("conductor of Heegner point must be 1")
        i = self.__E.heegner_index(self.discriminant(), *args, **kwds)
        lower = i.lower().round()
        upper = i.upper().round()
        if lower == upper:
            return lower
        # here we would say raise precision somehow.
        raise NotImplementedError("unable to compute index")

    def curve(self):
        """
        Return the elliptic curve on which this is a Heegner point.

        EXAMPLES::

            sage: E = EllipticCurve('389a'); P = E.heegner_point(-7, 5)
            sage: P.curve()
            Elliptic Curve defined by y^2 + y = x^3 + x^2 - 2*x over Rational Field
            sage: P.curve() is E
            True
        """
        return self.__E

    @cached_method
    def quadratic_form(self):
        """
        Return the integral primitive positive definite binary
        quadratic form associated to this Heegner point.

        EXAMPLES::

            sage: EllipticCurve('389a').heegner_point(-7, 5).quadratic_form()
            389*x^2 + 147*x*y + 14*y^2

            sage: P = EllipticCurve('389a').heegner_point(-7, 5, (778,925,275)); P
            Heegner point of discriminant -7 and conductor 5 on elliptic curve
             of conductor 389
            sage: P.quadratic_form()
            778*x^2 + 925*x*y + 275*y^2
        """
        return self.__x.quadratic_form()

    @cached_method
    def numerical_approx(self, prec=53, algorithm=None):
        """
        Return a numerical approximation to this Heegner point
        computed using a working precision of prec bits.

        .. WARNING::

            The answer is *not* provably correct to prec bits!  A
            priori, due to rounding and other errors, it is possible that
            not a single digit is correct.

        INPUT:

        - prec     -- (default: ``None``) the working precision

        EXAMPLES::

            sage: E = EllipticCurve('37a'); P = E.heegner_point(-7); P
            Heegner point of discriminant -7 on elliptic curve of conductor 37
            sage: P.numerical_approx()  # abs tol 1e-15
            (-3.36910401903861e-16 - 2.22076195576076e-16*I
             : 3.33066907387547e-16 + 2.22076195576075e-16*I : 1.00000000000000)
            sage: P.numerical_approx(10)  # expect random digits
            (0.0030 - 0.0028*I : -0.0030 + 0.0028*I : 1.0)
            sage: P.numerical_approx(100)[0]  # expect random digits
            8.4...e-31 + 6.0...e-31*I
            sage: E = EllipticCurve('37a'); P = E.heegner_point(-40); P
            Heegner point of discriminant -40 on elliptic curve of conductor 37
            sage: P.numerical_approx()  # abs tol 1e-14
            (-3.15940603400359e-16 + 1.41421356237309*I
             : 1.00000000000000 - 1.41421356237309*I : 1.00000000000000)

        A rank 2 curve, where all Heegner points of conductor 1 are 0::

            sage: E = EllipticCurve('389a'); E.rank()
            2
            sage: P = E.heegner_point(-7); P
            Heegner point of discriminant -7 on elliptic curve of conductor 389
            sage: P.numerical_approx()
            (0.000000000000000 : 1.00000000000000 : 0.000000000000000)

        However, Heegner points of bigger conductor are often nonzero::

            sage: E = EllipticCurve('389a'); P = E.heegner_point(-7, 5); P
            Heegner point of discriminant -7 and conductor 5 on elliptic curve
             of conductor 389
            sage: numerical_approx(P)
            (0.675507556926807 + 0.344749649302635*I
             : -0.377142931401887 + 0.843366227137146*I : 1.00000000000000)
            sage: P.numerical_approx()
            (0.6755075569268... + 0.3447496493026...*I
             : -0.3771429314018... + 0.8433662271371...*I : 1.00000000000000)
            sage: E.heegner_point(-7, 11).numerical_approx()
            (0.1795583794118... + 0.02035501750912...*I
             : -0.5573941377055... + 0.2738940831635...*I : 1.00000000000000)
            sage: E.heegner_point(-7, 13).numerical_approx()
            (1.034302915374... - 3.302744319777...*I
             : 1.323937875767... + 6.908264226850...*I : 1.00000000000000)

        We find (probably) the defining polynomial of the
        `x`-coordinate of `P`, which defines a class field.  The shape of
        the discriminant below is strong confirmation -- but not proof
        -- that this polynomial is correct::

            sage: f = P.numerical_approx(70)[0].algdep(6); f
            1225*x^6 + 1750*x^5 - 21675*x^4 - 380*x^3 + 110180*x^2 - 129720*x + 48771
            sage: f.discriminant().factor()
            2^6 * 3^2 * 5^11 * 7^4 * 13^2 * 19^6 * 199^2 * 719^2 * 26161^2
        """
        tau = ComplexField(prec)(self.tau())
        E = self.curve()
        return E.modular_parametrization()(tau)

    def tau(self):
        r"""
        Return `\tau` in the upper half plane that maps via the
        modular parametrization to this Heegner point.

        EXAMPLES::

            sage: E = EllipticCurve('389a'); P = E.heegner_point(-7, 5)
            sage: P.tau()
            5/778*sqrt_minus_7 - 147/778
        """
        return self.heegner_point_on_X0N().tau()

    @cached_method
    def x_poly_exact(self, prec=53, algorithm='lll'):
        """
        Return irreducible polynomial over the rational numbers
        satisfied by the `x` coordinate of this Heegner point.  A
        :exc:`ValueError` is raised if the precision is clearly insignificant
        to define a point on the curve.

        .. WARNING::

            It is in theory possible for this function to not raise a
            :exc:`ValueError`, find a polynomial, but via some very unlikely
            coincidence that point is not actually this Heegner point.

        INPUT:

        - ``prec`` -- integer (default: 53)

        - ``algorithm`` -- 'conjugates' or 'lll' (default); if
          'conjugates', compute numerically all the
          conjugates ``y[i]`` of the Heegner point and construct
          the characteristic polynomial as the product
          `f(X)=(X-y[i])`.  If 'lll', compute only one of the
          conjugates ``y[0]``, then uses the LLL algorithm to
          guess `f(X)`.

        EXAMPLES:

        We compute some `x`-coordinate polynomials of some conductor 1
        Heegner points::

            sage: E = EllipticCurve('37a')
            sage: v = E.heegner_discriminants_list(10)
            sage: [E.heegner_point(D).x_poly_exact() for D in v]
            [x, x, x^2 + 2, x^5 - x^4 + x^3 + x^2 - 2*x + 1, x - 6,
             x^7 - 2*x^6 + 9*x^5 - 10*x^4 - x^3 + 8*x^2 - 5*x + 1,
             x^3 + 5*x^2 + 10*x + 4, x^4 - 10*x^3 + 10*x^2 + 12*x - 12,
             x^8 - 5*x^7 + 7*x^6 + 13*x^5 - 10*x^4 - 4*x^3 + x^2 - 5*x + 7,
             x^6 - 2*x^5 + 11*x^4 - 24*x^3 + 30*x^2 - 16*x + 4]

        We compute `x`-coordinate polynomials for some Heegner points
        of conductor bigger than 1 on a rank 2 curve::

            sage: E = EllipticCurve('389a'); P = E.heegner_point(-7, 5); P
            Heegner point of discriminant -7 and conductor 5
             on elliptic curve of conductor 389
            sage: P.x_poly_exact()
            Traceback (most recent call last):
            ...
            ValueError: insufficient precision to determine Heegner point
            (fails discriminant test)
            sage: P.x_poly_exact(120)
            x^6 + 10/7*x^5 - 867/49*x^4 - 76/245*x^3 + 3148/35*x^2 - 25944/245*x + 48771/1225
            sage: E.heegner_point(-7, 11).x_poly_exact(500)
            x^10 + 282527/52441*x^9 + 27049007420/2750058481*x^8 - 22058564794/2750058481*x^7
             - 140054237301/2750058481*x^6 + 696429998952/30250643291*x^5
             + 2791387923058/30250643291*x^4 - 3148473886134/30250643291*x^3
             + 1359454055022/30250643291*x^2 - 250620385365/30250643291*x
             + 181599685425/332757076201

        Here we compute a Heegner point of conductor 5 on a rank 3 curve::

            sage: E = EllipticCurve('5077a'); P = E.heegner_point(-7,5); P
            Heegner point of discriminant -7 and conductor 5
             on elliptic curve of conductor 5077
            sage: P.x_poly_exact(500)
            x^6 + 1108754853727159228/72351048803252547*x^5 + 88875505551184048168/1953478317687818769*x^4 - 2216200271166098662132/3255797196146364615*x^3 + 14941627504168839449851/9767391588439093845*x^2 - 3456417460183342963918/3255797196146364615*x + 1306572835857500500459/5426328660243941025

        See :issue:`34121`::

            sage: E = EllipticCurve('11a1')
            sage: P = E.heegner_point(-7)
            sage: PE = P.point_exact()
            sage: PE
            (a : -4*a + 3 : 1)
            sage: all(c.parent().disc() == -7 for c in PE)
            True
        """
        L = self.ring_class_field()
        n = L.absolute_degree()

        if algorithm == 'lll':
            P = self.numerical_approx(prec)
            f = P[0].algdep(n)
            if f.is_irreducible() and self._check_poly_discriminant(f):
                return f.monic()
            else:
                raise ValueError("insufficient precision to determine Heegner point (fails discriminant test)")
        else:
            raise NotImplementedError("'lll' is the only algorithm implemented for Heegner points")

    def _check_poly_discriminant(self, f):
        """
        Return ``True`` if the prime to `Dc` part of the discriminant of
        each factor of the polynomial `f` is plus or minus a square.
        This is used for verifying that a polynomial is likely to
        define a subfield of a specific ring class field.

        INPUT:

        - ``f`` -- a polynomial

        EXAMPLES::

            sage: E = EllipticCurve('389a'); P = E.heegner_point(-7, 5); P
            Heegner point of discriminant -7 and conductor 5
             on elliptic curve of conductor 389
            sage: R.<x> = QQ[]
            sage: P._check_poly_discriminant(x^2 - 5)
            True
            sage: P._check_poly_discriminant(x^2 - 19)
            False
            sage: P._check_poly_discriminant((x^2 - 19)*(x^2-5))
            False
        """
        if f.is_irreducible():
            disc = f.discriminant()
            D, c = self.discriminant(), self.conductor()
            for p in D.prime_divisors() + c.prime_divisors():
                disc = disc // (p**disc.valuation(p))
            if disc < 0:
                disc = -disc
            return disc.is_square()

        return all(self._check_poly_discriminant(g) for g,_ in f.factor())

    def point_exact(self, prec=53, algorithm='lll', var='a', optimize=False):
        """
        Return exact point on the elliptic curve over a number field
        defined by computing this Heegner point to the given number of
        bits of precision. A :exc:`ValueError` is raised if the precision
        is clearly insignificant to define a point on the curve.

        .. WARNING::

            It is in theory possible for this function to not raise a
            :exc:`ValueError`, find a point on the curve, but via some very
            unlikely coincidence that point is not actually this Heegner point.

        .. WARNING::

            Currently we make an arbitrary choice of `y`-coordinate for
            the lift of the `x`-coordinate.

        INPUT:

        - ``prec`` -- integer (default: 53)

        - ``algorithm`` -- see the description of the algorithm
          parameter for the ``x_poly_exact`` method

        - ``var`` -- string (default: ``'a'``)

        - ``optimize`` -- boolean (default: ``False``); if ``True``, try to
          optimize defining polynomial for the number field that
          the point is defined over.  Off by default, since this
          can be very expensive.

        EXAMPLES::

            sage: E = EllipticCurve('389a'); P = E.heegner_point(-7, 5); P
            Heegner point of discriminant -7 and conductor 5
             on elliptic curve of conductor 389
            sage: z = P.point_exact(200, optimize=True)
            sage: z[1].charpoly()
            x^12 + 6*x^11 + 90089/1715*x^10 + 71224/343*x^9 + 52563964/588245*x^8 - 483814934/588245*x^7 - 156744579/16807*x^6 - 2041518032/84035*x^5 + 1259355443184/14706125*x^4 + 3094420220918/14706125*x^3 + 123060442043827/367653125*x^2 + 82963044474852/367653125*x + 211679465261391/1838265625
            sage: f = P.numerical_approx(500)[1].algdep(12); f / f.leading_coefficient()
            x^12 + 6*x^11 + 90089/1715*x^10 + 71224/343*x^9 + 52563964/588245*x^8 - 483814934/588245*x^7 - 156744579/16807*x^6 - 2041518032/84035*x^5 + 1259355443184/14706125*x^4 + 3094420220918/14706125*x^3 + 123060442043827/367653125*x^2 + 82963044474852/367653125*x + 211679465261391/1838265625

            sage: E = EllipticCurve('5077a')
            sage: P = E.heegner_point(-7)
            sage: P.point_exact(prec=100)
            (0 : 1 : 0)
        """
        E = self.__E
        if self.numerical_approx(prec)[-1] == 0:
            return E(0)
        f = self.x_poly_exact(prec, algorithm=algorithm)
        if f.degree() == 1:
            v = E.lift_x(-f[0], all=True)
            if v:
                return v[0]

        g, d = make_monic(f)
        K = NumberField(g, var)
        x = K.gen() / d
        if optimize:
            KO, from_KO, to_KO = K.optimized_representation()
            K = KO
            x = to_KO(x)
            if K.degree() < 2 * self.ring_class_field().degree_over_K():
                M = QuadraticField(self.discriminant(),'b')
                KD = K.composite_fields(M, names='a')[0]
                phi = K.embeddings(KD)[0]
                x = phi(x)
                K = KD.change_names(names=var)
            x = K.structure()[1](x)
        a1, a2, a3, a4, a6 = E.a_invariants()
        R = K['Y']
        Y = R.gen()
        g = Y**2 + a1*x*Y + a3*Y - (x**3 + a2*x**2 + a4*x + a6)
        F = g.factor()   # this takes a long time
        if len(F) == 1 and F[0][0] == 2:
            # reducible -- 1 factor squared
            y = F[0][0]
            L = K
        elif len(F) == 2:
            # reducible -- 2 factors
            y0 = -F[0][0][0]
            # y1 = -F[1][0][0]
            # Figure out which of y0 or y1 is right by
            # P = self.numerical_approx(prec)
            # TODO: finish this -- have to do some thing numerical
            y = y0
            L = K
        else:
            # TODO -- is there an issue with choice of root?
            # irreducible
            gg, dd = make_monic(g)
            M = K.extension(gg, names='b')
            y = M.gen()/dd
            x = M(x)
            L = M.absolute_field(names=var)
            phi = L.structure()[1]
            x = phi(x)
            y = phi(y)

        EL = E.change_ring(L)
        P = EL.point((x,y,L(1)), check=False)
        return P

    @cached_method
    def conjugates_over_K(self):
        r"""
        Return the `Gal(K_c/K)` conjugates of this Heegner point.

        EXAMPLES::

            sage: E = EllipticCurve('77a')
            sage: y = E.heegner_point(-52,5); y
            Heegner point of discriminant -52 and conductor 5
             on elliptic curve of conductor 77
            sage: print([z.quadratic_form() for z in y.conjugates_over_K()])
            [77*x^2 + 52*x*y + 13*y^2, 154*x^2 + 206*x*y + 71*y^2, 539*x^2 + 822*x*y + 314*y^2,
             847*x^2 + 1284*x*y + 487*y^2, 1001*x^2 + 52*x*y + y^2, 1078*x^2 + 822*x*y + 157*y^2,
             1309*x^2 + 360*x*y + 25*y^2, 1309*x^2 + 2054*x*y + 806*y^2,
             1463*x^2 + 976*x*y + 163*y^2, 2233*x^2 + 2824*x*y + 893*y^2,
             2387*x^2 + 2054*x*y + 442*y^2, 3619*x^2 + 3286*x*y + 746*y^2]
            sage: y.quadratic_form()
            77*x^2 + 52*x*y + 13*y^2
        """
        H = heegner_points(self.level(), self.discriminant(), self.conductor())
        E = self.curve()
        beta = self.quadratic_form()[1]
        return tuple([z.map_to_curve(E) for z in H.points(beta)])

    def _numerical_approx_conjugates_over_QQ(self, prec=53):
        """
        Return a list v of the numerical approximations to precision
        prec of the conjugates of this Heegner point, and their
        complex conjugates.

        INPUT:

        - ``prec`` -- positive integer (default: 53)

        EXAMPLES::

            sage: E = EllipticCurve('37a')
            sage: y = E.heegner_point(-7,3); y
            Heegner point of discriminant -7 and conductor 3 on elliptic curve of conductor 37
            sage: y._numerical_approx_conjugates_over_QQ()
            [(-1.89564392373896 - 0.444771808762067*I : -1.50000000000000 + 2.13102976222246*I : 1.00000000000000), ...]
            sage: y._numerical_approx_conjugates_over_QQ(prec=10)
            [(-1.9 - 0.44*I : -1.5 + 2.1*I : 1.0), ...
             (1.4 + 0.0024*I : -1.7 - 0.0046*I : 1.0)]
        """
        v = []
        for z in self.conjugates_over_K():
            m = z.numerical_approx(prec)
            v.append(m)
            v.append(m.curve().point([w.conjugate() for w in m], check=False))
        v.sort()
        return v

    def _numerical_approx_xy_poly(self, prec=53):
        r"""
        Return polynomials with real floating point coefficients got
        by taking the real part of the product of `X - \alpha` over
        the numerical approximations `\alpha` to the conjugates of
        this Heegner point.  The first polynomial runs through the
        `x`-coordinates and the second through the `y`-coordinates.

        INPUT:

        - ``prec`` -- positive integer (default: 53)

        OUTPUT: 2-tuple of polynomials with floating point coefficients

        EXAMPLES::

            sage: E = EllipticCurve('37a')
            sage: y = E.heegner_point(-7,3); y
            Heegner point of discriminant -7 and conductor 3 on elliptic curve of conductor 37
            sage: y._numerical_approx_xy_poly()  # rel tol 1e-14
            (X^8 + 6.00000000000000*X^7 + 8.99999999999998*X^6 - 12.0000000000000*X^5 - 42.0000000000000*X^4 - 17.9999999999999*X^3 + 36.0000000000001*X^2 + 35.9999999999999*X + 8.99999999999995,
             X^8 + 12.0000000000000*X^7 + 72.0000000000000*X^6 + 270.000000000000*X^5 + 678.000000000001*X^4 + 1152.00000000000*X^3 + 1269.00000000000*X^2 + 810.000000000002*X + 225.000000000001)
        """
        v = self._numerical_approx_conjugates_over_QQ(prec)
        R = ComplexField(prec)['X']
        S = RealField(prec)['X']
        X = R.gen()
        fx = prod(X-a[0] for a in v)
        fx = S([b.real() for b in fx])
        fy = prod(X-c[1] for c in v)
        fy = S([d.real() for d in fy])
        return fx, fy

    def _xy_poly_nearby(self, prec=53, max_error=10**(-10)):
        """
        Return polynomials with rational coefficients that for sufficiently
        tight bounds are the characteristic polynomial of the x and y
        coordinate of this Heegner point.

        INPUT:

        - ``prec`` -- positive integer (default: 53)

        - ``max_error`` -- very small floating point number

        OUTPUT: 2-tuple of polynomials with rational coefficients

        EXAMPLES::

            sage: E = EllipticCurve('37a')
            sage: y = E.heegner_point(-7,3); y
            Heegner point of discriminant -7 and conductor 3 on elliptic curve of conductor 37
            sage: y._xy_poly_nearby()
            [X^8 + 6*X^7 + 9*X^6 - 12*X^5 - 42*X^4 - 18*X^3 + 36*X^2 + 36*X + 9,
             X^8 + 12*X^7 + 72*X^6 + 270*X^5 + 678*X^4 + 1152*X^3 + 1269*X^2 + 810*X + 225]
        """
        v = self._numerical_approx_xy_poly(prec)
        return [nearby_rational_poly(g, max_error=max_error) for g in v]

    def _xy_poly_simplest(self, prec=53, prec2=None):
        """
        Return polynomials with rational coefficients that for
        sufficiently tight bounds are the characteristic polynomial of
        the x and y coordinate of this Heegner point.

        INPUT:

        - ``prec`` -- positive integer (default: 53)

        - ``prec2`` -- passed into simplest_rational_poly function

        EXAMPLES::

            sage: E = EllipticCurve('37a'); y = E.heegner_point(-7,3)
            sage: y._xy_poly_simplest()
            [X^8 + 6*X^7 + 9*X^6 - 12*X^5 - 42*X^4 - 18*X^3 + 36*X^2 + 36*X + 9,
             X^8 + 12*X^7 + 72*X^6 + 270*X^5 + 678*X^4 + 1152*X^3 + 1269*X^2 + 810*X + 225]
        """
        v = self._numerical_approx_xy_poly(prec)
        if prec2 is None:
            prec2 = max(2, prec - 20)
        return [simplest_rational_poly(g, prec2) for g in v]

    @cached_method
    def _square_roots_mod_2N_of_D_mod_4N(self):
        """
        Return the square roots of `D` modulo `4N` all reduced mod `2N`,
        without multiplicity.

        EXAMPLES::

            sage: E = EllipticCurve('37a'); P = E.heegner_point(-40); P
            Heegner point of discriminant -40 on elliptic curve of conductor 37
            sage: P._square_roots_mod_2N_of_D_mod_4N()
            [16, 58]
            sage: parent(P._square_roots_mod_2N_of_D_mod_4N()[0])
            Ring of integers modulo 74
        """
        N = self.__E.conductor()
        R = Integers(4*N)
        m = 2*N
        return sorted( {a % m for a in R(self.discriminant()).sqrt(all=True)} )

    def _trace_numerical_conductor_1(self, prec=53):
        """
        Return numerical approximation using ``prec`` terms of working
        precision to the trace down to the quadratic imaginary field
        `K` of this Heegner point.

        INPUT:

        - ``prec`` -- bits precision (default: 53)

        EXAMPLES::

            sage: E = EllipticCurve('57a1')
            sage: P = E.heegner_point(-8); P
            Heegner point of discriminant -8 on elliptic curve of conductor 57
            sage: P._trace_numerical_conductor_1() # approx. (1 : 0 : 1)
            (1.00000000000000 + ...e-16*I : ...e-16 - ...e-16*I : 1.00000000000000)
            sage: P = E(2,1) # a generator
            sage: E([1,0]).height()
            0.150298370947295
            sage: P.height()
            0.0375745927368238
            sage: E.heegner_index(-8)
            2.0000?
            sage: E.torsion_order()
            1
            sage: 2*P
            (1 : 0 : 1)
        """
        if self.conductor() != 1:
            raise ValueError("conductor must be 1")
        R, U = self._good_tau_representatives()
        E = self.__E
        phi = E.modular_parametrization()
        C = ComplexField(prec)
        F = E.change_ring(C)
        s = 0
        for u, weight in U:
            P = phi(C(self._qf_to_tau(u)))
            z = F.point(list(P),check=False)
            if abs(weight) == 2:
                t = F.point(z,check=False) + F.point(tuple([x.conjugate() for x in z]), check=False)
                if weight < 0:
                    s -= t
                else:
                    s += t
            else:
                if weight < 0:
                    s -= z
                else:
                    s += z
        return s

    @cached_method
    def _good_tau_representatives(self):
        """
        Return good upper half plane representatives for Heegner points.

        ALGORITHM: This is Algorithm 3.5 in Watkins's paper.

        EXAMPLES::

            sage: P = EllipticCurve('389a1').heegner_point(-7)
            sage: P._good_tau_representatives()
            ([(1, 1, 2)], [((389, 185, 22), 1)])
        """
        if self.conductor() != 1:
            raise NotImplementedError
        E = self.__E
        SDN = self._square_roots_mod_2N_of_D_mod_4N()
        beta = SDN[0]
        U = []
        R = []
        N = self.__E.conductor()
        D = self.discriminant()
        h = self.ring_class_field().degree_over_K()
        divs = D.gcd(N).divisors()
        a = 1
        while True:
            for b in SDN:
                b = b.lift()
                # todo (optimize) -- replace for over all s with for over solution
                # set that can be found quickly.
                y = ZZ((b*b - D)/(4*N))
                for s in Integers(a):
                    if N*s*s + b*s + y == 0:
                        s = s.lift()
                        f = (a*N, b+2*N*s, ZZ( ((b + 2*N*s)**2 - D)/(4*a*N)) )
                        for d in divs:
                            Q = d * prod(p**k for p,k in N.factor() if (b-beta) % (p**k) != 0)
                            g = self._qf_atkin_lehner_act(Q, f)
                            gbar = (ZZ(g[0]/N), -g[1], g[2]*N)
                            g = self._qf_reduce(g)
                            gbar = self._qf_reduce(gbar)
                            if g in R or gbar in R:
                                continue
                            R.append(g)
                            if g != gbar:
                                R.append(gbar)
                            epsilon_Q = prod([E.root_number(q) for q in Q.prime_divisors()])
                            if g == gbar:
                                # weight is epsilon_Q
                                weight = epsilon_Q
                            else:
                                # weight is 2*epsilon_Q
                                weight = 2*epsilon_Q
                            U.append((f,weight))
                            if len(R) == h:
                                return R, U
                            assert len(R) < h, "bug -- too many quadratic forms"
            a += 1

    def _qf_to_tau(self, f):
        r"""
        Function used internally that given a quadratic form
        `f=(A,B,C)`, return `\tau` in the upper half plane with
        `A\tau^2 + B \tau + C = 0`.  Here `A>0` and `\gcd(A,B,C)=1`.
        Also, `\tau` has discriminant `D=B^2-4AC`.  In fact, `\tau =
        (-B + \sqrt{D})/(2A)`.

        INPUT:

        - ``f`` -- binary quadratic form

        EXAMPLES::

            sage: P = EllipticCurve('57a1').heegner_point(-8)
            sage: R, U = P._good_tau_representatives()
            sage: f = U[0][0]; f
            (57, 26, 3)
            sage: P._qf_to_tau(f)
            1/114*sqrt_minus_8 - 13/57
        """
        c = self.conductor()
        A,B,_ = f
        alpha = c * self.quadratic_field().gen()   # this is sqrt(D) = sqrt(c^2*disc(K))
        return (-B + alpha)/(2*A)

    def _qf_from_tau(self, tau):
        r"""
        Return quadratic form associated to a given `\tau` in the upper
        half plane.

        INPUT:

        - `\tau` -- quadratic element of the upper half plane

        EXAMPLES::

            sage: P = EllipticCurve('57a1').heegner_point(-8)
            sage: R, U = P._good_tau_representatives()
            sage: f = U[0][0]; f
            (57, 26, 3)
            sage: tau = P._qf_to_tau(f); tau
            1/114*sqrt_minus_8 - 13/57
            sage: P._qf_from_tau(tau)
            (57, 26, 3)
        """
        g = tau.minpoly()
        g *= g.denominator()
        return (ZZ(g[2]), ZZ(g[1]), ZZ(g[0]))

    def _qf_atkin_lehner_act(self, Q, f):
        r"""
        Given a positive integer `Q` with `Q | N` and `\gcd(Q, N/Q) =
        1`, we compute the quadratic form corresponding to the image
        of the `tau` corresponding to `f` under the Atkin-Lehner
        operator `W_Q`.

        We do this by letting `u,v` be integers such that
        `u Q^2 - v N = Q`, and using that `W_Q` sends `\tau`
        to `( (u Q \tau + v) / (N \tau + Q) ) / Q`.

        INPUT:

        - ``Q`` -- integer that divides the level `N`

        - ``f`` -- quadratic form

        OUTPUT: quadratic form

        EXAMPLES::

            sage: P = EllipticCurve('57a1').heegner_point(-8)
            sage: R, U = P._good_tau_representatives()
            sage: f = U[0][0]; f
            (57, 26, 3)
            sage: P._qf_atkin_lehner_act(3, f)
            (1938, 1204, 187)
            sage: g = P._qf_atkin_lehner_act(19, f); g
            (114, -64, 9)
            sage: h = P._qf_atkin_lehner_act(19, g); h
            (7353, -4762, 771)
            sage: BinaryQF(f).reduced_form() == BinaryQF(h).reduced_form()
            True
        """
        N = self.__E.conductor()
        g, u, v = xgcd(Q*Q, -N)
        assert g == Q
        tau = self._qf_to_tau(f)
        tau2 = ((u*Q*tau + v) / (N*tau + Q))
        return self._qf_from_tau(tau2)

    def _qf_reduce(self, f):
        """
        Given a binary quadratic form `f` represented as a 3-tuple
        (A,B,C), return the reduced binary quadratic form equivalent
        to `f`, represented in the same way.

        EXAMPLES::

            sage: P = EllipticCurve('57a1').heegner_point(-8)
            sage: R, U = P._good_tau_representatives()
            sage: f = U[0][0]; f
            (57, 26, 3)
            sage: P._qf_reduce(f)
            (1, 0, 2)
        """
        return tuple(BinaryQF(f).reduced_form())

    def kolyvagin_cohomology_class(self, n=None):
        """
        Return the Kolyvagin class associated to this Heegner point.

        INPUT:

        - ``n`` -- positive integer that divides the gcd of `a_p`
          and `p+1` for all `p` dividing the conductor.  If `n` is
          ``None``, choose the largest valid `n`.

        EXAMPLES::

            sage: y = EllipticCurve('389a').heegner_point(-7,5)
            sage: y.kolyvagin_cohomology_class(3)
            Kolyvagin cohomology class c(5) in H^1(K,E[3])
        """
        return KolyvaginCohomologyClassEn(self.kolyvagin_point(), n)

#########################################################################################
# Kolyvagin Points P_c
#########################################################################################


class KolyvaginPoint(HeegnerPoint):
    """
    A Kolyvagin point.

    EXAMPLES:

    We create a few Kolyvagin points::

        sage: EllipticCurve('11a1').kolyvagin_point(-7)
        Kolyvagin point of discriminant -7 on elliptic curve of conductor 11
        sage: EllipticCurve('37a1').kolyvagin_point(-7)
        Kolyvagin point of discriminant -7 on elliptic curve of conductor 37
        sage: EllipticCurve('37a1').kolyvagin_point(-67)
        Kolyvagin point of discriminant -67 on elliptic curve of conductor 37
        sage: EllipticCurve('389a1').kolyvagin_point(-7, 5)
        Kolyvagin point of discriminant -7 and conductor 5
         on elliptic curve of conductor 389

    One can also associated a Kolyvagin point to a Heegner point::

        sage: y = EllipticCurve('37a1').heegner_point(-7); y
        Heegner point of discriminant -7 on elliptic curve of conductor 37
        sage: y.kolyvagin_point()
        Kolyvagin point of discriminant -7 on elliptic curve of conductor 37

    TESTS::

        sage: y = EllipticCurve('37a1').heegner_point(-7)
        sage: type(y)
        <class 'sage.schemes.elliptic_curves.heegner.HeegnerPointOnEllipticCurve'>
        sage: loads(dumps(y)) == y
        True
    """
    def __init__(self, heegner_point):
        """
        Create a Kolyvagin point.

        INPUT:

            - ``heegner_point`` -- a Heegner point on some elliptic curve

        EXAMPLES:

        We directly construct a Kolyvagin point from the KolyvaginPoint class::

            sage: y = EllipticCurve('37a1').heegner_point(-7)
            sage: sage.schemes.elliptic_curves.heegner.KolyvaginPoint(y)
            Kolyvagin point of discriminant -7 on elliptic curve of conductor 37
        """
        if not heegner_point.satisfies_kolyvagin_hypothesis():
            raise ValueError("Heegner point does not satisfy Kolyvagin hypothesis")
        self.__heegner_point = heegner_point
        HeegnerPoint.__init__(self, heegner_point.level(), heegner_point.discriminant(),
                              heegner_point.conductor())

    def satisfies_kolyvagin_hypothesis(self, n=None):
        r"""
        Return ``True`` if this Kolyvagin point satisfies the Heegner
        hypothesis for `n`, so that it defines a Galois equivariant
        element of `E(K_c)/n E(K_c)`.

        EXAMPLES::

            sage: y = EllipticCurve('389a').heegner_point(-7,5); P = y.kolyvagin_point()
            sage: P.kolyvagin_cohomology_class(3)
            Kolyvagin cohomology class c(5) in H^1(K,E[3])
            sage: P.satisfies_kolyvagin_hypothesis(3)
            True
            sage: P.satisfies_kolyvagin_hypothesis(5)
            False
            sage: P.satisfies_kolyvagin_hypothesis(7)
            False
            sage: P.satisfies_kolyvagin_hypothesis(11)
            False
        """
        return self.__heegner_point.satisfies_kolyvagin_hypothesis(n)

    def curve(self):
        r"""
        Return the elliptic curve over `\QQ` on which this Kolyvagin
        point sits.

        EXAMPLES::

            sage: E = EllipticCurve('37a1'); P = E.kolyvagin_point(-67, 3)
            sage: P.curve()
            Elliptic Curve defined by y^2 + y = x^3 - x over Rational Field
        """
        return self.__heegner_point.curve()

    def heegner_point(self):
        """
        This Kolyvagin point `P_c` is associated to some Heegner point
        `y_c` via Kolyvagin's construction.  This function returns that
        point `y_c`.

        EXAMPLES::

            sage: E = EllipticCurve('37a1')
            sage: P = E.kolyvagin_point(-67); P
            Kolyvagin point of discriminant -67 on elliptic curve of conductor 37
            sage: y = P.heegner_point(); y
            Heegner point of discriminant -67 on elliptic curve of conductor 37
            sage: y.kolyvagin_point() is P
            True
        """
        return self.__heegner_point

    def _repr_(self):
        """
        Return string representation of this Kolyvagin point.

        EXAMPLES::

            sage: E = EllipticCurve('37a1'); P = E.kolyvagin_point(-67,7); P._repr_()
            'Kolyvagin point of discriminant -67 and conductor 7 on elliptic curve of conductor 37'
        """
        s = repr(self.__heegner_point)
        return s.replace('Heegner','Kolyvagin')

    def index(self, *args, **kwds):
        """
        Return index of this Kolyvagin point in the full group of
        `K_c` rational points on `E`.

        When the conductor is 1, this is computed numerically using
        the Gross-Zagier formula and explicit point search, and it may
        be off by `2`. See the documentation for ``E.heegner_index``,
        where `E` is the curve attached to ``self``.

        EXAMPLES::

            sage: E = EllipticCurve('37a1'); P = E.kolyvagin_point(-67); P.index()
            6
        """
        if self.conductor() == 1:
            return self.__heegner_point._trace_index(*args, **kwds)
        raise NotImplementedError

    def numerical_approx(self, prec=53):
        """
        Return a numerical approximation to this Kolyvagin point using
        prec bits of working precision.

        INPUT:

        - ``prec`` -- precision in bits (default: 53)

        EXAMPLES::

            sage: P = EllipticCurve('37a1').kolyvagin_point(-7); P
            Kolyvagin point of discriminant -7 on elliptic curve of conductor 37
            sage: P.numerical_approx() # approx. (0 : 0 : 1)
            (...e-16 - ...e-16*I : ...e-16 + ...e-16*I : 1.00000000000000)
            sage: P.numerical_approx(100)[0].abs() < 2.0^-99
            True

            sage: P = EllipticCurve('389a1').kolyvagin_point(-7, 5); P
            Kolyvagin point of discriminant -7 and conductor 5
             on elliptic curve of conductor 389

        Numerical approximation is only implemented for points of conductor 1::

            sage: P.numerical_approx()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        if self.conductor() == 1:
            return self.__heegner_point._trace_numerical_conductor_1(prec)
        raise NotImplementedError

    def point_exact(self, prec=53):
        """
        INPUT:

        - ``prec`` -- precision in bits (default: 53)

        EXAMPLES:

        A rank 1 curve::

            sage: E = EllipticCurve('37a1'); P = E.kolyvagin_point(-67)
            sage: P.point_exact()
            (6 : -15 : 1)
            sage: P.point_exact(40)
            (6 : -15 : 1)
            sage: P.point_exact(20)
            Traceback (most recent call last):
            ...
            RuntimeError: insufficient precision to find exact point

        A rank 0 curve::

            sage: E = EllipticCurve('11a1'); P = E.kolyvagin_point(-7)
            sage: P.point_exact()
            (-1/2*sqrt_minus_7 + 1/2 : -2*sqrt_minus_7 - 2 : 1)

        A rank 2 curve::

            sage: E = EllipticCurve('389a1'); P = E.kolyvagin_point(-7)
            sage: P.point_exact()
            (0 : 1 : 0)
        """
        if self.conductor() == 1:
            # the result is a point defined over K in the conductor 1 case, which is easier.
            P = self.numerical_approx(prec)

            E = self.curve()
            if P[2] == 0:
                return E(0)

            if E.root_number() == -1:
                return self._recognize_point_over_QQ(P, 2*self.index())
            else:
                # root number +1.  We use algdep to recognize the x
                # coordinate, stick it in the appropriate quadratic
                # field, then make sure that we got the right
                # embedding, and if not fix things so we do.
                x = P[0]
                C = x.parent()
                f = x.algdep(2)
                K = self.quadratic_field()
                roots = [r[0] for r in f.roots(K)]
                if not roots:
                    raise RuntimeError("insufficient precision to find exact point")
                if len(roots) == 1:
                    X = roots[0]
                else:
                    d = [abs(C(r) - x) for r in roots]
                    if d[0] == d[1]:
                        raise RuntimeError("insufficient precision to distinguish roots")
                    if d[0] < d[1]:
                        X = roots[0]
                    else:
                        X = roots[1]
                F = E.change_ring(K)
                Q = F.lift_x(X, all=True)
                if len(Q) == 1:
                    return Q[0]
                if not Q:
                    raise RuntimeError("insufficient precision")
                y = P[1]
                d = [abs(C(r[1])-y) for r in Q]
                if d[0] == d[1]:
                    raise RuntimeError("insufficient precision to distinguish roots")
                if d[0] < d[1]:
                    return Q[0]
                else:
                    return Q[1]

        else:
            raise NotImplementedError

    def plot(self, prec=53, *args, **kwds):
        r"""
        Plot a Kolyvagin point `P_1` if it is defined over the
        rational numbers.

        EXAMPLES::

            sage: E = EllipticCurve('37a'); P = E.heegner_point(-11).kolyvagin_point()
            sage: P.plot(prec=30, pointsize=50, rgbcolor='red') + E.plot()              # needs sage.plot
            Graphics object consisting of 3 graphics primitives
        """
        if self.conductor() != 1:
            raise NotImplementedError

        E = self.curve()
        if E.root_number() == -1:
            P = self.numerical_approx(prec=prec)
            from sage.plot.all import point, Graphics
            if not P:
                # point at infinity
                return Graphics()
            return point((P[0].real(), P[1].real()), *args, **kwds)
        else:
            raise NotImplementedError

    @cached_method
    def trace_to_real_numerical(self, prec=53):
        """
        Return the trace of this Kolyvagin point down to the real
        numbers, computed numerically using prec bits of working
        precision.

        EXAMPLES::

            sage: E = EllipticCurve('37a1'); P = E.kolyvagin_point(-67)
            sage: PP = P.numerical_approx()
            sage: [c.real() for c in PP]
            [6.00000000000000, -15.0000000000000, 1.00000000000000]
            sage: all(c.imag().abs() < 1e-14 for c in PP)
            True
            sage: P.trace_to_real_numerical()
            (1.61355529131986 : -2.18446840788880 : 1.00000000000000)
            sage: P.trace_to_real_numerical(prec=80)  # abs tol 1e-21
            (1.6135552913198573127230 : -2.1844684078888023289187 : 1.0000000000000000000000)
        """
        # Compute numerical approximation of P in E(K).
        P = self.numerical_approx(prec=prec)
        # Trace this numerical approximation down to E(Q) (numerically).
        E = P.curve()
        if self.curve().root_number() == -1:
            R = 2*P
        else:
            R = P + E.point([x.conjugate() for x in P],check=False)
        F = self.curve().change_ring(RealField(prec))
        return F.point([x.real() for x in R], check=False)

    @cached_method
    def _trace_exact_conductor_1(self, prec=53):
        r"""
        Return the trace from `K` to `\QQ` of this Kolyvagin point in
        the case of conductor 1, computed using prec bits of
        precision, then approximated using some algorithm (e.g.,
        continued fractions).  If the precision is not enough to
        determine a point on the curve, then a :exc:`RuntimeError` is raised.
        Even if the precision determines a point, there is no guarantee
        that it is correct.

        EXAMPLES:

        A Kolyvagin point on a rank 1 curve::

            sage: E = EllipticCurve('37a1'); P = E.kolyvagin_point(-67)
            sage: P.trace_to_real_numerical()
            (1.61355529131986 : -2.18446840788880 : 1.00000000000000)
            sage: P._trace_exact_conductor_1()  # the actual point we're reducing
            (1357/841 : -53277/24389 : 1)
            sage: (P._trace_exact_conductor_1().height() / E.regulator()).sqrt()
            12.0000000000000
        """
        if not self.conductor() == 1:
            raise ValueError("the conductor must be 1")

        P = self.trace_to_real_numerical(prec)
        return self._recognize_point_over_QQ(P, 2*self.index())

    def _recognize_point_over_QQ(self, P, n):
        r"""
        Used internally when computing an exact point on an elliptic curve.

        INPUT:

        - ``P`` -- numerical approximation for a point on `E`

        - ``n`` -- upper bound on divisibility index of `P` in group `E(\QQ)`

        EXAMPLES::

            sage: E = EllipticCurve('43a'); P = E.heegner_point(-20).kolyvagin_point()
            sage: PP = P.numerical_approx(); PP
            (0.000000000000000 : -1.00000000000000 : 1.00000000000000)
            sage: P._recognize_point_over_QQ(PP, 4)
            (0 : -1 : 1)
        """
        # Here is where we *should* implement the "method of Cremona
        # etc" mentioned in Watkins' article... which involves local
        # heights.
        E = self.curve()  # over Q
        v = sum([list(n*w) for w in E.gens()] + [list(w) for w in E.torsion_points()], [])
        # note -- we do not claim to prove anything, so making up a factor of 100 is fine.
        max_denominator = 100*max([z.denominator() for z in v])
        try:
            # the coercion below also checks if point is on elliptic curve
            return E([x.real().nearby_rational(max_denominator=max_denominator) for x in P])
        except TypeError:
            raise RuntimeError("insufficient precision to find exact point")

    def mod(self, p, prec=53):
        r"""
        Return the trace of the reduction `Q` modulo a prime over `p` of this
        Kolyvagin point as an element of `E(\GF{p})`, where
        `p` is any prime that is inert in `K` that is coprime to `NDc`.

        The point `Q` is only well defined up to an element of
        `(p+1) E(\GF{p})`, i.e., it gives a well defined element
        of the abelian group `E(\GF{p}) / (p+1) E(\GF{p})`.

        See [St2011b]_, Proposition 5.4 for a proof of the above
        well-definedness assertion.

        EXAMPLES:

        A Kolyvagin point on a rank 1 curve::

            sage: E = EllipticCurve('37a1'); P = E.kolyvagin_point(-67)
            sage: P.mod(2)
            (1 : 1 : 1)
            sage: P.mod(3)
            (1 : 0 : 1)
            sage: P.mod(5)
            (2 : 2 : 1)
            sage: P.mod(7)
            (6 : 0 : 1)
            sage: P.trace_to_real_numerical()
            (1.61355529131986 : -2.18446840788880 : 1.00000000000000)
            sage: P._trace_exact_conductor_1()  # the actual point we're reducing
            (1357/841 : -53277/24389 : 1)
            sage: (P._trace_exact_conductor_1().height() / E.regulator()).sqrt()
            12.0000000000000

        Here the Kolyvagin point is a torsion point (since `E` has
        rank 1), and we reduce it modulo several primes.::

            sage: E = EllipticCurve('11a1'); P = E.kolyvagin_point(-7)
            sage: P.mod(3,70)  # long time (4s on sage.math, 2013)
            (1 : 2 : 1)
            sage: P.mod(5,70)
            (1 : 4 : 1)
            sage: P.mod(7,70)
            Traceback (most recent call last):
            ...
            ValueError: p must be coprime to conductors and discriminant
            sage: P.mod(11,70)
            Traceback (most recent call last):
            ...
            ValueError: p must be coprime to conductors and discriminant
            sage: P.mod(13,70)
            (3 : 4 : 1)
        """
        # check preconditions
        p = ZZ(p)
        if not p.is_prime():
            raise ValueError("p must be prime")
        E = self.curve()
        D = self.discriminant()
        if (E.conductor() * D * self.conductor()) % p == 0:
            raise ValueError("p must be coprime to conductors and discriminant")
        K = self.heegner_point().quadratic_field()
        if len(K.factor(p)) != 1:
            raise ValueError("p must be inert")

        # do actual calculation
        if self.conductor() == 1:
            P = self._trace_exact_conductor_1(prec=prec)
            return E.change_ring(GF(p))(P)
        raise NotImplementedError

##     def congruent_rational_point(self, n, prec=53):
##         r"""
##         Let `P` be this Kolyvagin point.  Determine whether there is a
##         point `z` in `E(\QQ)` such that `z - P \in n E(K_c)`, where `K_c`
##         is the ring class field over which this Kolyvagin point is defined.
##         If `z` exists return `z`.  Otherwise return None.
##
##         INPUT:
##
##            - ``n`` -- positive integer
##
##            - ``prec`` -- positive integer (default: 53)
##
##
##         EXAMPLES::
##
##         """
##         raise NotImplementedError

    def kolyvagin_cohomology_class(self, n=None):
        """
        INPUT:

        - ``n`` -- positive integer that divides the gcd of `a_p`
          and `p+1` for all `p` dividing the conductor.  If `n` is
          ``None``, choose the largest valid `n`.

        EXAMPLES::

            sage: y = EllipticCurve('389a').heegner_point(-7, 5)
            sage: P = y.kolyvagin_point()
            sage: P.kolyvagin_cohomology_class(3)
            Kolyvagin cohomology class c(5) in H^1(K,E[3])

            sage: y = EllipticCurve('37a').heegner_point(-7, 5).kolyvagin_point()
            sage: y.kolyvagin_cohomology_class()
            Kolyvagin cohomology class c(5) in H^1(K,E[2])
        """
        return KolyvaginCohomologyClassEn(self, n)


class KolyvaginCohomologyClass(SageObject):
    """
    A Kolyvagin cohomology class in `H^1(K,E[n])` or `H^1(K,E)[n]`
    attached to a Heegner point.

    EXAMPLES::

        sage: y = EllipticCurve('37a').heegner_point(-7)
        sage: c = y.kolyvagin_cohomology_class(3); c
        Kolyvagin cohomology class c(1) in H^1(K,E[3])
        sage: type(c)
        <class 'sage.schemes.elliptic_curves.heegner.KolyvaginCohomologyClassEn'>
        sage: loads(dumps(c)) == c
        True
        sage: y.kolyvagin_cohomology_class(5)
        Kolyvagin cohomology class c(1) in H^1(K,E[5])
    """
    def __init__(self, kolyvagin_point, n):
        """

        EXAMPLES::

            sage: y = EllipticCurve('389a').heegner_point(-7, 5)
            sage: y.kolyvagin_cohomology_class(3)
            Kolyvagin cohomology class c(5) in H^1(K,E[3])
        """
        if n is None:
            c = kolyvagin_point.conductor()
            E = kolyvagin_point.curve()
            n = gcd([(p+1).gcd(E.ap(p)) for p in c.prime_divisors()])

        if not kolyvagin_point.satisfies_kolyvagin_hypothesis(n):
            raise ValueError("Kolyvagin point does not satisfy Kolyvagin hypothesis for %s" % n)
        self.__kolyvagin_point = kolyvagin_point
        self.__n = n

    def __eq__(self, other):
        """
        EXAMPLES::

            sage: y = EllipticCurve('37a').heegner_point(-7)
            sage: c = y.kolyvagin_cohomology_class(3)
            sage: c == y.kolyvagin_cohomology_class(3)
            True
            sage: c == y.kolyvagin_cohomology_class(5)
            False

        This does not mean that c is nonzero (!) -- it just means c is not the number 0::

            sage: c == 0
            False
        """
        return isinstance(other, KolyvaginCohomologyClass) and \
               self.__kolyvagin_point == other.__kolyvagin_point and \
               self.__n == other.__n

    def __ne__(self, other):
        """
        EXAMPLES::

            sage: y = EllipticCurve('37a').heegner_point(-7)
            sage: c = y.kolyvagin_cohomology_class(3)
            sage: c != y.kolyvagin_cohomology_class(3)
            False
            sage: c != y.kolyvagin_cohomology_class(5)
            True
        """
        return not (self == other)

    def n(self):
        """
        Return the integer `n` so that this is a cohomology class in
        `H^1(K,E[n])` or `H^1(K,E)[n]`.

        EXAMPLES::

            sage: y = EllipticCurve('37a').heegner_point(-7)
            sage: t = y.kolyvagin_cohomology_class(3); t
            Kolyvagin cohomology class c(1) in H^1(K,E[3])
            sage: t.n()
            3
        """
        return self.__n

    def conductor(self):
        r"""
        Return the integer `c` such that this cohomology class is associated
        to the Heegner point `y_c`.

        EXAMPLES::

            sage: y = EllipticCurve('37a').heegner_point(-7, 5)
            sage: t = y.kolyvagin_cohomology_class()
            sage: t.conductor()
            5
        """
        return self.__kolyvagin_point.conductor()

    def kolyvagin_point(self):
        """
        Return the Kolyvagin point `P_c` to which this cohomology
        class is associated.

        EXAMPLES::

            sage: y = EllipticCurve('37a').heegner_point(-7, 5)
            sage: t = y.kolyvagin_cohomology_class()
            sage: t.kolyvagin_point()
            Kolyvagin point of discriminant -7 and conductor 5
             on elliptic curve of conductor 37
        """
        return self.__kolyvagin_point

    def heegner_point(self):
        """
        Return the Heegner point `y_c` to which this cohomology class
        is associated.

        EXAMPLES::

            sage: y = EllipticCurve('37a').heegner_point(-7, 5)
            sage: t = y.kolyvagin_cohomology_class()
            sage: t.heegner_point()
            Heegner point of discriminant -7 and conductor 5
             on elliptic curve of conductor 37
        """
        return self.__kolyvagin_point.heegner_point()


class KolyvaginCohomologyClassEn(KolyvaginCohomologyClass):
    def _repr_(self):
        """

        EXAMPLES::

            sage: y = EllipticCurve('37a').heegner_point(-7, 5)
            sage: t = y.kolyvagin_cohomology_class()
            sage: t._repr_()
            'Kolyvagin cohomology class c(5) in H^1(K,E[2])'
        """
        return "Kolyvagin cohomology class c(%s) in H^1(K,E[%s])" % (
            self.conductor(), self.n())


#############################################################################
# Reduction of Heegner points using Quaternion Algebras
#
# This section contains implementations of algorithms for computing
# information about reduction modulo primes of Heegner points using
# quaternion algebras.  Some of this code could later be moved to the
# quaternion algebras code, but it is too immature and not general
# enough at present for that.
#############################################################################

class HeegnerQuatAlg(SageObject):
    r"""
    Heegner points viewed as supersingular points on the modular curve
    `X_0(N)/\mathbf{F}_{\ell}`.

    EXAMPLES::

        sage: H = heegner_points(11).reduce_mod(13); H
        Heegner points on X_0(11) over F_13
        sage: type(H)
        <class 'sage.schemes.elliptic_curves.heegner.HeegnerQuatAlg'>
        sage: loads(dumps(H)) == H
        True
    """
    def __init__(self, level, ell):
        r"""
        INPUT:

        - ``level`` -- the level (a positive integer)

        - `\ell` -- the characteristic, a prime coprime to the level

        EXAMPLES::

            sage: sage.schemes.elliptic_curves.heegner.HeegnerQuatAlg(11, 13)
            Heegner points on X_0(11) over F_13
        """
        level = ZZ(level)
        ell = ZZ(ell)
        if not ell.is_prime():
            raise ValueError("ell must be prime")
        if level.gcd(ell) != 1:
            raise ValueError("level and ell must be coprime")
        self.__level = level
        self.__ell = ell

    def __eq__(self, other):
        """
        EXAMPLES::

            sage: H = heegner_points(11).reduce_mod(3)
            sage: H == heegner_points(11).reduce_mod(3)
            True
            sage: H == heegner_points(11).reduce_mod(5)
            False
            sage: H == 0
            False
        """
        return isinstance(other, HeegnerQuatAlg) and self.__level == other.__level \
               and self.__ell == other.__ell

    def __ne__(self, other):
        """
        EXAMPLES::

            sage: H = heegner_points(11).reduce_mod(3)
            sage: H != heegner_points(11).reduce_mod(3)
            False
            sage: H != heegner_points(11).reduce_mod(5)
            True
            sage: H != 0
            True
        """
        return not (self == other)

    def _repr_(self):
        """
        Return string representation.

        EXAMPLES::

            sage: heegner_points(11).reduce_mod(13)._repr_()
            'Heegner points on X_0(11) over F_13'
        """
        return "Heegner points on X_0(%s) over F_%s" % (
            self.__level, self.__ell)

    def level(self):
        """
        Return the level.

        EXAMPLES::

            sage: heegner_points(11).reduce_mod(3).level()
            11
        """
        return self.__level

    def ell(self):
        r"""
        Return the prime `\ell` modulo which we are working.

        EXAMPLES::

            sage: heegner_points(11).reduce_mod(3).ell()
            3
        """
        return self.__ell

    def satisfies_heegner_hypothesis(self, D, c=ZZ(1)):
        r"""
        The fundamental discriminant `D` must be coprime to `N\ell`,
        and must define a quadratic imaginary field `K` in which `\ell`
        is inert.  Also, all primes dividing `N` must split in `K`,
        and `c` must be squarefree and coprime to `ND\ell`.

        INPUT:

        - ``D`` -- negative integer

        - ``c`` -- positive integer (default: 1)

        OUTPUT: boolean

        EXAMPLES::

            sage: H = heegner_points(11).reduce_mod(7)
            sage: H.satisfies_heegner_hypothesis(-5)
            False
            sage: H.satisfies_heegner_hypothesis(-7)
            False
            sage: H.satisfies_heegner_hypothesis(-8)
            True
            sage: [D for D in [-1,-2..-100] if H.satisfies_heegner_hypothesis(D)]
            [-8, -39, -43, -51, -79, -95]
        """
        D = ZZ(D)
        c = ZZ(c)
        if (c * D).gcd(self.__level * self.__ell) != 1 or c.gcd(D) != 1:
            return False
        if not satisfies_weak_heegner_hypothesis(self.__level, D):
            return False
        if not is_inert(D, self.__ell):
            return False
        return True

    def heegner_discriminants(self, n=5):
        r"""
        Return the first `n` negative fundamental discriminants
        coprime to `N\ell` such that `\ell` is inert in the
        corresponding quadratic imaginary field and that field
        satisfies the Heegner hypothesis, and `N` is the level.

        INPUT:

        - ``n`` -- positive integer (default: 5)

        OUTPUT: list

        EXAMPLES::

            sage: H = heegner_points(11).reduce_mod(3)
            sage: H.heegner_discriminants()
            [-7, -19, -40, -43, -52]
            sage: H.heegner_discriminants(10)
            [-7, -19, -40, -43, -52, -79, -127, -139, -151, -184]
        """
        v = []
        D = ZZ(-5)
        while len(v) < n:
            if self.satisfies_heegner_hypothesis(D):
                v.append(D)
            D -= 1
        return v

    def heegner_conductors(self, D, n=5):
        r"""
        Return the first `n` negative fundamental discriminants
        coprime to `N\ell` such that `\ell` is inert in the
        corresponding quadratic imaginary field and that field
        satisfies the Heegner hypothesis.

        INPUT:

        - ``D`` -- negative integer; a fundamental Heegner discriminant

        - ``n`` -- positive integer (default: 5)

        OUTPUT: list

        EXAMPLES::

            sage: H = heegner_points(11).reduce_mod(3)
            sage: H.heegner_conductors(-7)
            [1, 2, 4, 5, 8]
            sage: H.heegner_conductors(-7, 10)
            [1, 2, 4, 5, 8, 10, 13, 16, 17, 19]
        """
        v = [ZZ(1)]
        c = ZZ(2)
        while len(v) < n:
            if self.satisfies_heegner_hypothesis(D, c):
                v.append(c)
            c += 1
        return v

    def optimal_embeddings(self, D, c, R):
        """
        INPUT:

        - ``D`` -- negative fundamental discriminant

        - ``c`` -- integer coprime

        - ``R`` -- Eichler order

        EXAMPLES::

            sage: H = heegner_points(11).reduce_mod(3)
            sage: R = H.left_orders()[0]
            sage: H.optimal_embeddings(-7, 1, R)
            [Embedding sending sqrt(-7) to i - j - k,
             Embedding sending sqrt(-7) to -i + j + k]
            sage: H.optimal_embeddings(-7, 2, R)
            [Embedding sending 2*sqrt(-7) to 5*i - k,
             Embedding sending 2*sqrt(-7) to -5*i + k,
             Embedding sending 2*sqrt(-7) to 2*i - 2*j - 2*k,
             Embedding sending 2*sqrt(-7) to -2*i + 2*j + 2*k]
        """
        Q, G = R.ternary_quadratic_form(include_basis=True)
        n = -D*c*c
        reps = Q.representation_vector_list(n+1)[-1]

        # The representatives give elements in terms of the
        # subspace's basis such that the embedding is given by
        #     phi(c*sqrt(D)) = beta
        E = []
        for r in reps:
            beta = sum(G[i]*r[i] for i in range(len(G)))
            phi = HeegnerQuatAlgEmbedding(D, c, R, beta)
            E.append(phi)
        return E

    @cached_method
    def brandt_module(self):
        """
        Return the Brandt module of right ideal classes that we
        used to represent the set of supersingular points on
        the modular curve.

        EXAMPLES::

            sage: heegner_points(11).reduce_mod(3).brandt_module()
            Brandt module of dimension 2 of level 3*11 of weight 2 over Rational Field
        """
        from sage.modular.quatalg.all import BrandtModule
        return BrandtModule(self.__ell, self.__level)

    @cached_method
    def quaternion_algebra(self):
        """
        Return the rational quaternion algebra used to implement ``self``.

        EXAMPLES::

            sage: heegner_points(389).reduce_mod(7).quaternion_algebra()
            Quaternion Algebra (-1, -7) with base ring Rational Field
        """
        return self.brandt_module().quaternion_algebra()

    def right_ideals(self):
        """
        Return representative right ideals in the Brandt module.

        EXAMPLES::

            sage: heegner_points(11).reduce_mod(3).right_ideals()
            (Fractional ideal (4, 44*i, 2 + 8*i + 2*j, 34*i + 2*k),
             Fractional ideal (8, 88*i, 2 + 52*i + 2*j, 4 + 78*i + 2*k))
        """
        return self.brandt_module().right_ideals()

    @cached_method
    def left_orders(self):
        """
        Return the left orders associated to the representative right
        ideals in the Brandt module.

        EXAMPLES::

            sage: heegner_points(11).reduce_mod(3).left_orders()
            [Order of Quaternion Algebra (-1, -3) with base ring Rational Field
              with basis (1/2 + 1/2*j + 7*k, 1/2*i + 13/2*k, j + 3*k, 11*k),
             Order of Quaternion Algebra (-1, -3) with base ring Rational Field
              with basis (1/2 + 1/2*j + 7*k, 1/4*i + 1/2*j + 63/4*k, j + 14*k, 22*k)]
        """
        return [I.left_order() for I in self.right_ideals()]

    @cached_method
    def heegner_divisor(self, D, c=ZZ(1)):
        r"""
        Return Heegner divisor as an element of the Brandt module
        corresponding to the discriminant `D` and conductor `c`, which
        both must be coprime to `N\ell`.

        More precisely, we compute the sum of the reductions of the
        `\textrm{Gal}(K_1/K)`-conjugates of each choice of `y_1`,
        where the choice comes from choosing the ideal `\mathcal{N}`.
        Then we apply the Hecke operator `T_c` to this sum.

        INPUT:

        - ``D`` -- discriminant (negative integer)

        - ``c`` -- conductor (positive integer)

        OUTPUT: a Brandt module element

        EXAMPLES::

            sage: H = heegner_points(11).reduce_mod(7)
            sage: H.heegner_discriminants()
            [-8, -39, -43, -51, -79]
            sage: H.heegner_divisor(-8)
            (1, 0, 0, 1, 0, 0)
            sage: H.heegner_divisor(-39)
            (1, 2, 2, 1, 2, 0)
            sage: H.heegner_divisor(-43)
            (1, 0, 0, 1, 0, 0)
            sage: H.heegner_divisor(-51)
            (1, 0, 0, 1, 0, 2)
            sage: H.heegner_divisor(-79)
            (3, 2, 2, 3, 0, 0)

            sage: sum(H.heegner_divisor(-39).element())
            8
            sage: QuadraticField(-39,'a').class_number()
            4
        """
        if not self.satisfies_heegner_hypothesis(D, c):
            raise ValueError("D and c must be coprime to N and ell")

        B = self.brandt_module()

        if c > 1:
            # Just apply T_c to divisor for c=1
            z = self.heegner_divisor(D)
            return B.hecke_operator(c)(z)

        n = -D
        v = [0]*B.degree()
        for i, R in enumerate(self.left_orders()):
            Q = R.ternary_quadratic_form()
            a = Q.theta_series(n+1)[n]
            if a > 0:
                reps = Q.representation_vector_list(n+1)[-1]
                k = len([r for r in reps if gcd(r) == 1])
                assert k % 2 == 0
                v[i] += k // 2
        return B(v)

    @cached_method
    def modp_splitting_data(self, p):
        r"""
        Return mod `p` splitting data for the quaternion algebra at the
        unramified prime `p`.  This is a pair of `2\times 2` matrices
        `A`, `B` over the finite field `\GF{p}` such that if the
        quaternion algebra has generators `i, j, k`, then the
        homomorphism sending `i` to `A` and `j` to `B` maps any
        maximal order homomorphically onto the ring of `2\times 2` matrices.

        Because of how the homomorphism is defined, we must assume that the
        prime `p` is odd.

        INPUT:

        - ``p`` -- unramified odd prime

        OUTPUT: a 2-tuple of matrices over finite field

        EXAMPLES::

            sage: H = heegner_points(11).reduce_mod(7)
            sage: H.quaternion_algebra()
            Quaternion Algebra (-1, -7) with base ring Rational Field
            sage: I, J = H.modp_splitting_data(13)
            sage: I
            [ 0 12]
            [ 1  0]
            sage: J
            [7 3]
            [3 6]
            sage: I^2
            [12  0]
            [ 0 12]
            sage: J^2
            [6 0]
            [0 6]
            sage: I*J == -J*I
            True

        The following is a good test because of the asserts in the code::

            sage: v = [H.modp_splitting_data(p) for p in primes(13,200)]

        Some edge cases::

            sage: H.modp_splitting_data(11)
            (
            [ 0 10]  [6 1]
            [ 1  0], [1 5]
            )

        Proper error handling::

            sage: H.modp_splitting_data(7)
            Traceback (most recent call last):
            ...
            ValueError: p (=7) must be an unramified prime

            sage: H.modp_splitting_data(2)
            Traceback (most recent call last):
            ...
            ValueError: p must be odd
        """
        p = ZZ(p)
        if not p.is_prime():
            raise ValueError("p (=%s) must be prime" % p)
        if p == 2:
            raise ValueError("p must be odd")
        Q = self.quaternion_algebra()
        if Q.discriminant() % p == 0:
            raise ValueError("p (=%s) must be an unramified prime" % p)
        i, j, k = Q.gens()
        F = GF(p)
        i2 = F(i*i)
        j2 = F(j*j)
        M = MatrixSpace(F, 2)
        I = M([0,i2,1,0])
        i2inv = 1/i2
        a = None
        #for b in reversed(list(F)):
        for b in list(F):
            if not b:
                continue
            c = j2 + i2inv * b*b
            if c.is_square():
                a = -c.sqrt()
                break
        assert a is not None, "bug in that no splitting solution found"
        J = M([a,b,(j2-a*a)/b, -a])
        assert I*J == -J*I, "bug in that I,J do not skew commute"
        return I, J

    def modp_splitting_map(self, p):
        r"""
        Return (algebra) map from the (`p`-integral) quaternion algebra to
        the set of `2\times 2` matrices over `\GF{p}`.

        INPUT:

        - ``p`` -- prime number

        EXAMPLES::

            sage: H = heegner_points(11).reduce_mod(7)
            sage: f = H.modp_splitting_map(13)
            sage: B = H.quaternion_algebra(); B
            Quaternion Algebra (-1, -7) with base ring Rational Field
            sage: i, j, k = H.quaternion_algebra().gens()
            sage: a = 2 + i - j + 3*k; b = 7 + 2*i - 4*j + k
            sage: f(a*b)
            [12  3]
            [10  5]
            sage: f(a)*f(b)
            [12  3]
            [10  5]
        """
        I, J = self.modp_splitting_data(p)
        K = I * J
        F = I.base_ring()

        def phi(q):
            v = [F(a) for a in q.coefficient_tuple()]
            return v[0] + I * v[1] + J * v[2] + K * v[3]
        return phi

    def cyclic_subideal_p1(self, I, c):
        r"""
        Compute dictionary mapping 2-tuples that defined normalized
        elements of `P^1(\ZZ/c\ZZ)`

        INPUT:

        - ``I`` -- right ideal of Eichler order or in quaternion algebra

        - ``c`` -- square free integer (currently must be odd prime
          and coprime to level, discriminant, characteristic, etc.

        OUTPUT:

        - dictionary mapping 2-tuples (u,v) to ideals

        EXAMPLES::

            sage: H = heegner_points(11).reduce_mod(7)
            sage: I = H.brandt_module().right_ideals()[0]
            sage: sorted(H.cyclic_subideal_p1(I, 3).items())
            [((0, 1), Fractional ideal (12, 132*i, 10 + 76*i + 2*j, 4 + 86*i + 2*k)),
             ((1, 0), Fractional ideal (12, 132*i, 2 + 32*i + 2*j, 8 + 130*i + 2*k)),
             ((1, 1), Fractional ideal (12, 132*i, 10 + 32*i + 2*j, 8 + 86*i + 2*k)),
             ((1, 2), Fractional ideal (12, 132*i, 2 + 76*i + 2*j, 4 + 130*i + 2*k))]
            sage: len(H.cyclic_subideal_p1(I, 17))
            18
        """
        c = ZZ(c)
        if not c.is_prime():
            raise NotImplementedError("currently c must be prime")
        if c == 2:
            raise NotImplementedError("currently c must be odd")
        phi = self.modp_splitting_map(c)
        B = self.brandt_module()
        P1 = P1List(c)
        ans = {}
        # Actually they are submodules despite the name below.
        for J in B.cyclic_submodules(I, c):
            B = J.basis()
            V = phi(B[0]).kernel()
            for i in [1,2,3]:
                V = V.intersection(phi(B[i]).kernel())
            b = V.basis()
            assert len(b) == 1, "common kernel must have dimension 1"
            uv = P1.normalize(ZZ(b[0][0]) % c, ZZ(b[0][1]) % c)
            ans[uv] = J
        assert len(ans) == c+1
        return ans

    @cached_method
    def galois_group_over_hilbert_class_field(self, D, c):
        """
        Return the Galois group of the extension of ring class fields
        `K_c` over the Hilbert class field `K_{1}` of the quadratic
        imaginary field of discriminant `D`.

        INPUT:

        - ``D`` -- fundamental discriminant

        - ``c`` -- conductor (square-free integer)

        EXAMPLES::

            sage: N = 37; D = -7; ell = 17; c = 41; p = 3
            sage: H = heegner_points(N).reduce_mod(ell)
            sage: H.galois_group_over_hilbert_class_field(D, c)
            Galois group of Ring class field extension of QQ[sqrt(-7)] of conductor 41
             over Hilbert class field of QQ[sqrt(-7)]
        """
        Kc = heegner_points(self.level(), D, c).ring_class_field()
        K1 = heegner_points(self.level(), D, 1).ring_class_field()
        return Kc.galois_group(K1)

    @cached_method
    def galois_group_over_quadratic_field(self, D, c):
        """
        Return the Galois group of the extension of ring class fields
        `K_c` over the quadratic imaginary field `K` of discriminant `D`.

        INPUT:

        - ``D`` -- fundamental discriminant

        - ``c`` -- conductor (square-free integer)

        EXAMPLES::

            sage: N = 37; D = -7; ell = 17; c = 41; p = 3
            sage: H = heegner_points(N).reduce_mod(ell)
            sage: H.galois_group_over_quadratic_field(D, c)
            Galois group of Ring class field extension of QQ[sqrt(-7)] of conductor 41
             over Number Field in sqrt_minus_7 with defining polynomial x^2 + 7
              with sqrt_minus_7 = 2.645751311064591?*I
        """
        Kc = heegner_points(self.level(), D, c).ring_class_field()
        return Kc.galois_group(Kc.quadratic_field())

    @cached_method
    def quadratic_field(self, D):
        """
        Return our fixed choice of quadratic imaginary field of
        discriminant `D`.

        INPUT:

        - ``D`` -- fundamental discriminant

        OUTPUT: a quadratic number field

        EXAMPLES::

            sage: H = heegner_points(389).reduce_mod(5)
            sage: H.quadratic_field(-7)
            Number Field in sqrt_minus_7 with defining polynomial x^2 + 7
             with sqrt_minus_7 = 2.645751311064591?*I
        """
        Kc = heegner_points(self.level(), D, 1).ring_class_field()
        return Kc.quadratic_field()

    @cached_method
    def kolyvagin_cyclic_subideals(self, I, p, alpha_quaternion):
        r"""
        Return list of pairs `(J, n)` where `J` runs through the
        cyclic subideals of `I` of index `(\ZZ/p\ZZ)^2`, and `J \sim
        \alpha^n(J_0)` for some fixed choice of cyclic subideal `J_0`.

        INPUT:

        - ``I`` -- right ideal of the quaternion algebra

        - ``p`` -- prime number

        - ``alpha_quaternion`` -- image in the quaternion algebra
          of generator `\alpha` for `(\mathcal{O}_K / c\mathcal{O}_K)^* / (\ZZ/c\ZZ)^*`

        OUTPUT: list of 2-tuples

        EXAMPLES::

            sage: N = 37; D = -7; ell = 17; c = 5
            sage: H = heegner_points(N).reduce_mod(ell)
            sage: I = H.brandt_module().right_ideals()[49]
            sage: f = H.optimal_embeddings(D, 1, I.left_order())[1]
            sage: g = H.kolyvagin_generators(f.domain().number_field(), c)
            sage: alpha_quaternion = f(g[0]); alpha_quaternion
            1 - 77/192*i - 5/128*j - 137/384*k
            sage: H.kolyvagin_cyclic_subideals(I, 5, alpha_quaternion)
            [(Fractional ideal (2560, 1280 + 47360*i, 1146 + 37678*i + 4*j, 212 + 54664/3*i + 2*j + 2/3*k), 0),
             (Fractional ideal (2560, 1280 + 47360*i, 2426 + 9262*i + 4*j, 2004 + 83080/3*i + 2*j + 2/3*k), 1),
             (Fractional ideal (2560, 1280 + 47360*i, 1914 + 9262*i + 4*j, 1748 + 111496/3*i + 2*j + 2/3*k), 2),
             (Fractional ideal (2560, 1280 + 47360*i, 2170 + 18734*i + 4*j, 212 + 111496/3*i + 2*j + 2/3*k), 3),
             (Fractional ideal (2560, 1280 + 47360*i, 890 + 28206*i + 4*j, 1748 + 54664/3*i + 2*j + 2/3*k), 4),
             (Fractional ideal (2560, 1280 + 47360*i, 634 + 37678*i + 4*j, 2516 + 83080/3*i + 2*j + 2/3*k), 5)]
        """
        X = I.cyclic_right_subideals(p, alpha_quaternion)
        return [(J, i) for i, J in enumerate(X)]

    @cached_method
    def kolyvagin_generator(self, K, p):
        r"""
        Return element in `K` that maps to the multiplicative generator
        for the quotient group

           `(\mathcal{O}_K / p \mathcal{O}_K)^* / (\ZZ/p\ZZ)^*`

        of the form `\sqrt{D}+n` with `n\geq 1` minimal.

        INPUT:

        - ``K`` -- quadratic imaginary field

        - ``p`` -- inert prime

        EXAMPLES::

            sage: N = 37; D = -7; ell = 17; p = 5
            sage: H = heegner_points(N).reduce_mod(ell)
            sage: I = H.brandt_module().right_ideals()[49]
            sage: f = H.optimal_embeddings(D, 1, I.left_order())[0]
            sage: H.kolyvagin_generator(f.domain().number_field(), 5)
            a + 1

        This function requires that `p` be prime, but ``kolyvagin_generators`` works in general::

            sage: H.kolyvagin_generator(f.domain().number_field(), 5*17)
            Traceback (most recent call last):
            ...
            NotImplementedError: p must be prime
            sage: H.kolyvagin_generators(f.domain().number_field(), 5*17)
            [-34*a + 1, 35*a + 106]
        """
        p = ZZ(p)
        if not p.is_prime():
            raise NotImplementedError("p must be prime")
        if K.discriminant() % p == 0:
            raise ValueError("p must be unramified")
        if len(K.factor(p)) != 1:
            raise ValueError("p must be inert")

        F = K.residue_field(p)
        a = F.gen()
        assert a*a == K.discriminant(), "bug: we assumed generator of finite field must be square root of discriminant, but for some reason this is not true"
        for n in range(1,p):
            if (a + n).multiplicative_order() % (p*p-1) == 0:
                return K.gen() + n
        raise RuntimeError("there is a bug in kolyvagin_generator")

    @cached_method
    def kolyvagin_generators(self, K, c):
        r"""
        Return elements in `\mathcal{O}_K` that map to multiplicative generators
        for the factors of the quotient group

           `(\mathcal{O}_K / c \mathcal{O}_K)^* / (\ZZ/c\ZZ)^*`

        corresponding to the prime divisors of c.  Each generator is
        of the form `\sqrt{D}+n` with `n\geq 1` minimal.

        INPUT:

        - ``K`` -- quadratic imaginary field

        - ``c`` -- square free product of inert prime

        EXAMPLES::

            sage: N = 37; D = -7; ell = 17; p = 5
            sage: H = heegner_points(N).reduce_mod(ell)
            sage: I = H.brandt_module().right_ideals()[49]
            sage: f = H.optimal_embeddings(D, 1, I.left_order())[0]
            sage: H.kolyvagin_generators(f.domain().number_field(), 5*17)
            [-34*a + 1, 35*a + 106]
        """
        v = []
        F = ZZ(c).factor()
        from sage.rings.integer_ring import crt_basis
        B = crt_basis([x[0] for x in F])
        for i, (p, e) in enumerate(F):
            if e > 1:
                raise ValueError("c must be square free")
            alpha = self.kolyvagin_generator(K, p)
            # Now we use the Chinese Remainder Theorem to make an element
            # of O_K that equals alpha modulo p and equals 1 modulo
            # all other prime divisors of c.
            Z = [1]*len(B)
            Z[i] = alpha[0]
            a0 = sum([Z[j]*B[j] for j in range(len(B))])
            Z = [0]*len(B)
            Z[i] = alpha[1]
            a1 = sum([Z[j]*B[j] for j in range(len(B))])
            v.append(alpha.parent()([a0,a1]))
        return v

    @cached_method
    def kolyvagin_sigma_operator(self, D, c, r, bound=None):
        """
        Return the action of the Kolyvagin sigma operator on the `r`-th
        basis vector.

        INPUT:

        - ``D`` -- fundamental discriminant

        - ``c`` -- conductor (square-free integer, need not be prime)

        - ``r`` -- nonnegative integer

        - ``bound`` -- (default: ``None``), if given, controls
          precision of computation of theta series, which could
          impact performance, but does not impact correctness

        EXAMPLES:

        We first try to verify Kolyvagin's conjecture for a rank 2
        curve by working modulo 5, but we are unlucky with `c=17`::

            sage: N = 389; D = -7; ell = 5; c = 17; q = 3
            sage: H = heegner_points(N).reduce_mod(ell)
            sage: E = EllipticCurve('389a')
            sage: V = H.modp_dual_elliptic_curve_factor(E, q, 5)  # long time (4s on sage.math, 2012)
            sage: k118 = H.kolyvagin_sigma_operator(D, c, 118)
            sage: k104 = H.kolyvagin_sigma_operator(D, c, 104)
            sage: [b.dot_product(k104.element().change_ring(GF(3)))  # long time
            ....:  for b in V.basis()]
            [0, 0]
            sage: [b.dot_product(k118.element().change_ring(GF(3)))  # long time
            ....:  for b in V.basis()]
            [0, 0]

        Next we try again with `c=41` and this does work, in that we
        get something nonzero, when dotting with V::

            sage: c = 41
            sage: k118 = H.kolyvagin_sigma_operator(D, c, 118)
            sage: k104 = H.kolyvagin_sigma_operator(D, c, 104)
            sage: [b.dot_product(k118.element().change_ring(GF(3)))  # long time
            ....:  for b in V.basis()]
            [2, 0]
            sage: [b.dot_product(k104.element().change_ring(GF(3)))  # long time
            ....:  for b in V.basis()]
            [1, 0]

        By the way, the above is the first ever provable verification
        of Kolyvagin's conjecture for any curve of rank at least 2.

        Another example, but where the curve has rank 1::

            sage: N = 37; D = -7; ell = 17; c = 41; q = 3
            sage: H = heegner_points(N).reduce_mod(ell)
            sage: H.heegner_divisor(D,1).element().nonzero_positions()
            [49, 51]
            sage: k49 = H.kolyvagin_sigma_operator(D, c, 49); k49
            (79, 32, 31, 11, 53, 37, 1, 23, 15, 7, 0, 0, 0, 64, 32, 34, 53, 0, 27, 27, 0, 0, 0, 26, 0, 0, 18, 0, 22, 0, 53, 19, 27, 10, 0, 0, 0, 30, 35, 38, 0, 0, 0, 53, 0, 0, 4, 0, 0, 0, 0, 0)
            sage: k51 = H.kolyvagin_sigma_operator(D, c, 51); k51
            (20, 12, 57, 0, 0, 0, 0, 52, 23, 15, 0, 7, 0, 0, 19, 4, 0, 73, 11, 0, 104, 31, 0, 38, 31, 0, 0, 31, 5, 47, 0, 27, 35, 0, 57, 32, 24, 10, 0, 8, 0, 31, 41, 0, 0, 0, 16, 0, 0, 0, 0, 0)
            sage: V = H.modp_dual_elliptic_curve_factor(EllipticCurve('37a'), q, 5); V
            Vector space of degree 52 and dimension 2 over Ring of integers modulo 3
            Basis matrix:
            2 x 52 dense matrix over Ring of integers modulo 3
            sage: [b.dot_product(k49.element().change_ring(GF(q))) for b in V.basis()]
            [1, 1]
            sage: [b.dot_product(k51.element().change_ring(GF(q))) for b in V.basis()]
            [1, 1]

        An example with `c` a product of two primes::

            sage: N = 389; D = -7; ell = 5; q = 3
            sage: H = heegner_points(N).reduce_mod(ell)
            sage: V = H.modp_dual_elliptic_curve_factor(EllipticCurve('389a'), q, 5)
            sage: k = H.kolyvagin_sigma_operator(D, 17*41, 104)     # long time
            sage: k                                                 # long time
            (990, 656, 219, ..., 246, 534, 1254)
            sage: [b.dot_product(k.element().change_ring(GF(3))) for b in V.basis()]   # long time (but only because depends on something slow)
            [0, 0]
        """
        B = self.brandt_module()
        RI = B.right_ideals()

        f = self.optimal_embeddings(D, 1, RI[r].left_order())[0]
        alphas = self.kolyvagin_generators(f.domain().number_field(), c)
        alpha_quaternions = [f(x) for x in alphas]

        if bound is None:
            bound = B.dimension() // 2 + 5
        theta_dict = B._theta_dict(bound)

        c = ZZ(c)
        J_lists = []
        F = c.factor()
        I = RI[r]
        for i, (p, e) in enumerate(F):
            if e > 1:
                raise ValueError("c must be square free")
            X = I.cyclic_right_subideals(p, alpha_quaternions[i])
            J_lists.append(dict(enumerate(X)))

        ans = [0] * B.dimension()
        for v in product(*[range(1, p + 1) for p, _ in F]):
            J = J_lists[0][v[0]]
            for i in range(1, len(J_lists)):
                J = J.intersection(J_lists[i][v[i]])
            J_theta = tuple(J.theta_series_vector(bound))
            d = theta_dict[J_theta]
            j = None
            if len(d) == 1:
                j = d[0]
            else:
                for z in d:
                    if RI[z].is_equivalent(J, 0):
                        j = z
                        # we found the right j
                        break
            if j is None:
                raise RuntimeError("bug finding equivalent ideal")
            ans[j] += prod(v)
        return B(ans)

    @cached_method
    def modp_dual_elliptic_curve_factor(self, E, p, bound=10):
        """
        Return the factor of the Brandt module space modulo `p`
        corresponding to the elliptic curve `E`, cut out using
        Hecke operators up to ``bound``.

        INPUT:

        - ``E`` -- elliptic curve of conductor equal to the level of ``self``

        - ``p`` -- prime number

        - ``bound`` -- positive integer (default: 10)

        EXAMPLES::

            sage: N = 37; D = -7; ell = 17; c = 41; q = 3
            sage: H = heegner_points(N).reduce_mod(ell)
            sage: V = H.modp_dual_elliptic_curve_factor(EllipticCurve('37a'), q, 5); V
            Vector space of degree 52 and dimension 2 over Ring of integers modulo 3
             Basis matrix: 2 x 52 dense matrix over Ring of integers modulo 3
        """
        if E.conductor() != self.level():
            raise ValueError("conductor of E must equal level of self")
        p = ZZ(p)
        if not p.is_prime():
            raise ValueError("p (=%s) must be prime" % p)
        bad = self.__level * self.__ell

        V = None
        q = ZZ(2)
        B = self.brandt_module()
        F = GF(p)
        while q <= bound and (V is None or V.dimension() > 2):
            verbose("q = %s" % q)
            if bad % q != 0:
                T = B._compute_hecke_matrix_directly(q).change_ring(F).transpose()
                if V is None:
                    V = (T - E.ap(q)).kernel()
                else:
                    t = T.restrict(V)
                    W = (t - E.ap(q)).kernel()
                    V = (W.basis_matrix() * V.basis_matrix()).row_space()
            q = q.next_prime()
        return V

    @cached_method
    def rational_kolyvagin_divisor(self, D, c):
        r"""
        Return the Kolyvagin divisor as an element of the Brandt module
        corresponding to the discriminant `D` and conductor `c`, which
        both must be coprime to `N\ell`.

        INPUT:

        - ``D`` -- discriminant (negative integer)

        - ``c`` -- conductor (positive integer)

        OUTPUT: Brandt module element (or tuple of them)

        EXAMPLES::

            sage: N = 389; D = -7; ell = 5; c = 17; q = 3
            sage: H = heegner_points(N).reduce_mod(ell)
            sage: k = H.rational_kolyvagin_divisor(D, c); k  # long time (5s on sage.math, 2013)
            (2, 0, 0, 0, 0, 0, 16, 0, 0, 0, 0, 4, 0, 0, 9, 11, 0, 6, 0, 0, 7, 0, 0, 0, 0, 14, 12, 13, 15, 17, 0, 0, 0, 0, 8, 0, 0, 0, 0, 10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
            sage: V = H.modp_dual_elliptic_curve_factor(EllipticCurve('389a'), q, 2)
            sage: [b.dot_product(k.element().change_ring(GF(q))) for b in V.basis()]  # long time
            [0, 0]
            sage: k = H.rational_kolyvagin_divisor(D, 59)
            sage: [b.dot_product(k.element().change_ring(GF(q))) for b in V.basis()]
            [2, 0]
        """
        if not self.satisfies_heegner_hypothesis(D, c):
            raise ValueError("D and c must be coprime to N and ell")

        hd = self.heegner_divisor(D)
        v = hd.element()
        if class_number(D) != 1:
            raise NotImplementedError("class number greater than 1 not implemented")
        i = min(v.nonzero_positions())
        return self.kolyvagin_sigma_operator(D, c, i)

    @cached_method
    def kolyvagin_point_on_curve(self, D, c, E, p, bound=10):
        r"""
        Compute image of the Kolyvagin divisor `P_c` in
        `E(\GF{\ell^2}) / p E(\GF{\ell^2})`.

        Note that this image is by definition only well defined up to
        scalars.  However, doing multiple computations will always
        yield the same result, and working modulo different `\ell` is
        compatible (since we always choose the same generator for
        `\textrm{Gal}(K_c/K_1)`).

        INPUT:

        - ``D`` -- fundamental negative discriminant

        - ``c`` -- conductor

        - ``E`` -- elliptic curve of conductor the level of self

        - ``p`` -- odd prime number such that we consider image in
          `E(\GF{\ell^2}) / p E(\GF{\ell^2})`

        - ``bound`` -- integer (default: 10)

        EXAMPLES::

            sage: N = 37; D = -7; ell = 17; c = 41; p = 3
            sage: H = heegner_points(N).reduce_mod(ell)
            sage: H.kolyvagin_point_on_curve(D, c, EllipticCurve('37a'), p)
            [1, 1]
        """
        k = self.rational_kolyvagin_divisor(D, c)
        V = self.modp_dual_elliptic_curve_factor(E, p, bound)
        return [b.dot_product(k.element().change_ring(GF(p))) for b in V.basis()]


def kolyvagin_reduction_data(E, q, first_only=True):
    r"""
    Given an elliptic curve of positive rank and a prime `q`, this
    function returns data about how to use Kolyvagin's `q`-torsion
    Heegner point Euler system to do computations with this curve.
    See the precise description of the output below.

    INPUT:

    - ``E`` -- elliptic curve over `\QQ` of rank 1 or 2

    - ``q`` -- an odd prime that does not divide the order of the
      rational torsion subgroup of `E`

    - ``first_only`` -- boolean (default: ``True``); whether two only return
      the first prime that one can work modulo to get data about the Euler system

    OUTPUT in the rank 1 case or when the default flag ``first_only=True``:

    - `\ell` -- first good odd prime satisfying the Kolyvagin condition that
      `q` divides \gcd(a_{\ell},\ell+1)` and the reduction map is surjective to
      `E(\GF{\ell}) / q E(\GF{\ell})`

    - ``D`` -- discriminant of the first quadratic imaginary field
       `K` that satisfies the Heegner hypothesis for `E` such that
       both `\ell` is inert in `K`, and the twist `E^D` has analytic
       rank `\leq 1`

    - ``h_D`` -- the class number of `K`

    - the dimension of the Brandt module `B(\ell,N)`, where `N` is
      the conductor of `E`

    OUTPUT in the rank 2 case:

    - `\ell_1` -- first prime (as above in the rank 1 case) where
      reduction map is surjective

    - `\ell_2` -- second prime (as above) where reduction map is
      surjective

    - ``D`` -- discriminant of the first quadratic imaginary field
       `K` that satisfies the Heegner hypothesis for `E` such that
       both `\ell_1` and `\ell_2` are simultaneously inert in `K`,
       and the twist `E^D` has analytic rank `\leq 1`

    - ``h_D`` -- the class number of `K`

    - the dimension of the Brandt module `B(\ell_1,N)`, where `N` is
      the conductor of `E`

    - the dimension of the Brandt module `B(\ell_2,N)`

    EXAMPLES:

    Import this function::

        sage: from sage.schemes.elliptic_curves.heegner import kolyvagin_reduction_data

    A rank 1 example::

        sage: kolyvagin_reduction_data(EllipticCurve('37a1'), 3)
        (17, -7, 1, 52)

    A rank 3 example::

        sage: kolyvagin_reduction_data(EllipticCurve('5077a1'), 3)
        (11, -47, 5, 4234)
        sage: H = heegner_points(5077, -47)
        sage: [c for c in H.kolyvagin_conductors(2, 10, EllipticCurve('5077a1'), 3)
        ....:    if c % 11]
        [667, 943, 1189, 2461]
        sage: factor(667)
        23 * 29

    A rank 4 example (the first Kolyvagin class that we could try to
    compute would be `P_{23\cdot 29\cdot 41}`, and would require
    working in a space of dimension 293060 (so prohibitive at
    present)::

        sage: E = elliptic_curves.rank(4)[0]
        sage: kolyvagin_reduction_data(E,3)              # long time
        (11, -71, 7, 293060)
        sage: H = heegner_points(293060, -71)
        sage: H.kolyvagin_conductors(1,4,E,3)
        [11, 17, 23, 41]

    The first rank 2 example::

        sage: kolyvagin_reduction_data(EllipticCurve('389a'), 3)
        (5, -7, 1, 130)
        sage: kolyvagin_reduction_data(EllipticCurve('389a'), 3, first_only=False)
        (5, 17, -7, 1, 130, 520)

    A large `q = 7`::

        sage: kolyvagin_reduction_data(EllipticCurve('1143c1'), 7, first_only=False)
        (13, 83, -59, 3, 1536, 10496)

    Additive reduction::

        sage: kolyvagin_reduction_data(EllipticCurve('2350g1'), 5, first_only=False)
        (19, 239, -311, 19, 6480, 85680)
    """
    from .ell_generic import EllipticCurve_generic
    if not isinstance(E, EllipticCurve_generic):
        raise TypeError("E must be an elliptic curve")

    q = ZZ(q)
    if not q.is_prime():
        raise ValueError("q must be a prime")

    if q.gcd(E.torsion_order()) != 1:
        raise NotImplementedError("q must be coprime to torsion")

    N = E.conductor()
    r = E.rank()

    if r == 0:
        raise ValueError("E must have positive rank")

    if E.rank() == 1:
        first_only = True

    from sage.modular.quatalg.all import BrandtModule

    def twist_is_minimal(D):
        # return True if the quadratic twist E^D has analytic rank <= 1
        return E.quadratic_twist(D).analytic_rank() <= 1

    def red(P, ell):
        # reduce the point P on the elliptic curve modulo ell
        w = list(P)
        d = lcm([a.denominator() for a in w])
        return E.change_ring(GF(ell))([d*a for a in w])

    def best_heegner_D(ell_1, ell_2):
        # return the first Heegner D satisfy all hypothesis such that
        # both ell_1 and ell_2 are inert
        D = ZZ(-5)
        while True:
            if D.is_fundamental_discriminant() and \
               D % ell_1 and D % ell_2 and \
               E.satisfies_heegner_hypothesis(D) and \
               is_inert(D, ell_1) and is_inert(D, ell_2) and \
               twist_is_minimal(D):
                return D
            D -= 1

    if first_only:
        # find first prime ell with various conditions
        # such that reduction is surjective to E(F_ell)/q.
        ell = ZZ(3)
        while True:
            while N % ell == 0 or gcd(ell+1,E.ap(ell)) % q != 0:
                ell = ell.next_prime()
            # determine if mod ell reduction is surjective, using
            # partly that it is a lemma that E(F_ell)/q is cyclic.
            m = ZZ(E.Np(ell) / q)
            for P in E.gens():
                if red(P,ell) * m != 0:
                    # bingo, is surjective
                    D = best_heegner_D(ell,ell)
                    return (ell, D, class_number(D), BrandtModule(ell,N).dimension())
            # end for
            ell = ell.next_prime()

    if E.rank() != 2:
        raise ValueError("if first_only is not True, then the curve E must have rank 1 or 2")

    P, Q = E.gens()

    def kernel_of_reduction(ell):
        # return list of reps for the kernel as a subgroup of the map
        # E(Q) / q E(Q)  ---->  E(F_ell) / q E(F_ell)
        m = ZZ(E.Np(ell) / q)
        A = (a * P + b * Q for a in range(q) for b in range(q))
        return [z for z in A if red(z, ell) * m == 0]

    # compute first good odd prime
    ell_1 = ZZ(3)
    while True:
        while N % ell_1 == 0 or gcd(ell_1+1,E.ap(ell_1)) % q != 0:
            ell_1 = ell_1.next_prime()
        # compute kernel of reduction modulo ell_1
        G1 = set(kernel_of_reduction(ell_1))
        if len(G1) == q:
            break
        ell_1 = ell_1.next_prime()

    # compute next good odd prime with distinct kernel of order q
    ell_2 = ell_1.next_prime()
    while True:
        while N % ell_2 == 0 or gcd(ell_2+1,E.ap(ell_2)) % q != 0:
            ell_2 = ell_2.next_prime()
        G2 = set(kernel_of_reduction(ell_2))
        if G1 != G2 and len(G2) == q:
            break
        ell_2 = ell_2.next_prime()

    # Find smallest D where both ell_1 and ell_2 are inert
    D = best_heegner_D(ell_1, ell_2)
    return (ell_1, ell_2, D, class_number(D),
            BrandtModule(ell_1,N).dimension(),
            BrandtModule(ell_2,N).dimension())


class HeegnerQuatAlgEmbedding(SageObject):
    r"""
    The homomorphism `\mathcal{O} \to R`, where `\mathcal{O}` is the
    order of conductor `c` in the quadratic field of discriminant `D`,
    and `R` is an Eichler order in a quaternion algebra.

    EXAMPLES::

        sage: H = heegner_points(11).reduce_mod(3); R = H.left_orders()[0]
        sage: f = H.optimal_embeddings(-7, 2, R)[1]; f
        Embedding sending 2*sqrt(-7) to -5*i + k
        sage: type(f)
        <class 'sage.schemes.elliptic_curves.heegner.HeegnerQuatAlgEmbedding'>
        sage: loads(dumps(f)) == f
        True
    """
    def __init__(self, D, c, R, beta):
        r"""
        INPUT:

        - ``D`` -- negative fundamental discriminant

        - ``c`` -- positive integer coprime to `D`

        - ``R`` -- Eichler order in a rational quaternion algebra

        - ``beta`` -- element of `R` such that the homomorphism
          sends `c\sqrt{D}` to `\beta`

        EXAMPLES::

            sage: H = heegner_points(11).reduce_mod(3); R = H.left_orders()[0]
            sage: i,j,k = H.quaternion_algebra().gens()
            sage: import sage.schemes.elliptic_curves.heegner as heegner
            sage: heegner.HeegnerQuatAlgEmbedding(-7, 2, R, -5*i+k)
            Embedding sending 2*sqrt(-7) to -5*i + k
        """
        self.__D = D
        self.__c = c
        self.__R = R
        self.__beta = beta

    def __eq__(self, other):
        """
        EXAMPLES::

            sage: H = heegner_points(11).reduce_mod(3); R = H.left_orders()[0]
            sage: f = H.optimal_embeddings(-7, 2, R)[0]
            sage: f == H.optimal_embeddings(-7, 2, R)[0]
            True
            sage: f == H.optimal_embeddings(-7, 2, R)[1]
            False
            sage: f == 0
            False
        """
        return isinstance(other, HeegnerQuatAlgEmbedding) and \
               self.__D == other.__D and \
               self.__c == other.__c and \
               self.__R == other.__R and \
               self.__beta == other.__beta

    def __ne__(self, other):
        """
        EXAMPLES::

            sage: H = heegner_points(11).reduce_mod(3); R = H.left_orders()[0]
            sage: f = H.optimal_embeddings(-7, 2, R)[0]
            sage: f != H.optimal_embeddings(-7, 2, R)[0]
            False
            sage: f != H.optimal_embeddings(-7, 2, R)[1]
            True
            sage: f != 0
            True
        """
        return not (self == other)

    def __call__(self, x):
        """
        Return image of `x` under this embedding.

        INPUT:

        - ``x`` -- element of the quadratic order

        EXAMPLES::

            sage: H = heegner_points(11).reduce_mod(3); R = H.left_orders()[0]
            sage: f = H.optimal_embeddings(-7, 1, R)[1]; f
            Embedding sending sqrt(-7) to -i + j + k
            sage: a = f.domain_gen(); a^2
            -7
            sage: f(2 + 3*a)
            2 - 3*i + 3*j + 3*k
            sage: 2 + 3*f(a)
            2 - 3*i + 3*j + 3*k
            sage: f(a)^2
            -7
        """
        v = self.domain().number_field()(x).vector()
        w = v * self.matrix()
        z = self.codomain().quaternion_algebra()(w)
        # There is no notion of an "element of an order" implemented
        # for quaternion algebras right now.  All elements are
        # elements of the ambient rational quaternion algebra.
        return z

    @cached_method
    def matrix(self):
        r"""
        Return matrix over `\QQ` of this morphism, with respect to the
        basis 1, `c\sqrt{D}` of the domain and the basis `1,i,j,k` of
        the ambient rational quaternion algebra (which contains the
        domain).

        EXAMPLES::

            sage: H = heegner_points(11).reduce_mod(3); R = H.left_orders()[0]
            sage: f = H.optimal_embeddings(-7, 1, R)[1]; f
            Embedding sending sqrt(-7) to -i + j + k
            sage: f.matrix()
            [ 1  0  0  0]
            [ 0 -1  1  1]
            sage: f.conjugate().matrix()
            [ 1  0  0  0]
            [ 0  1 -1 -1]
        """
        return matrix(QQ,2,4,[[1,0,0,0], self.__beta.coefficient_tuple()])

    @cached_method
    def domain(self):
        """
        Return the domain of this embedding.

        EXAMPLES::

            sage: H = heegner_points(11).reduce_mod(3); R = H.left_orders()[0]
            sage: H.optimal_embeddings(-7, 2, R)[0].domain()
            Order of conductor 4 generated by 2*a
             in Number Field in a with defining polynomial x^2 + 7
             with a = 2.645751311064591?*I
        """
        R, a = quadratic_order(self.__D, self.__c)

        # The following assumption is used, e.g., in the __call__
        # method.  I know that it is satisfied by the current
        # implementation.  But somebody might someday annoying change
        # the implementation, and we want to catch that if it were to
        # ever happen.

        assert R.basis() == [1, a], "an assumption about construction of orders is violated"
        self.__domain_gen = a
        return R

    def domain_gen(self):
        r"""
        Return the specific generator `c \sqrt{D}` for the domain
        order.

        EXAMPLES::

            sage: H = heegner_points(11).reduce_mod(3); R = H.left_orders()[0]
            sage: f = H.optimal_embeddings(-7, 2, R)[0]
            sage: f.domain_gen()
            2*a
            sage: f.domain_gen()^2
            -28
        """
        self.domain()
        return self.__domain_gen

    def domain_conductor(self):
        """
        Return the conductor of the domain.

        EXAMPLES::

            sage: H = heegner_points(11).reduce_mod(3); R = H.left_orders()[0]
            sage: H.optimal_embeddings(-7, 2, R)[0].domain_conductor()
            2
        """
        return self.__c

    def beta(self):
        r"""
        Return the element `\beta` in the quaternion algebra order
        that `c\sqrt{D}` maps to.

        EXAMPLES::

            sage: H = heegner_points(11).reduce_mod(3); R = H.left_orders()[0]
            sage: H.optimal_embeddings(-7, 2, R)[1].beta()
            -5*i + k
        """
        return self.__beta

    def codomain(self):
        """
        Return the codomain of this embedding.

        EXAMPLES::

            sage: H = heegner_points(11).reduce_mod(3); R = H.left_orders()[0]
            sage: H.optimal_embeddings(-7, 2, R)[0].codomain()
            Order of Quaternion Algebra (-1, -3) with base ring Rational Field
             with basis (1/2 + 1/2*j + 7*k, 1/2*i + 13/2*k, j + 3*k, 11*k)
        """
        return self.__R

    @cached_method
    def _repr_(self):
        """
        Return string representation of this embedding.

        EXAMPLES::

            sage: H = heegner_points(11).reduce_mod(3); R = H.left_orders()[0]
            sage: f = H.optimal_embeddings(-7, 2, R)[1]; f._repr_()
            'Embedding sending 2*sqrt(-7) to -5*i + k'
        """
        a = '%ssqrt(%s)' % ('%s*' % self.__c if self.__c > 1 else '', self.__D)
        return "Embedding sending %s to %s" % (a, self.__beta)

    def conjugate(self):
        """
        Return the conjugate of this embedding, which is also an
        embedding.

        EXAMPLES::

            sage: H = heegner_points(11).reduce_mod(3); R = H.left_orders()[0]
            sage: f = H.optimal_embeddings(-7, 2, R)[1]
            sage: f.conjugate()
            Embedding sending 2*sqrt(-7) to 5*i - k
            sage: f
            Embedding sending 2*sqrt(-7) to -5*i + k
        """
        return HeegnerQuatAlgEmbedding(self.__D, self.__c,
                                       self.__R, self.__beta.conjugate())


#############################################################################
# Utility Functions
#############################################################################


def quadratic_order(D, c, names='a'):
    r"""
    Return order of conductor `c` in quadratic field with fundamental
    discriminant `D`.

    INPUT:

    - ``D`` -- fundamental discriminant

    - ``c`` -- conductor

    - ``names`` -- string (default: ``'a'``)

    OUTPUT:

    - order `R` of conductor `c` in an imaginary quadratic field

    - the element `c\sqrt{D}` as an element of `R`

    The generator for the field is named 'a' by default.

    EXAMPLES::

        sage: sage.schemes.elliptic_curves.heegner.quadratic_order(-7,3)
        (Order of conductor 6 generated by 3*a in Number Field in a
         with defining polynomial x^2 + 7 with a = 2.645751311064591?*I,
         3*a)
        sage: sage.schemes.elliptic_curves.heegner.quadratic_order(-7,3,'alpha')
        (Order of conductor 6 generated by 3*alpha in Number Field in alpha
         with defining polynomial x^2 + 7 with alpha = 2.645751311064591?*I,
         3*alpha)
    """
    K = QuadraticField(D, names)
    sqrtD = K.gen(0)
    t = sqrtD * c
    R = K.order([t])
    return R, R(t)


def class_number(D):
    """
    Return the class number of the quadratic field with fundamental
    discriminant `D`.

    INPUT:

    - ``D`` -- integer

    EXAMPLES::

        sage: sage.schemes.elliptic_curves.heegner.class_number(-20)
        2
        sage: sage.schemes.elliptic_curves.heegner.class_number(-23)
        3
        sage: sage.schemes.elliptic_curves.heegner.class_number(-163)
        1

    A :exc:`ValueError` is raised when `D` is not a fundamental
    discriminant::

        sage: sage.schemes.elliptic_curves.heegner.class_number(-5)
        Traceback (most recent call last):
        ...
        ValueError: D (=-5) must be a fundamental discriminant
    """
    if not D.is_fundamental_discriminant():
        raise ValueError("D (=%s) must be a fundamental discriminant" % D)
    return D.class_number()


def is_inert(D, p):
    r"""
    Return ``True`` if p is an inert prime in the field `\QQ(\sqrt{D})`.

    INPUT:

    - ``D`` -- fundamental discriminant

    - ``p`` -- prime integer

    EXAMPLES::

        sage: sage.schemes.elliptic_curves.heegner.is_inert(-7,3)
        True
        sage: sage.schemes.elliptic_curves.heegner.is_inert(-7,7)
        False
        sage: sage.schemes.elliptic_curves.heegner.is_inert(-7,11)
        False
    """
    K = QuadraticField(D,'a')
    F = K.factor(p)
    return len(F) == 1 and F[0][1] == 1


def is_split(D, p):
    r"""
    Return ``True`` if p is a split prime in the field `\QQ(\sqrt{D})`.

    INPUT:

    - ``D`` -- fundamental discriminant

    - ``p`` -- prime integer

    EXAMPLES::

        sage: sage.schemes.elliptic_curves.heegner.is_split(-7,3)
        False
        sage: sage.schemes.elliptic_curves.heegner.is_split(-7,7)
        False
        sage: sage.schemes.elliptic_curves.heegner.is_split(-7,11)
        True
    """
    K = QuadraticField(D,'a')
    F = K.factor(p)
    return len(F) == 2


def is_ramified(D, p):
    r"""
    Return ``True`` if p is a ramified prime in the field `\QQ(\sqrt{D})`.

    INPUT:

    - ``D`` -- fundamental discriminant

    - ``p`` -- prime integer

    EXAMPLES::

        sage: sage.schemes.elliptic_curves.heegner.is_ramified(-7,2)
        False
        sage: sage.schemes.elliptic_curves.heegner.is_ramified(-7,7)
        True
        sage: sage.schemes.elliptic_curves.heegner.is_ramified(-1,2)
        True
    """
    return QuadraticField(D,'a').discriminant() % p == 0


def nearby_rational_poly(f, **kwds):
    r"""
    Return a polynomial whose coefficients are rational numbers close
    to the coefficients of `f`.

    INPUT:

    - ``f`` -- polynomial with real floating point entries

    - ``**kwds`` -- passed on to ``nearby_rational`` method

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.heegner import nearby_rational_poly
        sage: R.<x> = RR[]
        sage: nearby_rational_poly(2.1*x^2 + 3.5*x - 1.2, max_error=10e-16)
        21/10*X^2 + 7/2*X - 6/5
        sage: nearby_rational_poly(2.1*x^2 + 3.5*x - 1.2, max_error=10e-17)
        4728779608739021/2251799813685248*X^2 + 7/2*X - 5404319552844595/4503599627370496
        sage: RR(4728779608739021/2251799813685248  - 21/10)
        8.88178419700125e-17
    """
    R = QQ['X']
    return R([a.nearby_rational(**kwds) for a in f])


def simplest_rational_poly(f, prec):
    """
    Return a polynomial whose coefficients are as simple as possible
    rationals that are also close to the coefficients of f.

    INPUT:

    - ``f`` -- polynomial with real floating point entries

    - ``prec`` -- positive integer

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.heegner import simplest_rational_poly
        sage: R.<x> = RR[]
        sage: simplest_rational_poly(2.1*x^2 + 3.5*x - 1.2, 53)
        21/10*X^2 + 7/2*X - 6/5
    """
    R = QQ['X']
    Z = RealField(prec)
    return R([Z(a).simplest_rational() for a in f])


def satisfies_weak_heegner_hypothesis(N, D):
    r"""
    Check that `D` satisfies the weak Heegner hypothesis relative to `N`.
    This is all that is needed to define Heegner points.

    The condition is that `D<0` is a fundamental discriminant and that
    each unramified prime dividing `N` splits in `K=\QQ(\sqrt{D})` and
    each ramified prime exactly divides `N`.  We also do not require
    that `D<-4`.

    INPUT:

    - ``N`` -- positive integer

    - ``D`` -- negative integer

    EXAMPLES::

        sage: s = sage.schemes.elliptic_curves.heegner.satisfies_weak_heegner_hypothesis
        sage: s(37,-7)
        True
        sage: s(37,-37)
        False
        sage: s(37,-37*4)
        True
        sage: s(100,-4)
        False
        sage: [D for D in [-1,-2,..,-40] if s(37,D)]
        [-3, -4, -7, -11, -40]
        sage: [D for D in [-1,-2,..,-100] if s(37,D)]
        [-3, -4, -7, -11, -40, -47, -67, -71, -83, -84, -95]
        sage: EllipticCurve('37a').heegner_discriminants_list(10)
        [-7, -11, -40, -47, -67, -71, -83, -84, -95, -104]
    """
    if not D.is_fundamental_discriminant():
        return False
    if D >= 0:
        return False
    for p, e in N.factor():
        if D % p == 0:
            if e > 1:
                return False
        elif D.kronecker(p) != 1:
            return False
    return True


def make_monic(f):
    r"""
    Return a monic integral polynomial `g` and an integer `d` such
    that if `\alpha` is a root of `g`, then `\alpha/d` is a root of `f`.
    In other words, `c f(x) = g(d x)` for some scalar `c`.

    INPUT:

    - ``f`` -- polynomial over the rational numbers

    OUTPUT: a monic integral polynomial and an integer

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.heegner import make_monic
        sage: R.<x> = QQ[]
        sage: make_monic(3*x^3 + 14*x^2 - 7*x + 5)
        (x^3 + 14*x^2 - 21*x + 45, 3)

    In this example we verify that ``make_monic`` does what we claim it does::

        sage: K.<a> = NumberField(x^3 + 17*x - 3)
        sage: f = (a/7+2/3).minpoly(); f
        x^3 - 2*x^2 + 247/147*x - 4967/9261
        sage: g, d = make_monic(f); (g, d)
        (x^3 - 42*x^2 + 741*x - 4967, 21)
        sage: K.<b> = NumberField(g)
        sage: (b/d).minpoly()
        x^3 - 2*x^2 + 247/147*x - 4967/9261

    TESTS::

        sage: f = x^5 + x^3/4 + 5
        sage: make_monic(f)
        (x^5 + x^3 + 160, 2)

    Scalar factors do not matter, the result is always monic::

        sage: make_monic(f * 1000000)
        (x^5 + x^3 + 160, 2)
        sage: make_monic(f / 1000000)
        (x^5 + x^3 + 160, 2)
    """
    R = f.parent()
    n = f.degree()
    lc = f[n]
    d = ZZ.one()
    for i in range(n):
        expo = n - i
        # We require that (d^expo * f[i] / lc) is an integer
        den = (d**expo * f[i] / lc).denominator()
        for p, e in factor_trial_division(den, 1000000):
            # Round up e/expo
            d *= p ** ((e + expo - 1) // expo)
    g = R([d**(n-i) * f[i] / lc for i in range(n+1)])
    return g, d


#####################################################################
# Elliptic curve methods
# Everywhere self below is an elliptic curve over QQ.
#####################################################################

def ell_heegner_point(self, D, c=ZZ(1), f=None, check=True):
    r"""
    Return the Heegner point on this curve associated to the
    quadratic imaginary field `K=\QQ(\sqrt{D})`.

    If the optional parameter `c` is given, returns the higher Heegner
    point associated to the order of conductor `c`.

    INPUT:

    - ``D`` -- a Heegner discriminant

    - ``c`` -- (default: 1) conductor, must be coprime to `DN`

    - ``f`` -- binary quadratic form or 3-tuple `(A,B,C)` of coefficients
      of `AX^2 + BXY + CY^2`

    - ``check`` -- boolean (default: ``True``)

    OUTPUT: the Heegner point `y_c`

    EXAMPLES::

        sage: E = EllipticCurve('37a')
        sage: E.heegner_discriminants_list(10)
        [-7, -11, -40, -47, -67, -71, -83, -84, -95, -104]
        sage: P = E.heegner_point(-7); P                          # indirect doctest
        Heegner point of discriminant -7 on elliptic curve of conductor 37
        sage: z = P.point_exact(); z == E(0, 0, 1)  or -z == E(0, 0, 1)
        True
        sage: P.curve()
        Elliptic Curve defined by y^2 + y = x^3 - x over Rational Field
        sage: P = E.heegner_point(-40).point_exact(); P
        (a : -a + 1 : 1)
        sage: P = E.heegner_point(-47).point_exact(); P
        (a : a^4 + a - 1 : 1)
        sage: P[0].parent()
        Number Field in a with defining polynomial x^5 - x^4 + x^3 + x^2 - 2*x + 1

    Working out the details manually::

        sage: P = E.heegner_point(-47).numerical_approx(prec=200)
        sage: f = algdep(P[0], 5); f
        x^5 - x^4 + x^3 + x^2 - 2*x + 1
        sage: f.discriminant().factor()
        47^2

    The Heegner hypothesis is checked::

        sage: E = EllipticCurve('389a'); P = E.heegner_point(-5,7);
        Traceback (most recent call last):
        ...
        ValueError: N (=389) and D (=-5) must satisfy the Heegner hypothesis

    We can specify the quadratic form::

        sage: P = EllipticCurve('389a').heegner_point(-7, 5, (778,925,275)); P
        Heegner point of discriminant -7 and conductor 5
         on elliptic curve of conductor 389
        sage: P.quadratic_form()
        778*x^2 + 925*x*y + 275*y^2
    """
    y = HeegnerPointOnX0N(self.conductor(), D, c, f, check=check)
    return y.map_to_curve(self)


def kolyvagin_point(self, D, c=ZZ(1), check=True):
    r"""
    Return the Kolyvagin point on this curve associated to the
    quadratic imaginary field `K=\QQ(\sqrt{D})` and conductor `c`.

    INPUT:

    - ``D`` -- a Heegner discriminant

    - ``c`` -- (default: 1) conductor, must be coprime to `DN`

    - ``check`` -- boolean (default: ``True``)

    OUTPUT: the Kolyvagin point `P` of conductor `c`

    EXAMPLES::

        sage: E = EllipticCurve('37a1')
        sage: P = E.kolyvagin_point(-67); P
        Kolyvagin point of discriminant -67 on elliptic curve of conductor 37
        sage: P.numerical_approx()  # abs tol 1e-14
        (6.00000000000000 : -15.0000000000000 : 1.00000000000000)
        sage: P.index()
        6
        sage: g = E((0,-1,1)) # a generator
        sage: E.regulator() == E.regulator_of_points([g])
        True
        sage: 6*g
        (6 : -15 : 1)
    """
    return self.heegner_point(D,c,check=check).kolyvagin_point()


def ell_heegner_discriminants(self, bound):
    """
    Return the list of ``self``'s Heegner discriminants between -1 and
    -bound.

    INPUT:

    - ``bound`` -- integer; upper bound for -discriminant

    OUTPUT: the list of Heegner discriminants between -1 and -bound for
    the given elliptic curve

    EXAMPLES::

        sage: E=EllipticCurve('11a')
        sage: E.heegner_discriminants(30)                     # indirect doctest
        [-7, -8, -19, -24]
    """
    return [ZZ(-D) for D in range(1, bound)
            if self.satisfies_heegner_hypothesis(-D)]


def ell_heegner_discriminants_list(self, n):
    """
    Return the list of ``self``'s first `n` Heegner discriminants smaller
    than -5.

    INPUT:

    - ``n`` -- integer; the number of discriminants to compute

    OUTPUT: the list of the first `n` Heegner discriminants smaller than
    `-5` for the given elliptic curve

    EXAMPLES::

        sage: E=EllipticCurve('11a')
        sage: E.heegner_discriminants_list(4)                     # indirect doctest
        [-7, -8, -19, -24]
    """
    v = []
    D = ZZ(-5)
    while len(v) < n:
        while not self.satisfies_heegner_hypothesis(D):
            D -= 1
        v.append(D)
        D -= 1
    return v


def heegner_point_height(self, D, prec=2, check_rank=True):
    r"""
    Use the Gross-Zagier formula to compute the Neron-Tate canonical
    height over `K` of the Heegner point corresponding to `D`, as an
    interval (it is computed to some precision using `L`-functions).

    If the curve has rank at least 2, then the returned height is the
    exact Sage integer 0.

    INPUT:

    - ``D`` -- integer; fundamental discriminant (=/= -3, -4)

    - ``prec`` -- integer (default: 2); use `prec \cdot \sqrt(N) + 20`
      terms of `L`-series in computations, where `N` is the
      conductor

    - ``check_rank`` -- whether to check if the rank is at least 2 by
      computing the Mordell-Weil rank directly

    OUTPUT: interval that contains the height of the Heegner point

    EXAMPLES::

        sage: E = EllipticCurve('11a')
        sage: E.heegner_point_height(-7)
        0.22227?

    Some higher rank examples::

        sage: E = EllipticCurve('389a')
        sage: E.heegner_point_height(-7)
        0
        sage: E = EllipticCurve('5077a')
        sage: E.heegner_point_height(-7)
        0
        sage: E.heegner_point_height(-7, check_rank=False)
        0.0000?
    """

    if not self.satisfies_heegner_hypothesis(D):
        raise ArithmeticError("Discriminant (=%s) must be a fundamental discriminant that satisfies the Heegner hypothesis." % D)

    if check_rank and self.rank() >= 2:
        return ZZ(0)

    if D == -3 or D == -4:
        raise ArithmeticError("Discriminant (=%s) must not be -3 or -4." % D)
    eps = self.root_number()
    L1_vanishes = self.lseries().L1_vanishes()

    IR = RealIntervalField(20)    # TODO: why 20 bits here?

    if eps == 1 and L1_vanishes:
        return IR(0) # rank even hence >= 2, so Heegner point is torsion.

    RR = RealField()
    from math import sqrt

    alpha = RR(sqrt(abs(D)))/(2*self.period_lattice().complex_area())
    F = self.quadratic_twist(D)
    E = self
    k_E = prec*sqrt(E.conductor()) + 20
    k_F = prec*sqrt(F.conductor()) + 20

    MIN_ERR = RR('1e-6')  # we assume that regulator and
                         # discriminant, etc., computed to this accuracy (which is easily the case).
                         # this should be made more intelligent / rigorous relative
                         # to the rest of the system.

    if eps == 1:   # E has even rank
        LF1, err_F = F.lseries().deriv_at1(k_F)
        LE1, err_E = E.lseries().at1(k_E)
        err_F = max(err_F, MIN_ERR)
        err_E = max(err_E, MIN_ERR)
        return IR(alpha-MIN_ERR,alpha+MIN_ERR) * IR(LE1-err_E,LE1+err_E) * IR(LF1-err_F,LF1+err_F)

    else:          # E has odd rank
        LE1, err_E = E.lseries().deriv_at1(k_E)
        LF1, err_F = F.lseries().at1(k_F)
        err_F = max(err_F, MIN_ERR)
        err_E = max(err_E, MIN_ERR)
        return IR(alpha-MIN_ERR,alpha+MIN_ERR) * IR(LE1-err_E,LE1+err_E) * IR(LF1-err_F,LF1+err_F)


def heegner_index(self, D, min_p=2, prec=5, descent_second_limit=12,
                  verbose_mwrank=False, check_rank=True):
    r"""
    Return an interval that contains the index of the Heegner
    point `y_K` in the group of `K`-rational points modulo torsion
    on this elliptic curve, computed using the Gross-Zagier
    formula and/or a point search, or possibly half the index
    if the rank is greater than one.

    If the curve has rank > 1, then the returned index is infinity.

    .. NOTE::

        If ``min_p`` is bigger than 2 then the index can be off by
        any prime less than ``min_p``. This function returns the
        index divided by `2` exactly when the rank of `E(K)` is
        greater than 1 and `E(\QQ)_{/tor} \oplus E^D(\QQ)_{/tor}`
        has index `2` in `E(K)_{/tor}`, where the second factor
        undergoes a twist.

    INPUT:

    - ``D`` -- integer; Heegner discriminant

    - ``min_p`` -- integer (default: 2); only rule out primes
      = min_p dividing the index

    - ``verbose_mwrank`` -- boolean (default: ``False``); print lots of
      mwrank search status information when computing regulator

    - ``prec`` -- integer (default: 5);, use prec\*sqrt(N) +
      20 terms of `L`-series in computations, where N is the conductor

    - ``descent_second_limit`` -- (default: 12)- used in 2-descent
      when computing regulator of the twist

    - ``check_rank`` -- whether to check if the rank is at least 2 by
      computing the Mordell-Weil rank directly

    OUTPUT: an interval that contains the index, or half the index

    EXAMPLES::

        sage: E = EllipticCurve('11a')
        sage: E.heegner_discriminants(50)
        [-7, -8, -19, -24, -35, -39, -40, -43]
        sage: E.heegner_index(-7)
        1.00000?

    ::

        sage: E = EllipticCurve('37b')
        sage: E.heegner_discriminants(100)
        [-3, -4, -7, -11, -40, -47, -67, -71, -83, -84, -95]
        sage: E.heegner_index(-95)          # long time (1 second)
        2.00000?

    This tests doing direct computation of the Mordell-Weil group.

    ::

        sage: EllipticCurve('675b').heegner_index(-11)
        3.0000?

    Currently discriminants -3 and -4 are not supported::

        sage: E.heegner_index(-3)
        Traceback (most recent call last):
        ...
        ArithmeticError: Discriminant (=-3) must not be -3 or -4.

    The curve 681b returns the true index, which is `3`::

        sage: E = EllipticCurve('681b')
        sage: I = E.heegner_index(-8); I
        3.0000?

    In fact, whenever the returned index has a denominator of
    `2`, the true index is got by multiplying the returned
    index by `2`. Unfortunately, this is not an if and only if
    condition, i.e., sometimes the index must be multiplied by
    `2` even though the denominator is not `2`.

    This example demonstrates the ``descent_second_limit`` option,
    which can be used to fine tune the 2-descent used to compute
    the regulator of the twist::

        sage: E = EllipticCurve([1,-1,0,-1228,-16267])
        sage: E.heegner_index(-8)
        Traceback (most recent call last):
        ...
        RuntimeError: ...

    However when we search higher, we find the points we need::

        sage: E.heegner_index(-8, descent_second_limit=16, check_rank=False)  # long time
        2.00000?

    Two higher rank examples (of ranks 2 and 3)::

        sage: E = EllipticCurve('389a')
        sage: E.heegner_index(-7)
        +Infinity
        sage: E = EllipticCurve('5077a')
        sage: E.heegner_index(-7)
        +Infinity
        sage: E.heegner_index(-7, check_rank=False)
        0.001?
        sage: E.heegner_index(-7, check_rank=False).lower() == 0
        True
    """
    if not self.satisfies_heegner_hypothesis(D):
        raise ArithmeticError("Discriminant (=%s) must be a fundamental discriminant that satisfies the Heegner hypothesis." % D)

    if check_rank and self.rank() >= 2:
        return infinity

    # First compute upper bound on height of Heegner point.
    tm = verbose("computing heegner point height...")
    h0 = self.heegner_point_height(D, prec=prec, check_rank=check_rank)
    if h0 == 0:
        return infinity

    # We divide by 2 to get the height **over Q** of the
    # Heegner point on the twist.

    ht = h0/2
    verbose('Height of heegner point = %s' % ht, tm)

    if self.root_number() == 1:
        F = self.quadratic_twist(D)
    else:
        F = self
    # Now rank(F) > 0
    h = ht.upper()
    verbose("Heegner height bound = %s" % h)
    B = F.CPS_height_bound()
    verbose("CPS bound = %s" % B)
    c = h/(min_p**2) + B
    verbose("Search would have to be up to height = %s" % c)

    from .ell_rational_field import _MAX_HEIGHT

    IR = RealIntervalField(20)  # todo: 20?

    a = 1
    if c > _MAX_HEIGHT or F is self:
        verbose("Doing direct computation of MW group.")
        reg = F.regulator(descent_second_limit=descent_second_limit, verbose=verbose_mwrank)
        if F.rank(use_database=True) == 1:
            z = F.gens()[0]
            FK = F.base_extend(QuadraticField(D,'a'))
            z = FK(z)
            if z.is_divisible_by(2):
                a = 2
            else:
                FK_even_tor_pts = [T for T in FK.torsion_subgroup().gens() if T.order() % 2 == 0]
                if len(FK_even_tor_pts) == 2:
                    FK_even_tor_pts.append(sum(FK_even_tor_pts))
                for T in FK_even_tor_pts:
                    if (z + T).is_divisible_by(2):
                        a = 2
                        break
        return a*self._adjust_heegner_index(ht/IR(reg))

    # Do naive search to eliminate possibility that Heegner point
    # is divisible by p<min_p, without finding Heegner point.
    verbose("doing point search")
    P = F.point_search(c)
    verbose("done with point search")
    P = [x for x in P if x.order() == infinity]
    a = 1
    if len(P) == 0:
        return IR(1)
    elif len(P) == 1:
        z = P[0]
        FK = F.base_extend(QuadraticField(D,'a'))
        z = FK(z)
        if z.is_divisible_by(2):
            a = 2
        else:
            FK_even_tor_pts = [T for T in FK.torsion_subgroup().gens() if T.order() % 2 == 0]
            if len(FK_even_tor_pts) == 2:
                FK_even_tor_pts.append(sum(FK_even_tor_pts))
            for T in FK_even_tor_pts:
                if (z + T).is_divisible_by(2):
                    a = 2
                    break

    verbose("saturating")
    S, I, reg = F.saturation(P)
    verbose("done saturating")
    return a*self._adjust_heegner_index(ht/IR(reg))


def _adjust_heegner_index(self, a):
    r"""
    Take the square root of the interval that contains the Heegner
    index.

    EXAMPLES::

        sage: E = EllipticCurve('11a1')
        sage: a = RIF(sqrt(2)) - RIF(1.4142135623730951)
        sage: E._adjust_heegner_index(a)
        1.?e-8
    """
    if a.lower() < 0:
        IR = RealIntervalField(20)  # todo: 20?
        a = IR((0, a.upper()))
    return a.sqrt()


def heegner_index_bound(self, D=0, prec=5, max_height=None):
    r"""
    Assume ``self`` has rank 0.

    Return a list `v` of primes such that if an odd prime `p` divides
    the index of the Heegner point in the group of rational points
    modulo torsion, then `p` is in `v`.

    If 0 is in the interval of the height of the Heegner point
    computed to the given prec, then this function returns `v =
    0`. This does not mean that the Heegner point is torsion, just
    that it is very likely torsion.

    If we obtain no information from a search up to ``max_height``,
    e.g., if the Siksek et al. bound is bigger than ``max_height``,
    then we return `v = -1`.

    INPUT:

    - ``D`` -- integer (default: 0); Heegner discriminant; if
      0, use the first discriminant -4 that satisfies the Heegner
      hypothesis

    - ``verbose`` -- boolean (default: ``True``)

    - ``prec`` -- integer (default: 5); use `prec \cdot \sqrt(N) + 20`
      terms of `L`-series in computations, where `N` is the conductor

    - ``max_height (float)`` -- should be = 21; bound on
      logarithmic naive height used in point searches. Make smaller to
      make this function faster, at the expense of possibly obtaining a
      worse answer. A good range is between 13 and 21.

    OUTPUT:

    - ``v`` -- list or int (bad primes or 0 or -1)

    - ``D`` -- the discriminant that was used (this is
      useful if `D` was automatically selected)

    - ``exact`` -- either False, or the exact Heegner index
      (up to factors of 2)

    EXAMPLES::

        sage: E = EllipticCurve('11a1')
        sage: E.heegner_index_bound()
        ([2], -7, 2)
    """
    from .ell_rational_field import _MAX_HEIGHT
    if max_height is None:
        max_height = _MAX_HEIGHT
    else:
        max_height = min(float(max_height), _MAX_HEIGHT)
    if self.root_number() != 1:
        raise RuntimeError("The rank must be 0.")

    if D == 0:
        D = -5
        while not self.satisfies_heegner_hypothesis(D):
            D -= 1

    # First compute upper bound on Height of Heegner point.
    ht = self.heegner_point_height(D, prec=prec)
    if 0 in ht:
        return 0, D, False
    F = self.quadratic_twist(D)
    h = ht.upper()
    verbose("Heegner height bound = %s" % h)
    B = F.CPS_height_bound()
    verbose("CPS bound = %s" % B)
    if self.two_torsion_rank() == 0:
        H = h
    else:
        H = 4*h
    p = 3
    from sage.arith.misc import next_prime
    while True:
        c = H/(2*p**2) + B
        if c < max_height:
            break
        if p > 100:
            break
        p = next_prime(p)
    verbose("Using p = %s" % p)

    if c > max_height:
        verbose("No information by searching only up to max_height (=%s)." % c)
        return -1, D, False

    verbose("Searching up to height = %s" % c)
    eps = 10e-5

    def _bound(P):
        """
        We will use this function below in two places. It bounds the index
        using a nontrivial point.
        """
        assert len(P) == 1

        S, I, reg = F.saturation(P)

        IR = RealIntervalField(20)  # todo: 20?
        h = IR(reg-eps,reg+eps)
        ind2 = ht/(h/2)
        verbose("index squared = %s" % ind2)
        ind = ind2.sqrt()
        verbose("index = %s" % ind)
        # Compute upper bound on square root of index.
        if ind.absolute_diameter() < 1:
            t, i = ind.is_int()
            if t:   # unique integer in interval, so we've found exact index squared.
                return prime_divisors(i), D, i
        raise RuntimeError("Unable to compute bound for e=%s, D=%s (try increasing precision)" % (self, D))

    # First try a quick search, in case we get lucky and find
    # a generator.
    P = F.point_search(13, rank_bound=1)
    P = [x for x in P if x.order() == infinity]
    if len(P) > 0:
        return _bound(P)

    # Do search to eliminate possibility that Heegner point is
    # divisible by primes up to p, without finding Heegner point.
    P = F.point_search(c, rank_bound=1)
    P = [x for x in P if x.order() == infinity]
    if len(P) == 0:
        # We've eliminated the possibility of a divisor up to p.
        return prime_range(3, p), D, False
    else:
        return _bound(P)


#################################################################################
def _heegner_index_in_EK(self, D):
    r"""
    Return the index of the sum of `E(\QQ)/tor + E^D(\QQ)/tor` in `E(K)/tor`.

    INPUT:

    - ``D`` -- negative integer; the Heegner discriminant

    OUTPUT: a power of 2 -- the given index

    EXAMPLES:

    We compute the index for a rank 2 curve and found that it is 2::

        sage: E = EllipticCurve('389a')
        sage: E._heegner_index_in_EK(-7)
        2

    We explicitly verify in the above example that indeed that
    index is divisible by 2 by writing down a generator of
    `E(\QQ)/tor + E^D(\QQ)/tor` that is divisible by 2 in `E(K)`::

        sage: F = E.quadratic_twist(-7)
        sage: K = QuadraticField(-7,'a')
        sage: G = E.change_ring(K)
        sage: phi = F.change_ring(K).isomorphism_to(G)
        sage: P = G(E(-1,1)) + G((0,-1)) + G(phi(F(14,25))); P
        (-867/3872*a - 3615/3872 : -18003/170368*a - 374575/170368 : 1)
        sage: P.division_points(2)
        [(1/8*a + 5/8 : -5/16*a - 9/16 : 1)]
    """
    # check conditions, then use cache if possible.
    if not self.satisfies_heegner_hypothesis(D):
        raise ValueError("D (=%s) must satisfy the Heegner hypothesis" % D)
    try:
        return self.__heegner_index_in_EK[D]
    except AttributeError:
        self.__heegner_index_in_EK = {}
    except KeyError:
        pass

    #####################################################################
    # THE ALGORITHM:
    #
    # For an element P of an abelian group A, let [P] denote the
    # equivalence class of P in the quotient A/A_tor of A by
    # its torsion subgroup.   Then for P in E(Q) + E^D(QQ), we
    # have that [P] is divisible by 2 in E(K)/tor if and only
    # there is R in E(K) such that 2*[R] = [P], and this is
    # only if there is R in E(K) and t in E(K)_tor such that
    #          2*R = P + t.
    #
    # Using complex conjugation, one sees that the quotient
    # group E(K)/tor / ( E(Q)/tor + E^D(Q)/tor ) is killed by 2.
    # So to compute the order of this group we run through
    # representatives P for A/(2A) where A = E(Q)/tor + E^D(Q)/tor,
    # and for each we see whether there is a torsion point t in E(K)
    # such that P + t is divisible by 2.   Also, we have
    #    2 | P+t  <==> 2 | P+n*t for any odd integer n,
    # so we may assume t is of 2-power order.
    #####################################################################

    E = self  # nice shortcut
    F = E.quadratic_twist(D).minimal_model()
    K = QuadraticField(D, 'a')

    # Define a map phi that we'll use to put the points of E^D(QQ)
    # into E(K):
    G = E.change_ring(K)
    G2 = F.change_ring(K)
    phi = G2.isomorphism_to(G)

    # Basis for E(Q)/tor oplus E^D(QQ)/tor in E(K):
    basis = [G(z) for z in E.gens()] + [G(phi(z)) for z in F.gens()]
    # Make a list of the 2-power order torsion points in E(K), including 0.
    T = [G(z) for z in G.torsion_subgroup().list() if z.order() == 1 or
            (z.order() % 2 == 0 and len(z.order().factor()) == 1)]

    r = len(basis)   # rank
    V = QQ**r
    B = []

    # Iterate through reps for A/(2*A) creating vectors in (1/2)*ZZ^r
    for v in GF(2)**r:
        if not v:
            continue
        P = sum([basis[i] for i in range(r) if v[i]])
        w = V(v) / 2
        B.extend(w for t in T if (P + t).is_divisible_by(2))

    A = ZZ**r
    # Take span of our vectors in (1/2)*ZZ^r, along with ZZ^r.  This is E(K)/tor.
    W = V.span(B, ZZ) + A

    # Compute the index in E(K)/tor of A = E(Q)/tor + E^D(Q)/tor, cache, and return.
    index = A.index_in(W)
    self.__heegner_index_in_EK[D] = index
    return index


def heegner_sha_an(self, D, prec=53):
    r"""
    Return the conjectural (analytic) order of Sha for E over the field `K=\QQ(\sqrt{D})`.

    INPUT:

    - ``D`` -- negative integer; the Heegner discriminant

    - ``prec`` -- integer (default: 53); bits of precision to
      compute analytic order of Sha

    OUTPUT:

    (floating point number) an approximation to the conjectural order of Sha.

    .. NOTE::

        Often you'll want to do ``proof.elliptic_curve(False)`` when
        using this function, since often the twisted elliptic
        curves that come up have enormous conductor, and Sha is
        nontrivial, which makes provably finding the Mordell-Weil
        group using 2-descent difficult.

    EXAMPLES:

    An example where E has conductor 11::

        sage: E = EllipticCurve('11a')
        sage: E.heegner_sha_an(-7)                                  # long time
        1.00000000000000

    The cache works::

        sage: E.heegner_sha_an(-7) is E.heegner_sha_an(-7)          # long time
        True

    Lower precision::

        sage: E.heegner_sha_an(-7,10)                               # long time
        1.0

    Checking that the cache works for any precision::

        sage: E.heegner_sha_an(-7,10) is E.heegner_sha_an(-7,10)    # long time
        True

    Next we consider a rank 1 curve with nontrivial Sha over the
    quadratic imaginary field `K`; however, there is no Sha for `E`
    over `\QQ` or for the quadratic twist of `E`::

        sage: E = EllipticCurve('37a')
        sage: E.heegner_sha_an(-40)                                 # long time
        4.00000000000000
        sage: E.quadratic_twist(-40).sha().an()                     # long time
        1
        sage: E.sha().an()                                          # long time
        1

    A rank 2 curve::

        sage: E = EllipticCurve('389a')                             # long time
        sage: E.heegner_sha_an(-7)                                  # long time
        1.00000000000000

    If we remove the hypothesis that `E(K)` has rank 1 in Conjecture
    2.3 in [GZ1986]_ page 311, then that conjecture is
    false, as the following example shows::

        sage: # long time
        sage: E = EllipticCurve('65a')
        sage: E.heegner_sha_an(-56)
        1.00000000000000
        sage: E.torsion_order()
        2
        sage: E.tamagawa_product()
        1
        sage: E.quadratic_twist(-56).rank()
        2
    """
    # check conditions, then return from cache if possible.
    if not self.satisfies_heegner_hypothesis(D):
        raise ValueError("D (=%s) must satisfy the Heegner hypothesis" % D)
    try:
        return self.__heegner_sha_an[(D, prec)]
    except AttributeError:
        self.__heegner_sha_an = {}
    except KeyError:
        pass

    # Use the BSD conjecture over the quadratic imaginary K --
    # see page 311 of [GZ1986]_ for the formula.
    E = self  # notational convenience
    F = E.quadratic_twist(D).minimal_model()
    K = QuadraticField(D, 'a')

    # Compute each of the quantities in BSD
    #  - The torsion subgroup over K.
    T = E.change_ring(K).torsion_order()

    #  - The product of the Tamagawa numbers, which because D is
    #    coprime to N is just the square of the product of the
    #    Tamagawa numbers over QQ for E.  (we square below in the
    #    BSD formula)
    cqprod = E.tamagawa_product()

    #  - The leading term of the L-series, as a product of two
    #  other L-series.
    rE = E.rank()
    rF = F.rank()
    L_E = E.lseries().dokchitser(prec).derivative(1, rE)
    L_F = F.lseries().dokchitser(prec).derivative(1, rF)
    #    NOTE: The binomial coefficient in the following formula
    #    for the leading term in terms of the other two leading
    #    terms comes from the product rule for the derivative.
    #    You can think this through or just type something like
    #      f = function('f',x); g = function('g',x); diff(f*g,6)
    #    into Sage to be convinced.
    L = binomial(rE + rF, rE) * (L_E * L_F / factorial(rE+rF) )

    #  - ||omega||^2 -- the period.  It is twice the volume of the
    #    period lattice.  See the following paper for a derivation:
    #    "Verification of the Birch and Swinnerton-Dyer Conjecture
    #     for Specific Elliptic Curves", G. Grigorov, A. Jorza, S. Patrikis,
    #     C. Patrascu, W. Stein
    omega = 2 * abs(E.period_lattice().basis_matrix().det())

    #  - The regulator.
    #    First we compute the regulator of the subgroup E(QQ) + E^D(QQ)
    #    of E(K).   The factor of 2 in the regulator
    #    accounts for the fact that the height over K is twice the
    #    height over QQ, i.e., for P in E(QQ) we have h_K(P,P) =
    #    2*h_Q(P,P).  See, e.g., equation (6.4) on page 230 of
    #    [GZ1986]_.
    Reg_prod = 2**(rE + rF) * E.regulator(precision=prec) * F.regulator(precision=prec)
    #    Next we call off to the _heegner_index_in_EK function, which
    #    saturates the group E(QQ) + E^D(QQ) in E(K), given us the index,
    #    which must be a power of 2, since E(QQ) is the +1 eigenspace for
    #    complex conjugation, and E^D(QQ) is the -1 eigenspace.
    ind = self._heegner_index_in_EK(D)
    #    Finally, we know the regulator of E(K).
    Reg = Reg_prod / ind**2

    #  - Square root of the absolute value of the discriminant.  This is
    #    easy; we just make sure the D passed in is an integer, so we
    #    can call sqrt with the chosen precision.
    sqrtD = ZZ(abs(D)).sqrt(prec=prec)

    #  - Done: Finally, we plug everything into the BSD formula to get the
    #    analytic order of Sha.
    sha_an = (L * T**2 * sqrtD) / (omega * Reg * cqprod**2)

    #  - We cache and return the answer.
    self.__heegner_sha_an[(D, prec)] = sha_an
    return sha_an


def _heegner_forms_list(self, D, beta=None, expected_count=None):
    r"""
    Return a list of quadratic forms corresponding to Heegner points
    with discriminant `D` and a choice of `\beta` a square root of
    `D` mod `4N`. Specifically, given a quadratic form
    `f = Ax^2 + Bxy + Cy^2` we let `\tau_f` be a root of `Ax^2 + Bx + C`
    and the discriminant `\Delta(\tau_f) = \Delta(f) = D` must be
    invariant under multiplication by `N`, the conductor of ``self``.

        `\Delta(N\tau_f) = \Delta(\tau_f) = \Delta(f) = D`

    EXAMPLES::

        sage: E = EllipticCurve('37a')
        sage: E._heegner_forms_list(-7)
        [37*x^2 + 17*x*y + 2*y^2]
        sage: E._heegner_forms_list(-195)
        [37*x^2 + 29*x*y + 7*y^2, 259*x^2 + 29*x*y + y^2, 111*x^2 + 177*x*y + 71*y^2, 2627*x^2 + 177*x*y + 3*y^2]
        sage: E._heegner_forms_list(-195)[-1].discriminant()
        -195
        sage: len(E._heegner_forms_list(-195))
        4
        sage: QQ[sqrt(-195)].class_number()
        4

        sage: E = EllipticCurve('389a')
        sage: E._heegner_forms_list(-7)
        [389*x^2 + 185*x*y + 22*y^2]
        sage: E._heegner_forms_list(-59)
        [389*x^2 + 313*x*y + 63*y^2, 1167*x^2 + 313*x*y + 21*y^2, 3501*x^2 + 313*x*y + 7*y^2]
    """
    if expected_count is None:
        expected_count = D.class_number()
    N = self.conductor()
    if beta is None:
        beta = Integers(4*N)(D).sqrt(extend=False)
    else:
        assert beta**2 == Integers(4*N)(D)
    from sage.quadratic_forms.binary_qf import BinaryQF
    b = ZZ(beta) % (2*N)
    all = []
    seen = []
    # Note: This may give a sub-optimal list of forms.
    while True:
        R = (b**2-D)//(4*N)
        for d in R.divisors():
            f = BinaryQF([d*N, b, R//d])
            fr = f.reduced_form()
            if fr not in seen:
                seen.append(fr)
                all.append(f)
                if len(all) == expected_count:
                    return all
        b += 2*N


def _heegner_best_tau(self, D, prec=None):
    r"""
    Given a discriminant `D`, find the Heegner point `\tau` in the
    upper half plane with largest imaginary part (which is optimal
    for evaluating the modular parametrization). If the optional
    parameter ``prec`` is given, return `\tau` to ``prec`` bits of
    precision, otherwise return it exactly as a symbolic object.

    EXAMPLES::

        sage: E = EllipticCurve('37a')
        sage: E._heegner_best_tau(-7)
        1/74*sqrt(-7) - 17/74
        sage: EllipticCurve('389a')._heegner_best_tau(-11)
        1/778*sqrt(-11) - 355/778
        sage: EllipticCurve('389a')._heegner_best_tau(-11, prec=100)
        -0.45629820051413881748071979434 + 0.0042630138693514136878083968338*I
    """
    # We know that N|A, so A = N is optimal.
    N = self.conductor()
    b = ZZ(Integers(4*N)(D).sqrt(extend=False) % (2*N))
    # TODO: make sure a different choice of b is not better?
    return (-b + ZZ(D).sqrt(prec=prec)) / (2*N)


def satisfies_heegner_hypothesis(self, D):
    """
    Return ``True`` precisely when `D` is a fundamental discriminant that
    satisfies the Heegner hypothesis for this elliptic curve.

    EXAMPLES::

        sage: E = EllipticCurve('11a1')
        sage: E.satisfies_heegner_hypothesis(-7)
        True
        sage: E.satisfies_heegner_hypothesis(-11)
        False
    """
    D = ZZ(D)
    if D >= 0:
        return False
    if not D.is_fundamental_discriminant():
        return False
    if D.gcd(self.conductor()) != 1:
        return False
    for p, _ in self.conductor().factor():
        if D.kronecker(p) != 1:
            return False
    return True


#####################################################################
# End of elliptic curve methods.
#####################################################################
