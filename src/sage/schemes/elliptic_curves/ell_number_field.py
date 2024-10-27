# sage.doctest: needs sage.rings.number_field
r"""
Elliptic curves over number fields

An elliptic curve `E` over a number field `K` can be given
by a Weierstrass equation whose coefficients lie in `K` or by
using ``base_extend`` on an elliptic curve defined over a subfield.

One major difference to elliptic curves over `\QQ` is that there
might not exist a global minimal equation over `K`, when `K` does
not have class number one.
Another difference is the lack of understanding of modularity for
general elliptic curves over general number fields.

Currently Sage can obtain local information about `E/K_v` for finite places
`v`, it has an interface to Denis Simon's script for 2-descent, it can compute
the torsion subgroup of the Mordell-Weil group `E(K)`, and it can work with
isogenies defined over `K`.

EXAMPLES::

    sage: K.<i> = NumberField(x^2 + 1)
    sage: E = EllipticCurve([0, 4+i])
    sage: E.discriminant()
    -3456*i - 6480
    sage: P= E([i,2])
    sage: P+P
    (-2*i + 9/16 : -9/4*i - 101/64 : 1)

::

    sage: E.has_good_reduction(2 + i)
    True
    sage: E.local_data(4+i)
    Local data at Fractional ideal (i + 4):
      Reduction type: bad additive
      Local minimal model: Elliptic Curve defined by y^2 = x^3 + (i+4)
                           over Number Field in i with defining polynomial x^2 + 1
      Minimal discriminant valuation: 2
      Conductor exponent: 2
      Kodaira Symbol: II
      Tamagawa Number: 1
    sage: E.tamagawa_product_bsd()
    1

::

    sage: E.simon_two_descent()
    (1, 1, [(i : 2 : 1)])

::

    sage: E.torsion_order()
    1

::

    sage: E.isogenies_prime_degree(3)
    [Isogeny of degree 3
      from Elliptic Curve defined by y^2 = x^3 + (i+4)
           over Number Field in i with defining polynomial x^2 + 1
        to Elliptic Curve defined by y^2 = x^3 + (-27*i-108)
           over Number Field in i with defining polynomial x^2 + 1]

AUTHORS:

- Robert Bradshaw 2007

- John Cremona

- Chris Wuthrich

REFERENCE:

- [Sil] Silverman, Joseph H. The arithmetic of elliptic curves. Second edition. Graduate Texts in
  Mathematics, 106. Springer, 2009.

- [Sil2] Silverman, Joseph H. Advanced topics in the arithmetic of elliptic curves. Graduate Texts in
  Mathematics, 151. Springer, 1994.
"""

# ****************************************************************************
#       Copyright (C) 2007 Robert Bradshaw <robertwb@math.washington.edu>
#                          William Stein <wstein@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

import sage.rings.abc
from .ell_field import EllipticCurve_field
from .ell_generic import EllipticCurve_generic
from .ell_point import EllipticCurvePoint_number_field
from .constructor import EllipticCurve
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.rings.real_mpfr import RealField
from sage.rings.integer import Integer
from sage.misc.cachefunc import cached_method
from sage.misc.misc_c import prod


class EllipticCurve_number_field(EllipticCurve_field):
    r"""
    Elliptic curve over a number field.

    EXAMPLES::

        sage: K.<i> = NumberField(x^2 + 1)
        sage: EllipticCurve([i, i - 1, i + 1, 24*i + 15, 14*i + 35])
        Elliptic Curve defined by
         y^2 + i*x*y + (i+1)*y = x^3 + (i-1)*x^2 + (24*i+15)*x + (14*i+35)
         over Number Field in i with defining polynomial x^2 + 1
    """
    def __init__(self, K, ainvs):
        r"""
        EXAMPLES:

        A curve from the database of curves over `\QQ`, but over a larger field::

            sage: K.<i> = NumberField(x^2 + 1)
            sage: EllipticCurve(K,'389a1')
            Elliptic Curve defined by y^2 + y = x^3 + x^2 + (-2)*x
             over Number Field in i with defining polynomial x^2 + 1

        Making the field of definition explicitly larger::

            sage: EllipticCurve(K,[0,-1,1,0,0])
            Elliptic Curve defined by y^2 + y = x^3 + (-1)*x^2
             over Number Field in i with defining polynomial x^2 + 1
        """
        self._known_points = []
        EllipticCurve_field.__init__(self, K, ainvs)

    _point = EllipticCurvePoint_number_field

    def base_extend(self, R):
        """
        Return the base extension of ``self`` to `R`.

        EXAMPLES::

            sage: E = EllipticCurve('11a3')
            sage: K = QuadraticField(-5, 'a')
            sage: E.base_extend(K)
            Elliptic Curve defined by y^2 + y = x^3 + (-1)*x^2 over Number Field in a
             with defining polynomial x^2 + 5 with a = 2.236067977499790?*I

        Check that non-torsion points are remembered when extending
        the base field (see :issue:`16034`)::

            sage: E = EllipticCurve([1, 0, 1, -1751, -31352])
            sage: K.<d> = QuadraticField(5)
            sage: E.gens()
            [(52 : 111 : 1)]
            sage: EK = E.base_extend(K)
            sage: EK.gens()
            [(52 : 111 : 1)]
        """
        E = super().base_extend(R)
        if isinstance(E, EllipticCurve_number_field):
            E._known_points = [E([R(_) for _ in P.xy()]) for P in self._known_points if not P.is_zero()]
        return E

    def simon_two_descent(self, verbose=0, lim1=2, lim3=4, limtriv=2,
                          maxprob=20, limbigprime=30, known_points=None):
        r"""
        Return lower and upper bounds on the rank of the Mordell-Weil
        group `E(K)` and a list of points.

        This method is used internally by the :meth:`~rank`,
        :meth:`~rank_bounds` and :meth:`~gens` methods.

        INPUT:

        - ``self`` -- an elliptic curve `E` over a number field `K`

        - ``verbose`` -- 0, 1, 2, or 3 (default: 0); the verbosity level

        - ``lim1`` -- (default: 2) limit on trivial points on quartics

        - ``lim3`` -- (default: 4) limit on points on ELS quartics

        - ``limtriv`` -- (default: 2) limit on trivial points on `E`

        - ``maxprob`` -- (default: 20)

        - ``limbigprime`` -- (default: 30) to distinguish between
          small and large prime numbers. Use probabilistic tests for
          large primes. If 0, do not use probabilistic tests.

        - ``known_points`` -- (default: ``None``) list of known points on
          the curve

        OUTPUT: a triple ``(lower, upper, list)`` consisting of

        - ``lower`` -- integer; lower bound on the rank

        - ``upper`` -- integer; upper bound on the rank

        - ``list`` -- list of points in `E(K)`

        The integer ``upper`` is in fact an upper bound on the
        dimension of the 2-Selmer group, hence on the dimension of
        `E(K)/2E(K)`.  It is equal to the dimension of the 2-Selmer
        group except possibly if `E(K)[2]` has dimension 1.  In that
        case, ``upper`` may exceed the dimension of the 2-Selmer group
        by an even number, due to the fact that the algorithm does not
        perform a second descent.

        .. NOTE::

            For non-quadratic number fields, this code does return, but
            it takes a long time.

        ALGORITHM:

        Uses Denis Simon's PARI/GP scripts from
        https://simond.users.lmno.cnrs.fr/.

        EXAMPLES::

            sage: x = polygen(ZZ, 'x')
            sage: K.<a> = NumberField(x^2 + 23, 'a')
            sage: E = EllipticCurve(K, '37')
            sage: E == loads(dumps(E))
            True
            sage: E.simon_two_descent()
            (2, 2, [(0 : 0 : 1), (1/8*a + 5/8 : -3/16*a - 7/16 : 1)])
            sage: E.simon_two_descent(lim1=5, lim3=5, limtriv=10, maxprob=7, limbigprime=10)
            (2, 2, [(-1 : 0 : 1), (1/2*a - 5/2 : -1/2*a - 13/2 : 1)])

        ::

            sage: K.<a> = NumberField(x^2 + 7, 'a')
            sage: E = EllipticCurve(K, [0,0,0,1,a]); E
            Elliptic Curve defined by y^2 = x^3 + x + a
             over Number Field in a with defining polynomial x^2 + 7

            sage: v = E.simon_two_descent(verbose=1); v
             elliptic curve: Y^2 = x^3 + Mod(1, y^2 + 7)*x + Mod(y, y^2 + 7)
             Trivial points on the curve = [[1, 1, 0], [Mod(1/2*y + 3/2, y^2 + 7), Mod(-y - 2, y^2 + 7), 1]]
            #S(E/K)[2]    = 2
            #E(K)/2E(K)   = 2
            #III(E/K)[2]  = 1
            rank(E/K)     = 1
             listpoints = [[Mod(1/2*y + 3/2, y^2 + 7), Mod(-y - 2, y^2 + 7), 1]]
            (1, 1, [(1/2*a + 3/2 : -a - 2 : 1)])

        A curve with 2-torsion::

            sage: K.<a> = NumberField(x^2 + 7)
            sage: E = EllipticCurve(K, '15a')
            sage: E.simon_two_descent()  # long time (3s on sage.math, 2013), points can vary
            (1, 3, [...])

        Check that the bug reported in :issue:`15483` is fixed::

            sage: K.<s> = QuadraticField(229)
            sage: c4 = 2173 - 235*(1 - s)/2
            sage: c6 = -124369 + 15988*(1 - s)/2
            sage: E = EllipticCurve([-c4/48, -c6/864])
            sage: E.simon_two_descent()
            (0, 0, [])

            sage: R.<t> = QQ[]
            sage: L.<g> = NumberField(t^3 - 9*t^2 + 13*t - 4)
            sage: E1 = EllipticCurve(L, [1-g*(g-1), -g^2*(g-1), -g^2*(g-1), 0, 0])
            sage: E1.rank()  # long time (about 5 s)
            0

            sage: K = CyclotomicField(43).subfields(3)[0][0]
            sage: E = EllipticCurve(K, '37')
            sage: E.simon_two_descent()  # long time (4s on sage.math, 2013)
            (3,
             3,
             [(-1/2*zeta43_0^2 - 1/2*zeta43_0 + 7 : -3/2*zeta43_0^2 - 5/2*zeta43_0 + 18 : 1),
              (0 : 0 : 1)])
        """
        verbose = int(verbose)
        if known_points is None:
            known_points = self._known_points
        known_points = [self(P) for P in known_points]

        # We deliberately do not use known_points as a key in the
        # following caching code, so that calling E.gens() a second
        # time (when known_points may have increased) will not cause
        # another execution of simon_two_descent.
        try:
            result = self._simon_two_descent_data[lim1,lim3,limtriv,maxprob,limbigprime]
            if verbose == 0:
                return result
        except AttributeError:
            self._simon_two_descent_data = {}
        except KeyError:
            pass

        from .gp_simon import simon_two_descent
        t = simon_two_descent(self, verbose=verbose,
                              lim1=lim1, lim3=lim3, limtriv=limtriv,
                              maxprob=maxprob, limbigprime=limbigprime,
                              known_points=known_points)
        self._simon_two_descent_data[lim1,lim3,limtriv,maxprob,limbigprime] = t
        self._known_points.extend([P for P in t[2]
                                   if P not in self._known_points])
        return t

    def height_pairing_matrix(self, points=None, precision=None, normalised=True):
        r"""
        Return the height pairing matrix of the given points.

        INPUT:

        - ``points`` -- list or ``None`` (default); list of points
          on this curve, or ``None``, in which case ``self.gens()`` will
          be used

        - ``precision`` -- integer or ``None`` (default); number of bits
          of precision of result, or ``None``, for default RealField precision

        - ``normalised`` -- boolean (default: ``True``); if ``True``, use
          normalised heights which are independent of base change.
          Otherwise use the non-normalised Néron-Tate height, as
          required for the regulator in the BSD conjecture.

        EXAMPLES::

            sage: E = EllipticCurve([0, 0, 1, -1, 0])
            sage: E.height_pairing_matrix()
            [0.0511114082399688]

        For rank 0 curves, the result is a valid 0x0 matrix::

            sage: EllipticCurve('11a').height_pairing_matrix()
            []
            sage: E = EllipticCurve('5077a1')
            sage: E.height_pairing_matrix([E.lift_x(x) for x in [-2,-7/4,1]], precision=100)
            [  1.3685725053539301120518194471  -1.3095767070865761992624519454 -0.63486715783715592064475542573]
            [ -1.3095767070865761992624519454   2.7173593928122930896610589220   1.0998184305667292139777571432]
            [-0.63486715783715592064475542573   1.0998184305667292139777571432  0.66820516565192793503314205089]

            sage: E = EllipticCurve('389a1')
            sage: E = EllipticCurve('389a1')
            sage: P, Q = E.point([-1,1,1]), E.point([0,-1,1])
            sage: E.height_pairing_matrix([P,Q])
            [0.686667083305587 0.268478098806726]
            [0.268478098806726 0.327000773651605]

        Over a number field::

            sage: x = polygen(QQ)
            sage: K.<t> = NumberField(x^2 + 47)
            sage: EK = E.base_extend(K)
            sage: EK.height_pairing_matrix([EK(P),EK(Q)])
            [0.686667083305587 0.268478098806726]
            [0.268478098806726 0.327000773651605]

        ::

            sage: K.<i> = QuadraticField(-1)
            sage: E = EllipticCurve([0,0,0,i,i])
            sage: P = E(-9+4*i, -18-25*i)
            sage: Q = E(i,-i)
            sage: E.height_pairing_matrix([P,Q])
            [  2.16941934493768 -0.870059380421505]
            [-0.870059380421505  0.424585837470709]
            sage: E.regulator_of_points([P,Q])
            0.164101403936070

        When the parameter ``normalised`` is set to ``False``, each
        height is multiplied by the degree `d` of the base field, and
        the regulator of `r` points is multiplied by `d^r`::

            sage: E.height_pairing_matrix([P,Q], normalised=False)
            [ 4.33883868987537 -1.74011876084301]
            [-1.74011876084301 0.849171674941418]
            sage: E.regulator_of_points([P,Q], normalised=False)
            0.656405615744281
        """
        if points is None:
            points = self.gens()
        else:
            for P in points:
                assert P.curve() == self

        r = len(points)
        if precision is None:
            RR = RealField()
        else:
            RR = RealField(precision)
        from sage.matrix.matrix_space import MatrixSpace
        M = MatrixSpace(RR, r)
        mat = M()
        for j in range(r):
            mat[j, j] = points[j].height(precision=precision, normalised=normalised)
        for j in range(r):
            for k in range(j + 1, r):
                mat[j, k] = ((points[j] + points[k]).height(precision=precision, normalised=normalised) - mat[j, j] - mat[k, k]) / 2
                mat[k, j] = mat[j, k]
        return mat

    def regulator_of_points(self, points=[], precision=None, normalised=True):
        """Return the regulator of the given points on this curve.

        INPUT:

        - ``points`` -- (default: empty list) a list of points on this curve

        - ``precision`` -- integer or ``None`` (default); the
          precision in bits of the result (default: real precision if ``None``)

        - ``normalised`` -- boolean (default: ``True``); if ``True``, use
          normalised heights which are independent of base change.
          Otherwise use the non-normalised Néron-Tate height, as
          required for the regulator in the BSD conjecture

        EXAMPLES::

            sage: E = EllipticCurve('37a1')
            sage: P = E(0,0)
            sage: Q = E(1,0)
            sage: E.regulator_of_points([P,Q])
            0.000000000000000
            sage: 2*P == Q
            True

        ::

            sage: E = EllipticCurve('5077a1')
            sage: points = [E.lift_x(x) for x in [-2,-7/4,1]]
            sage: E.regulator_of_points(points)
            0.417143558758384
            sage: E.regulator_of_points(points, precision=100)
            0.41714355875838396981711954462

        ::

            sage: E = EllipticCurve('389a')
            sage: E.regulator_of_points()
            1.00000000000000
            sage: points = [P,Q] = [E(-1,1), E(0,-1)]
            sage: E.regulator_of_points(points)
            0.152460177943144
            sage: E.regulator_of_points(points, precision=100)
            0.15246017794314375162432475705
            sage: E.regulator_of_points(points, precision=200)
            0.15246017794314375162432475704945582324372707748663081784028
            sage: E.regulator_of_points(points, precision=300)
            0.152460177943143751624324757049455823243727077486630817840280980046053225683562463604114816

        Examples over number fields::

            sage: K.<a> = QuadraticField(97)
            sage: E = EllipticCurve(K, [1,1])
            sage: P = E(0,1)
            sage: P.height()
            0.476223106404866
            sage: E.regulator_of_points([P])
            0.476223106404866

        When the parameter ``normalised`` is set to ``False``, each
        height is multiplied by the degree `d` of the base field, and
        the regulator of `r` points is multiplied by `d^r`::

            sage: P.height(normalised=False)
            0.952446212809731
            sage: E.regulator_of_points([P], normalised=False)
            0.952446212809731

        ::

            sage: E = EllipticCurve('11a1')
            sage: x = polygen(QQ)
            sage: K.<t> = NumberField(x^2 + 47)
            sage: EK = E.base_extend(K)
            sage: T = EK(5, 5)
            sage: T.order()
            5
            sage: P = EK(-2, -1/2*t - 1/2)
            sage: P.order()
            +Infinity
            sage: EK.regulator_of_points([P,T]) # random very small output
            -1.23259516440783e-32
            sage: EK.regulator_of_points([P,T]).abs() < 1e-30
            True

        ::

            sage: E = EllipticCurve('389a1')
            sage: P,Q = E.gens()
            sage: E.regulator_of_points([P,Q])
            0.152460177943144
            sage: K.<t> = NumberField(x^2 + 47)
            sage: EK = E.base_extend(K)
            sage: EK.regulator_of_points([EK(P),EK(Q)])
            0.152460177943144

        ::

            sage: K.<i> = QuadraticField(-1)
            sage: E = EllipticCurve([0,0,0,i,i])
            sage: P = E(-9+4*i, -18-25*i)
            sage: Q = E(i, -i)
            sage: E.height_pairing_matrix([P,Q])
            [  2.16941934493768 -0.870059380421505]
            [-0.870059380421505  0.424585837470709]
            sage: E.regulator_of_points([P,Q])
            0.164101403936070
        """
        if points is None:
            points = []
        mat = self.height_pairing_matrix(points=points, precision=precision, normalised=normalised)
        return mat.det(algorithm='hessenberg')

    def is_local_integral_model(self, *P):
        r"""
        Test if ``self`` is integral at the prime ideal `P`, or at all the
        primes if `P` is a list or tuple.

        INPUT:

        - ``*P`` -- a prime ideal, or a list or tuple of primes

        EXAMPLES::

            sage: K.<i> = NumberField(x^2 + 1)
            sage: P1, P2 = K.primes_above(5)
            sage: E = EllipticCurve([i/5,i/5,i/5,i/5,i/5])
            sage: E.is_local_integral_model(P1, P2)
            False
            sage: Emin = E.local_integral_model(P1, P2)
            sage: Emin.is_local_integral_model(P1, P2)
            True
        """
        if len(P) == 1:
            P = P[0]
        if isinstance(P, (tuple, list)):
            return all(self.is_local_integral_model(x) for x in P)
        return all(x.valuation(P) >= 0 for x in self.ainvs())

    def local_integral_model(self,*P):
        r"""
        Return a model of ``self`` which is integral at the prime ideal `P`.

        .. NOTE::

            The integrality at other primes is not affected, even if
            `P` is non-principal.

        INPUT:

        - ``*P`` -- a prime ideal, or a list or tuple of primes

        EXAMPLES::

            sage: K.<i> = NumberField(x^2 + 1)
            sage: P1, P2 = K.primes_above(5)
            sage: E = EllipticCurve([i/5,i/5,i/5,i/5,i/5])
            sage: E.local_integral_model((P1,P2))
            Elliptic Curve defined by
             y^2 + (-i)*x*y + (-25*i)*y = x^3 + 5*i*x^2 + 125*i*x + 3125*i
             over Number Field in i with defining polynomial x^2 + 1
        """
        if len(P) == 1:
            P = P[0]
        if isinstance(P, (tuple, list)):
            E = self
            for Pi in P:
                E = E.local_integral_model(Pi)
            return E
        ai = self.a_invariants()
        e = min((ai[i].valuation(P) / [1, 2, 3, 4, 6][i])
                for i in range(5)).floor()
        pi = self.base_field().uniformizer(P, 'negative')
        return EllipticCurve([ai[i]/pi**(e*[1,2,3,4,6][i]) for i in range(5)])

    def is_global_integral_model(self):
        r"""
        Return whether ``self`` is integral at all primes.

        EXAMPLES::

            sage: K.<i> = NumberField(x^2 + 1)
            sage: E = EllipticCurve([i/5,i/5,i/5,i/5,i/5])
            sage: P1, P2 = K.primes_above(5)
            sage: Emin = E.global_integral_model()
            sage: Emin.is_global_integral_model()
            True
        """
        return all(x.is_integral() for x in self.a_invariants())

    def global_integral_model(self):
        r"""
        Return a model of ``self`` which is integral at all primes.

        EXAMPLES::

            sage: K.<i> = NumberField(x^2 + 1)
            sage: E = EllipticCurve([i/5,i/5,i/5,i/5,i/5])
            sage: P1, P2 = K.primes_above(5)
            sage: E.global_integral_model()
            Elliptic Curve defined by
             y^2 + (-i)*x*y + (-25*i)*y = x^3 + 5*i*x^2 + 125*i*x + 3125*i
             over Number Field in i with defining polynomial x^2 + 1

        :issue:`7935`::

            sage: K.<a> = NumberField(x^2 - 38)
            sage: E = EllipticCurve([a,1/2])
            sage: E.global_integral_model()
            Elliptic Curve defined by y^2 = x^3 + 1444*a*x + 27436
             over Number Field in a with defining polynomial x^2 - 38

        :issue:`9266`::

            sage: K.<s> = NumberField(x^2 - 5)
            sage: w = (1+s)/2
            sage: E = EllipticCurve(K, [2,w])
            sage: E.global_integral_model()
            Elliptic Curve defined by y^2 = x^3 + 2*x + (1/2*s+1/2)
             over Number Field in s with defining polynomial x^2 - 5

        :issue:`12151`::

            sage: K.<v> = NumberField(x^2 + 161*x - 150)
            sage: E = EllipticCurve([25105/216*v - 3839/36, 634768555/7776*v - 98002625/1296, 634768555/7776*v - 98002625/1296, 0, 0])
            sage: M = E.global_integral_model(); M # choice varies, not tested
            Elliptic Curve defined by
            y^2 + (2094779518028859*v-1940492905300351)*x*y + (477997268472544193101178234454165304071127500*v-442791377441346852919930773849502871958097500)*y = x^3 + (26519784690047674853185542622500*v-24566525306469707225840460652500)*x^2
             over Number Field in v with defining polynomial x^2 + 161*x - 150

        :issue:`14476`::

            sage: R.<t> = QQ[]
            sage: K.<g> = NumberField(t^4 - t^3 - 3*t^2 - t + 1)
            sage: E = EllipticCurve([ -43/625*g^3 + 14/625*g^2 - 4/625*g + 706/625, -4862/78125*g^3 - 4074/78125*g^2 - 711/78125*g + 10304/78125,  -4862/78125*g^3 - 4074/78125*g^2 - 711/78125*g + 10304/78125, 0,0])
            sage: E.global_integral_model()
            Elliptic Curve defined by
             y^2 + (15*g^3-48*g-42)*x*y + (-111510*g^3-162162*g^2-44145*g+37638)*y = x^3 + (-954*g^3-1134*g^2+81*g+576)*x^2
             over Number Field in g with defining polynomial t^4 - t^3 - 3*t^2 - t + 1

        TESTS:

        Check the skipped test from above::

            sage: K.<v> = NumberField(x^2 + 161*x - 150)
            sage: E = EllipticCurve([25105/216*v - 3839/36, 634768555/7776*v - 98002625/1296, 634768555/7776*v - 98002625/1296, 0, 0])
            sage: M = E.global_integral_model()
            sage: b = M.ainvs()
            sage: b[0] in (2094779518028859*v-1940492905300351, 33872485050625*v - 31078224284250)
            True
            sage: b[1] in (26519784690047674853185542622500*v - 24566525306469707225840460652500,
            ....:          6933305282258321342920781250*v - 6422644400723486559914062500)
            True
            sage: b[2] in (477997268472544193101178234454165304071127500*v -442791377441346852919930773849502871958097500,
            ....:          2020602604156076340058146664245468750000*v - 1871778534673615560803175189398437500000)
            True
            sage: b[3:]
            (0, 0)
        """
        K = self.base_field()
        ai = self.a_invariants()
        Ps = {ff[0] for a in ai if not a.is_integral()
              for ff in a.denominator_ideal().factor()}
        for P in Ps:
            pi = K.uniformizer(P, 'positive')
            e = min((ai[i].valuation(P)/[1,2,3,4,6][i])
                    for i in range(5)).floor()
            if e < 0:
                ai = [ai[i]/pi**(e*[1,2,3,4,6][i]) for i in range(5)]
            if all(a.is_integral() for a in ai):
                break
        for z in ai:
            assert z.is_integral(), "bug in global_integral_model: %s" % list(ai)
        return EllipticCurve(list(ai))

    integral_model = global_integral_model

    def _reduce_model(self):
        r"""
        Return a reduced model for this elliptic curve.

        Transforms the elliptic curve to a model which is optimally scaled
        with respect to units and in which `a_1`, `a_2`, `a_3` are
        reduced modulo 2, 3, 2 respectively.

        .. NOTE::

            This only works on integral models, i.e. it requires that
            `a_1`, `a_2` and `a_3` lie in the ring of integers of the base
            field.

        EXAMPLES::

            sage: x = polygen(ZZ, 'x')
            sage: K.<a> = NumberField(x^2 - 38)
            sage: E = EllipticCurve([a, -5*a + 19, -39*a + 237, 368258520200522046806318224*a - 2270097978636731786720858047, 8456608930180227786550494643437985949781*a - 52130038506835491453281450568107193773505])
            sage: E.ainvs()
            (a,
            -5*a + 19,
            -39*a + 237,
            368258520200522046806318224*a - 2270097978636731786720858047,
            8456608930180227786550494643437985949781*a - 52130038506835491453281450568107193773505)
            sage: E._reduce_model().ainvs()
            (a,
             a + 1,
             a + 1,
             368258520200522046806318444*a - 2270097978636731786720859345,
             8456608930173478039472018047583706316424*a - 52130038506793883217874390501829588391299)
            sage: EllipticCurve([101,202,303,404,505])._reduce_model().ainvs()
            (1, 1, 0, -2509254, 1528863051)
            sage: EllipticCurve([-101,-202,-303,-404,-505])._reduce_model().ainvs()
            (1, -1, 0, -1823195, 947995262)

            sage: E = EllipticCurve([a/4, 1])
            sage: E._reduce_model()
            Traceback (most recent call last):
            ...
            TypeError: _reduce_model() requires an integral model.
        """
        K = self.base_ring()
        ZK = K.maximal_order()
        try:
            a1, a2, a3, a4, a6 = (ZK(a) for a in self.a_invariants())
        except TypeError:
            raise TypeError("_reduce_model() requires an integral model.")

        # N.B. Must define s, r, t in the right order.
        if ZK.absolute_degree() == 1:
            s = ((-a1) / 2).round('up')
            r = ((-a2 + s * a1 + s * s) / 3).round()
            t = ((-a3 - r * a1) / 2).round('up')
        else:
            pariK = K.__pari__()
            s = K(pariK.nfeltdiveuc(-a1, 2))
            r = K(pariK.nfeltdiveuc(-a2 + s * a1 + s * s, 3))
            t = K(pariK.nfeltdiveuc(-a3 - r * a1, 2))

        return self.rst_transform(r, s, t)

    def _scale_by_units(self):
        r"""
        Return a model reduced with respect to scaling by units.

        OUTPUT:

        A model for this elliptic curve, optimally scaled with respect
        to scaling by units, with respect to the logarithmic embedding
        of `|c4|^(1/4)+|c6|^(1/6)`. No scaling by roots of unity is
        carried out, so there is no change when the unit rank is 0.

        EXAMPLES::

            sage: x = polygen(ZZ, 'x')
            sage: K.<a> = NumberField(x^2 - 10)
            sage: u = a + 3
            sage: u.is_unit()
            True
            sage: E = EllipticCurve([0, 0, 0, 4536*a + 14148, -163728*a - 474336])
            sage: E1 = E.scale_curve(u^5)
            sage: E1.ainvs()
            (0,
            0,
            0,
            28087920796764302856*a + 88821804456186580548,
            -77225139016967233228487820912*a - 244207331916752959911655344864)
            sage: E1._scale_by_units().ainvs()
            (0, 0, 0, 4536*a + 14148, -163728*a - 474336)

        A totally real cubic example::

            sage: K.<a> = NumberField(x^3-x^2-6*x+5)
            sage: E = EllipticCurve([a + 1, a^2 + a - 1, a + 1, 44*a^2 + a - 258, -215*a^2 + 53*a + 1340])
            sage: u1, u2 = K.units()
            sage: u = u1^2/u2^3
            sage: E1 = E.scale_curve(u)
            sage: E1._scale_by_units().ainvs() == E.ainvs()
            True

        A complex quartic example::

            sage: K.<a> = CyclotomicField(5)
            sage: E = EllipticCurve([a + 1, a^2 + a - 1, a + 1, 44*a^2 + a - 258, -215*a^2 + 53*a + 1340])
            sage: u = K.units()[0]
            sage: E1 = E.scale_curve(u^5)
            sage: E1._scale_by_units().ainvs() == E.ainvs()
            True

        TESTS:

        See :issue:`34174`.  This used to raise an error due to insufficient precision::

            sage: K.<a> = QuadraticField(4569)
            sage: j = 46969655/32768
            sage: E = EllipticCurve(j=K(j))
            sage: C = E.isogeny_class()                                                 # long time (9.5s)
        """
        K = self.base_field()
        r1, r2 = K.signature()
        if r1 + r2 == 1:  # unit rank is 0
            return self

        degs = [1]*r1 + [2]*r2
        fu = K.units()
        c4, c6 = self.c_invariants()

        from sage.matrix.constructor import Matrix
        from sage.modules.free_module_element import vector

        prec = 1000 # initial value, will be increased if necessary
        ok = False
        while not ok:
            embs = K.places(prec=prec)
            c4s = [e(c4) for e in embs]
            c6s = [e(c6) for e in embs]

            U = Matrix([[e(u).abs().log()*d for d,e in zip(degs,embs)] for u in fu])
            v = vector([(x4.abs().nth_root(4)+x6.abs().nth_root(6)).log()*d for x4,x6,d in zip(c4s,c6s,degs)])
            w = -(U*U.transpose()).inverse()*U*v
            try:
                es = [e.round() for e in w]
                ok = True
            except ValueError:
                prec *= 2

        u = prod([uj**ej for uj,ej in zip(fu,es)])
        return self.scale_curve(u)

    def local_data(self, P=None, proof=None, algorithm='pari', globally=False):
        r"""
        Local data for this elliptic curve at the prime `P`.

        INPUT:

        - ``P`` -- either ``None``, a prime ideal of the base field of ``self``,
          or an element of the base field that generates a prime ideal

        - ``proof`` -- whether to only use provably correct methods
          (default controlled by global proof module).  Note that the
          proof module is number_field, not elliptic_curves, since the
          functions that actually need the flag are in number fields.

        - ``algorithm`` -- string (default: ``'pari'``); ignored unless the
          base field is `\QQ`.  If ``'pari'``, use the PARI C-library
          :pari:`ellglobalred` implementation of Tate's algorithm over
          `\QQ`. If ``'generic'``, use the general number field
          implementation.

        - ``globally`` -- whether the local algorithm uses global generators
          for the prime ideals. Default is False, which will not require any
          information about the class group. If ``True``, a generator for `P`
          will be used if `P` is principal. Otherwise, or if ``globally``
          is False, the minimal model returned will preserve integrality
          at other primes, but not minimality.

        OUTPUT:

        If `P` is specified, returns the ``EllipticCurveLocalData``
        object associated to the prime `P` for this curve.  Otherwise,
        returns a list of such objects, one for each prime `P` in the
        support of the discriminant of this model.

        .. NOTE::

            The model is not required to be integral on input.

        EXAMPLES::

            sage: K.<i> = NumberField(x^2 + 1)
            sage: E = EllipticCurve([1 + i, 0, 1, 0, 0])
            sage: E.local_data()
            [Local data at Fractional ideal (-2*i - 1):
               Reduction type: bad non-split multiplicative
               Local minimal model: Elliptic Curve defined by y^2 + (i+1)*x*y + y = x^3
                                    over Number Field in i with defining polynomial x^2 + 1
               Minimal discriminant valuation: 1
               Conductor exponent: 1
               Kodaira Symbol: I1
               Tamagawa Number: 1,
             Local data at Fractional ideal (3*i + 2):
               Reduction type: bad split multiplicative
               Local minimal model: Elliptic Curve defined by y^2 + (i+1)*x*y + y = x^3
                                    over Number Field in i with defining polynomial x^2 + 1
               Minimal discriminant valuation: 2
               Conductor exponent: 1
               Kodaira Symbol: I2
               Tamagawa Number: 2]
            sage: E.local_data(K.ideal(3))
            Local data at Fractional ideal (3):
              Reduction type: good
              Local minimal model: Elliptic Curve defined by y^2 + (i+1)*x*y + y = x^3
                                   over Number Field in i with defining polynomial x^2 + 1
              Minimal discriminant valuation: 0
              Conductor exponent: 0
              Kodaira Symbol: I0
              Tamagawa Number: 1
            sage: E.local_data(2*i + 1)
            Local data at Fractional ideal (-2*i - 1):
              Reduction type: bad non-split multiplicative
              Local minimal model: Elliptic Curve defined by y^2 + (i+1)*x*y + y = x^3
                                   over Number Field in i with defining polynomial x^2 + 1
              Minimal discriminant valuation: 1
              Conductor exponent: 1
              Kodaira Symbol: I1
              Tamagawa Number: 1

        An example raised in :issue:`3897`::

            sage: E = EllipticCurve([1,1])
            sage: E.local_data(3)
            Local data at Principal ideal (3) of Integer Ring:
              Reduction type: good
              Local minimal model: Elliptic Curve defined by y^2 = x^3 + x + 1
                                   over Rational Field
              Minimal discriminant valuation: 0
              Conductor exponent: 0
              Kodaira Symbol: I0
              Tamagawa Number: 1
        """
        if proof is None:
            import sage.structure.proof.proof
            # We use the "number_field" flag because the actual proof dependence is in PARI's number field functions.
            proof = sage.structure.proof.proof.get_flag(None, "number_field")

        if P is None:
            primes = self.base_ring()(self.integral_model().discriminant()).support()
            return [self._get_local_data(pr, proof) for pr in primes]

        from sage.schemes.elliptic_curves.ell_local_data import check_prime
        P = check_prime(self.base_field(),P)

        return self._get_local_data(P,proof,algorithm,globally)

    def _get_local_data(self, P, proof, algorithm='pari', globally=False):
        r"""
        Internal function to create data for this elliptic curve at the prime `P`.

        This function handles the caching of local data.  It is called
        by local_data() which is the user interface and which parses
        the input parameters `P` and proof.

        INPUT:

        - ``P`` -- either ``None`` or a prime ideal of the base field of ``self``

        - ``proof`` -- whether to only use provably correct methods
          (default controlled by global proof module).  Note that the
          proof module is number_field, not elliptic_curves, since the
          functions that actually need the flag are in number fields.

        - ``algorithm`` -- string (default: ``'pari'``); ignored unless the
          base field is `\QQ`.  If ``'pari'``, use the PARI C-library
          :pari:`ellglobalred` implementation of Tate's algorithm over
          `\QQ`. If ``'generic'``, use the general number field
          implementation.

        - ``globally`` -- whether the local algorithm uses global generators
          for the prime ideals. Default is False, which will not require any
          information about the class group. If ``True``, a generator for `P`
          will be used if `P` is principal. Otherwise, or if ``globally``
          is False, the minimal model returned will preserve integrality
          at other primes, but not minimality.

        EXAMPLES::

            sage: K.<i> = NumberField(x^2 + 1)
            sage: E = EllipticCurve(K, [0,1,0,-160,308])
            sage: p = K.ideal(i + 1)
            sage: E._get_local_data(p, False)
            Local data at Fractional ideal (i + 1):
              Reduction type: good
              Local minimal model: Elliptic Curve defined by
                                   y^2 + x*y + y = x^3 + x^2 + (-10)*x + (-10)
                                   over Number Field in i with defining polynomial x^2 + 1
              Minimal discriminant valuation: 0
              Conductor exponent: 0
              Kodaira Symbol: I0
              Tamagawa Number: 1

        Verify that we cache based on the proof value and the algorithm choice::

            sage: E._get_local_data(p, False) is E._get_local_data(p, True)
            False

            sage: E._get_local_data(p, None, "pari") is E._get_local_data(p, None, "generic")
            False
        """
        try:
            return self._local_data[P, proof, algorithm, globally]
        except AttributeError:
            self._local_data = {}
        except KeyError:
            pass
        from sage.schemes.elliptic_curves.ell_local_data import EllipticCurveLocalData
        self._local_data[P, proof, algorithm, globally] = EllipticCurveLocalData(self, P, proof, algorithm, globally)
        return self._local_data[P, proof, algorithm, globally]

    def local_minimal_model(self, P, proof=None, algorithm='pari'):
        r"""
        Return a model which is integral at all primes and minimal at `P`.

        INPUT:

        - ``P`` -- either ``None`` or a prime ideal of the base field of ``self``

        - ``proof`` -- whether to only use provably correct methods
          (default controlled by global proof module).  Note that the
          proof module is number_field, not elliptic_curves, since the
          functions that actually need the flag are in number fields.

        - ``algorithm`` -- string (default: ``'pari'``); ignored unless the
          base field is `\QQ`.  If ``'pari'``, use the PARI C-library
          :pari:`ellglobalred` implementation of Tate's algorithm over
          `\QQ`. If ``'generic'``, use the general number field
          implementation.

        OUTPUT:

        A model of the curve which is minimal (and integral) at `P`.

        .. NOTE::

            The model is not required to be integral on input.

            For principal `P`, a generator is used as a uniformizer,
            and integrality or minimality at other primes is not
            affected.  For non-principal `P`, the minimal model
            returned will preserve integrality at other primes, but not
            minimality.

        EXAMPLES::

            sage: K.<a> = NumberField(x^2 - 5)
            sage: E = EllipticCurve([20, 225, 750, 1250*a + 6250, 62500*a + 15625])
            sage: P = K.ideal(a)
            sage: E.local_minimal_model(P).ainvs()
            (0, 1, 0, 2*a - 34, -4*a + 66)
        """
        if proof is None:
            import sage.structure.proof.proof
            # We use the "number_field" flag because the actual proof dependence is in PARI's number field functions.
            proof = sage.structure.proof.proof.get_flag(None, "number_field")

        return self.local_data(P, proof, algorithm).minimal_model()

    def has_good_reduction(self, P):
        r"""
        Return ``True`` if this elliptic curve has good reduction at the prime `P`.

        INPUT:

        - ``P`` -- a prime ideal of the base field of ``self``, or a field
          element generating such an ideal

        OUTPUT:

        boolean; ``True`` if the curve has good reduction at `P`, else ``False``.

        .. NOTE::

            This requires determining a local integral minimal model;
            we do not just check that the discriminant of the current
            model has valuation zero.

        EXAMPLES::

            sage: E = EllipticCurve('14a1')
            sage: [(p, E.has_good_reduction(p)) for p in prime_range(15)]
            [(2, False), (3, True), (5, True), (7, False), (11, True), (13, True)]

            sage: K.<a> = NumberField(x^3 - 2)
            sage: P17a, P17b = [P for P,e in K.factor(17)]
            sage: E = EllipticCurve([0,0,0,0,2*a+1])
            sage: [(p, E.has_good_reduction(p)) for p in [P17a,P17b]]
            [(Fractional ideal (4*a^2 - 2*a + 1), True),
             (Fractional ideal (2*a + 1), False)]
        """
        return self.local_data(P).has_good_reduction()

    def has_bad_reduction(self, P):
        r"""
        Return ``True`` if this elliptic curve has bad reduction at the prime `P`.

        INPUT:

        - ``P`` -- a prime ideal of the base field of ``self``, or a field
          element generating such an ideal

        OUTPUT:

        boolean; ``True`` if the curve has bad reduction at `P`, else ``False``.

        .. NOTE::

            This requires determining a local integral minimal model;
            we do not just check that the discriminant of the current
            model has valuation zero.

        EXAMPLES::

            sage: E = EllipticCurve('14a1')
            sage: [(p, E.has_bad_reduction(p)) for p in prime_range(15)]
            [(2, True), (3, False), (5, False), (7, True), (11, False), (13, False)]

            sage: K.<a> = NumberField(x^3 - 2)
            sage: P17a, P17b = [P for P,e in K.factor(17)]
            sage: E = EllipticCurve([0,0,0,0,2*a+1])
            sage: [(p, E.has_bad_reduction(p)) for p in [P17a,P17b]]
            [(Fractional ideal (4*a^2 - 2*a + 1), False),
             (Fractional ideal (2*a + 1), True)]
        """
        return self.local_data(P).has_bad_reduction()

    def has_multiplicative_reduction(self, P):
        r"""
        Return ``True`` if this elliptic curve has (bad) multiplicative
        reduction at the prime `P`.

        .. NOTE::

            See also ``has_split_multiplicative_reduction()`` and
            ``has_nonsplit_multiplicative_reduction()``.

        INPUT:

        - ``P`` -- a prime ideal of the base field of ``self``, or a field
          element generating such an ideal

        OUTPUT:

        boolean; ``True`` if the curve has multiplicative reduction at `P`,
        else ``False``.

        EXAMPLES::

            sage: E = EllipticCurve('14a1')
            sage: [(p, E.has_multiplicative_reduction(p)) for p in prime_range(15)]
            [(2, True), (3, False), (5, False), (7, True), (11, False), (13, False)]

            sage: K.<a> = NumberField(x^3 - 2)
            sage: P17a, P17b = [P for P,e in K.factor(17)]
            sage: E = EllipticCurve([0,0,0,0,2*a+1])
            sage: [(p, E.has_multiplicative_reduction(p)) for p in [P17a,P17b]]
            [(Fractional ideal (4*a^2 - 2*a + 1), False),
             (Fractional ideal (2*a + 1), False)]
        """
        return self.local_data(P).has_multiplicative_reduction()

    def has_split_multiplicative_reduction(self, P):
        r"""
        Return ``True`` if this elliptic curve has (bad) split multiplicative reduction at the prime `P`.

        INPUT:

        - ``P`` -- a prime ideal of the base field of ``self``, or a field
          element generating such an ideal

        OUTPUT:

        boolean; ``True`` if the curve has split multiplicative reduction at
        `P`, else ``False``.

        EXAMPLES::

            sage: E = EllipticCurve('14a1')
            sage: [(p, E.has_split_multiplicative_reduction(p)) for p in prime_range(15)]
            [(2, False), (3, False), (5, False), (7, True), (11, False), (13, False)]

            sage: K.<a> = NumberField(x^3 - 2)
            sage: P17a, P17b = [P for P,e in K.factor(17)]
            sage: E = EllipticCurve([0,0,0,0,2*a+1])
            sage: [(p, E.has_split_multiplicative_reduction(p)) for p in [P17a,P17b]]
            [(Fractional ideal (4*a^2 - 2*a + 1), False),
             (Fractional ideal (2*a + 1), False)]
        """
        return self.local_data(P).has_split_multiplicative_reduction()

    def has_nonsplit_multiplicative_reduction(self, P):
        r"""
        Return ``True`` if this elliptic curve has (bad) non-split
        multiplicative reduction at the prime `P`.

        INPUT:

        - ``P`` -- a prime ideal of the base field of ``self``, or a field
          element generating such an ideal

        OUTPUT:

        boolean; ``True`` if the curve has non-split multiplicative
        reduction at `P`, else ``False``.

        EXAMPLES::

            sage: E = EllipticCurve('14a1')
            sage: [(p, E.has_nonsplit_multiplicative_reduction(p)) for p in prime_range(15)]
            [(2, True), (3, False), (5, False), (7, False), (11, False), (13, False)]

            sage: K.<a> = NumberField(x^3 - 2)
            sage: P17a, P17b = [P for P,e in K.factor(17)]
            sage: E = EllipticCurve([0,0,0,0,2*a+1])
            sage: [(p, E.has_nonsplit_multiplicative_reduction(p)) for p in [P17a,P17b]]
            [(Fractional ideal (4*a^2 - 2*a + 1), False),
             (Fractional ideal (2*a + 1), False)]
        """
        return self.local_data(P).has_nonsplit_multiplicative_reduction()

    def has_additive_reduction(self, P):
        r"""
        Return ``True`` if this elliptic curve has (bad) additive reduction at
        the prime `P`.

        INPUT:

        - ``P`` -- a prime ideal of the base field of ``self``, or a field
          element generating such an ideal

        OUTPUT:

        boolean; ``True`` if the curve has additive reduction at `P`, else
        ``False``.

        EXAMPLES::

            sage: E = EllipticCurve('27a1')
            sage: [(p, E.has_additive_reduction(p)) for p in prime_range(15)]
            [(2, False), (3, True), (5, False), (7, False), (11, False), (13, False)]

            sage: K.<a> = NumberField(x^3 - 2)
            sage: P17a, P17b = [P for P,e in K.factor(17)]
            sage: E = EllipticCurve([0,0,0,0,2*a+1])
            sage: [(p, E.has_additive_reduction(p)) for p in [P17a,P17b]]
            [(Fractional ideal (4*a^2 - 2*a + 1), False),
             (Fractional ideal (2*a + 1), True)]
        """
        return self.local_data(P).has_additive_reduction()

    def tamagawa_number(self, P, proof=None):
        r"""
        Return the Tamagawa number of this elliptic curve at the prime `P`.

        INPUT:

        - ``P`` -- either ``None`` or a prime ideal of the base field of ``self``

        - ``proof`` -- whether to only use provably correct methods
          (default controlled by global proof module).  Note that the
          proof module is number_field, not elliptic_curves, since the
          functions that actually need the flag are in number fields.

        OUTPUT: positive integer; the Tamagawa number of the curve at `P`

        EXAMPLES::

            sage: x = polygen(ZZ, 'x')
            sage: K.<a> = NumberField(x^2 - 5)
            sage: E = EllipticCurve([20, 225, 750, 625*a + 6875, 31250*a + 46875])
            sage: [E.tamagawa_number(P) for P in E.discriminant().support()]
            [1, 1, 1, 1]
            sage: K.<a> = QuadraticField(-11)
            sage: E = EllipticCurve('11a1').change_ring(K)
            sage: [E.tamagawa_number(P) for P in K(11).support()]
            [10]
        """
        if proof is None:
            import sage.structure.proof.proof
            # We use the "number_field" flag because the actual proof dependence is in PARI's number field functions.
            proof = sage.structure.proof.proof.get_flag(None, "number_field")

        return self.local_data(P, proof).tamagawa_number()

    def tamagawa_numbers(self):
        """
        Return a list of all Tamagawa numbers for all prime divisors of the
        conductor (in order).

        EXAMPLES::

            sage: e = EllipticCurve('30a1')
            sage: e.tamagawa_numbers()
            [2, 3, 1]
            sage: vector(e.tamagawa_numbers())
            (2, 3, 1)
            sage: K.<a> = NumberField(x^2 + 3)
            sage: eK = e.base_extend(K)
            sage: eK.tamagawa_numbers()
            [4, 6, 1]
        """
        return [self.tamagawa_number(p)
                for p in self.conductor().prime_factors()]

    def tamagawa_exponent(self, P, proof=None):
        r"""
        Return the Tamagawa index of this elliptic curve at the prime `P`.

        INPUT:

        - ``P`` -- either ``None`` or a prime ideal of the base field of ``self``

        - ``proof`` -- whether to only use provably correct methods
          (default controlled by global proof module).  Note that the
          proof module is number_field, not elliptic_curves, since the
          functions that actually need the flag are in number fields.

        OUTPUT: positive integer; the Tamagawa index of the curve at P

        EXAMPLES::

            sage: K.<a> = NumberField(x^2 - 5)
            sage: E = EllipticCurve([20, 225, 750, 625*a + 6875, 31250*a + 46875])
            sage: [E.tamagawa_exponent(P) for P in E.discriminant().support()]
            [1, 1, 1, 1]
            sage: K.<a> = QuadraticField(-11)
            sage: E = EllipticCurve('11a1').change_ring(K)
            sage: [E.tamagawa_exponent(P) for P in K(11).support()]
            [10]
        """
        if proof is None:
            import sage.structure.proof.proof
            # We use the "number_field" flag because the actual proof dependence is in PARI's number field functions.
            proof = sage.structure.proof.proof.get_flag(None, "number_field")

        return self.local_data(P, proof).tamagawa_exponent()

    def tamagawa_product(self):
        r"""Return the product of the Tamagawa numbers `c_v` where `v` runs
        over all prime ideals of `K`.

        .. NOTE::

            See also tamagawa_product_bsd(), which includes an
            additional factor when the model is not globally minimal,
            as required by the BSD formula.

        OUTPUT: a positive integer

        EXAMPLES::

            sage: K.<i> = NumberField(x^2 + 1)
            sage: E = EllipticCurve([0, 2+i])
            sage: E.tamagawa_product()
            1

            sage: E = EllipticCurve([(2*i+1)^2, i*(2*i+1)^7])
            sage: E.tamagawa_product()
            4

        An example over `\QQ`::

            sage: E = EllipticCurve('30a')
            sage: E.tamagawa_product()
            6

        An example with everywhere good reduction, where the product
        is empty::

            sage: x = polygen(QQ)
            sage: K.<a> = NumberField(x^2 - 38)
            sage: E = EllipticCurve( [a, -a + 1, a + 1, -5*a + 15, -5*a + 21])
            sage: E.tamagawa_numbers()
            []
            sage: E.tamagawa_product()
            1
        """
        return prod([ld.tamagawa_number() for ld in self.local_data()], Integer(1))

    def tamagawa_product_bsd(self):
        r"""Given an elliptic curve `E` over a number field `K`, this function
        returns the integer `C(E/K)` that appears in the Birch and
        Swinnerton-Dyer conjecture accounting for the local
        information at finite places. If the model is a global minimal
        model then `C(E/K)` is simply the product of the Tamagawa
        numbers `c_v` where `v` runs over all prime ideals of
        `K`. Otherwise, if the model has to be changed at a place `v`
        a correction factor appears.  The definition is such that
        `C(E/K)` times the periods at the infinite places is invariant
        under change of the Weierstrass model. See [Tate1966]_ and
        [DD2010]_ for details.

        .. NOTE::

            This definition differs from the definition of
            ``tamagawa_product`` for curves defined over `\QQ`. Over
            the rational number it is always defined to be the product
            of the Tamagawa numbers, so the two definitions only agree
            when the model is global minimal.

        OUTPUT: a rational number

        EXAMPLES::

            sage: K.<i> = NumberField(x^2 + 1)
            sage: E = EllipticCurve([0, 2+i])
            sage: E.tamagawa_product_bsd()
            1

            sage: E = EllipticCurve([(2*i+1)^2, i*(2*i+1)^7])
            sage: E.tamagawa_product_bsd()
            4

        An example where the Neron model changes over K::

            sage: K.<t> = NumberField(x^5 - 10*x^3 + 5*x^2 + 10*x + 1)
            sage: E = EllipticCurve(K, '75a1')
            sage: E.tamagawa_product_bsd()
            5
            sage: da = E.local_data()
            sage: [dav.tamagawa_number() for dav in da]
            [1, 1]

        An example over `\QQ` (:issue:`9413`)::

            sage: E = EllipticCurve('30a')
            sage: E.tamagawa_product_bsd()
            6
        """
        pr = QQ(1)
        for v in self.local_data():
            pp = v.prime()
            cv = v.tamagawa_number()
            # uu is the quotient of the Neron differential at pp divided by
            # the differential associated to this particular model of E
            uu = self.isomorphism_to(v.minimal_model()).u
            if self.base_field().absolute_degree() == 1:
                p = pp.gens_reduced()[0]
                f = 1
                v = ZZ(uu).valuation(p)
            else:
                p = pp.smallest_integer()
                f = pp.residue_class_degree()
                v = uu.valuation(pp)
            uu_abs_val = p**(f*v)
            pr *= cv * uu_abs_val
        return pr

    def kodaira_symbol(self, P, proof=None):
        r"""
        Return the Kodaira Symbol of this elliptic curve at the prime `P`.

        INPUT:

        - ``P`` -- either ``None`` or a prime ideal of the base field of ``self``

        - ``proof`` -- whether to only use provably correct methods
          (default controlled by global proof module).  Note that the
          proof module is number_field, not elliptic_curves, since the
          functions that actually need the flag are in number fields.

        OUTPUT: the Kodaira Symbol of the curve at ``P``, represented as a string

        EXAMPLES::

            sage: K.<a> = NumberField(x^2 - 5)
            sage: E = EllipticCurve([20, 225, 750, 625*a + 6875, 31250*a + 46875])
            sage: bad_primes = E.discriminant().support(); bad_primes
            [Fractional ideal (-a),
             Fractional ideal (-7/2*a + 81/2),
             Fractional ideal (-a - 52),
             Fractional ideal (2)]
            sage: [E.kodaira_symbol(P) for P in bad_primes]
            [I0, I1, I1, II]
            sage: K.<a> = QuadraticField(-11)
            sage: E = EllipticCurve('11a1').change_ring(K)
            sage: [E.kodaira_symbol(P) for P in K(11).support()]
            [I10]
        """
        if proof is None:
            import sage.structure.proof.proof
            # We use the "number_field" flag because the actual proof dependence is in PARI's number field functions.
            proof = sage.structure.proof.proof.get_flag(None, "number_field")

        return self.local_data(P, proof).kodaira_symbol()

    def conductor(self):
        r"""
        Return the conductor of this elliptic curve as a fractional
        ideal of the base field.

        OUTPUT: fractional ideal; the conductor of the curve

        EXAMPLES::

            sage: K.<i> = NumberField(x^2 + 1)
            sage: EllipticCurve([i, i - 1, i + 1, 24*i + 15, 14*i + 35]).conductor()
            Fractional ideal (3*i + 21)
            sage: K.<a> = NumberField(x^2 - x + 3)
            sage: cond = EllipticCurve([1 + a, -1 + a, 1 + a, -11 + a, 5 - 9*a]).conductor(); cond # random
            Fractional ideal (6*a)
            sage: cond == K.ideal(6*a)
            True

        A not so well known curve with everywhere good reduction::

            sage: K.<a> = NumberField(x^2 - 38)
            sage: E = EllipticCurve([0,0,0, 21796814856932765568243810*a - 134364590724198567128296995, 121774567239345229314269094644186997594*a - 750668847495706904791115375024037711300])
            sage: E.conductor()
            Fractional ideal (1)

        An example which used to fail (see :issue:`5307`)::

            sage: K.<w> = NumberField(x^2 + x + 6)
            sage: E = EllipticCurve([w, -1, 0, -w-6, 0])
            sage: E.conductor()
            Fractional ideal (86304, w + 5898)

        An example raised in :issue:`11346`::

            sage: K.<g> = NumberField(x^2 - x - 1)
            sage: E1 = EllipticCurve(K, [0, 0, 0, -1/48, -161/864])
            sage: [(p.smallest_integer(), e) for p,e in E1.conductor().factor()]
            [(2, 4), (3, 1), (5, 1)]
        """
        try:
            return self._conductor
        except AttributeError:
            pass

        # Note: for number fields other than QQ we could initialize
        # N=K.ideal(1) or N=OK.ideal(1), which are the same, but for
        # K == QQ it has to be ZZ.ideal(1).
        K = self.base_field()
        N = ZZ.ideal(1) if K is QQ else K.fractional_ideal(1)
        self._conductor = prod([d.prime()**d.conductor_valuation()
                                for d in self.local_data()],
                               N)
        return self._conductor

    def minimal_discriminant_ideal(self):
        r"""
        Return the minimal discriminant ideal of this elliptic curve.

        OUTPUT:

        The integral ideal `D` whose valuation at every prime `P` is
        that of the local minimal model for `E` at `P`.  If `E` has a
        global minimal model, this will be the principal ideal
        generated by the discriminant of any such model, but otherwise
        it can be a proper divisor of the discriminant of any model.

        EXAMPLES::

            sage: K.<a> = NumberField(x^2 - x - 57)
            sage: K.class_number()
            3
            sage: E = EllipticCurve([a, -a, a, -5692-820*a, -259213-36720*a])
            sage: K.ideal(E.discriminant())
            Fractional ideal (90118662980*a + 636812084644)
            sage: K.ideal(E.discriminant()).factor()
            (Fractional ideal (2))^2 * (Fractional ideal (3, a + 2))^12

        Here the minimal discriminant ideal is principal but there is
        no global minimal model since the quotient is the 12th power
        of a non-principal ideal::

            sage: E.minimal_discriminant_ideal()
            Fractional ideal (4)
            sage: E.minimal_discriminant_ideal().factor()
            (Fractional ideal (2))^2

        If (and only if) the curve has everywhere good reduction the
        result is the unit ideal::

            sage: K.<a> = NumberField(x^2 - 26)
            sage: E = EllipticCurve([a, a-1, a+1, 4*a+10, 2*a+6])
            sage: E.conductor()
            Fractional ideal (1)
            sage: E.discriminant()
            -104030*a - 530451
            sage: E.minimal_discriminant_ideal()
            Fractional ideal (1)

        Over `\QQ`, the result returned is an ideal of `\ZZ` rather
        than a fractional ideal of `\QQ`::

            sage: E = EllipticCurve([1,2,3,4,5])
            sage: E.minimal_discriminant_ideal()
            Principal ideal (10351) of Integer Ring
        """
        dat = self.local_data()
        # we treat separately the case where there are no bad primes,
        # which cannot happen over QQ, since ideals of QQ behave
        # differently to (fractional) ideals of other number fields.
        if not dat:
            return self.base_field().ideal(1)
        return prod([d.prime()**d.discriminant_valuation() for d in dat])

    def non_minimal_primes(self):
        r"""
        Return a list of primes at which this elliptic curve is not minimal.

        OUTPUT:

        A list of prime ideals (or prime numbers when the base field
        is `\QQ`, empty if this is a global minimal model.

        EXAMPLES::

            sage: K.<a> = NumberField(x^2 - 10)
            sage: E = EllipticCurve([0, 0, 0, -22500, 750000*a])
            sage: E.non_minimal_primes()
            [Fractional ideal (2, a), Fractional ideal (5, a)]
            sage: K.ideal(E.discriminant()).factor()
            (Fractional ideal (2, a))^24 * (Fractional ideal (3, a + 1))^5 * (Fractional ideal (3, a + 2))^5 * (Fractional ideal (5, a))^24 * (Fractional ideal (7))
            sage: E.minimal_discriminant_ideal().factor()
            (Fractional ideal (2, a))^12 * (Fractional ideal (3, a + 1))^5 * (Fractional ideal (3, a + 2))^5 * (Fractional ideal (7))

        Over `\QQ`, the primes returned are integers, not ideals::

            sage: E = EllipticCurve([0, 0, 0, -3024, 46224])
            sage: E.non_minimal_primes()
            [2, 3]
            sage: Emin = E.global_minimal_model()
            sage: Emin.non_minimal_primes()
            []

        If the model is not globally integral, a :exc:`ValueError` is
        raised::

            sage: E = EllipticCurve([0, 0, 0, 1/2, 1/3])
            sage: E.non_minimal_primes()
            Traceback (most recent call last):
            ...
            ValueError: non_minimal_primes only defined for integral models
        """
        if not self.is_global_integral_model():
            raise ValueError("non_minimal_primes only defined for integral models")
        dat = self.local_data()
        D = self.discriminant()
        primes = [d.prime() for d in dat]
        if self.base_field() is QQ:
            primes = [P.gen() for P in primes]
        vals = [d.discriminant_valuation() for d in dat]
        return [P for P,v in zip(primes,vals) if D.valuation(P) > v]

    def is_global_minimal_model(self):
        r"""
        Return whether this elliptic curve is a global minimal model.

        OUTPUT:

        Boolean, False if E is not integral, or if E is non-minimal at
        some prime, else True.

        EXAMPLES::

            sage: K.<a> = NumberField(x^2 - 10)
            sage: E = EllipticCurve([0, 0, 0, -22500, 750000*a])
            sage: E.is_global_minimal_model()
            False
            sage: E.non_minimal_primes()
            [Fractional ideal (2, a), Fractional ideal (5, a)]

            sage: E = EllipticCurve([0, 0, 0, -3024, 46224])
            sage: E.is_global_minimal_model()
            False
            sage: E.non_minimal_primes()
            [2, 3]
            sage: Emin = E.global_minimal_model()
            sage: Emin.is_global_minimal_model()
            True

        A necessary condition to be a global minimal model is that the
        model must be globally integral::

            sage: E = EllipticCurve([0,0,0,1/2,1/3])
            sage: E.is_global_minimal_model()
            False
            sage: Emin.is_global_minimal_model()
            True
            sage: Emin.ainvs()
            (0, 1, 1, -2, 0)
        """
        if not self.is_global_integral_model():
            return False
        return self.non_minimal_primes() == []

    def global_minimality_class(self):
        r"""
        Return the obstruction to this curve having a global minimal model.

        OUTPUT:

        An ideal class of the base number field, which is trivial if
        and only if the elliptic curve has a global minimal model, and
        which can be used to find global and semi-global minimal
        models.

        EXAMPLES:

        A curve defined over a field of class number 2 with no global
        minimal model was a nontrivial minimality class::

            sage: K.<a> = NumberField(x^2 - 10)
            sage: K.class_number()
            2
            sage: E = EllipticCurve([0, 0, 0, -22500, 750000*a])
            sage: E.global_minimality_class()
            Fractional ideal class (10, 5*a)
            sage: E.global_minimality_class().order()
            2

        Over the same field, a curve defined by a non-minimal model
        has trivial class, showing that a global minimal model does
        exist::

            sage: K.<a> = NumberField(x^2 - 10)
            sage: E = EllipticCurve([0, 0, 0, 4536*a+14148, -163728*a-474336])
            sage: E.is_global_minimal_model()
            False
            sage: E.global_minimality_class()
            Trivial principal fractional ideal class

        Over a field of class number 1 the result is always the
        trivial class::

            sage: K.<a> = NumberField(x^2 - 5)
            sage: E = EllipticCurve([0, 0, 0, K(16), K(64)])
            sage: E.global_minimality_class()
            Trivial principal fractional ideal class

            sage: E = EllipticCurve([0, 0, 0, 16, 64])
            sage: E.base_field()
            Rational Field
            sage: E.global_minimality_class()
            1
        """
        K = self.base_field()
        Cl = K.class_group()
        if K.class_number() == 1:
            return Cl(1)
        D = self.discriminant()
        dat = self.local_data()
        primes = [d.prime() for d in dat]
        vals = [d.discriminant_valuation() for d in dat]
        I = prod([P**((D.valuation(P)-v)//12) for P,v in zip(primes,vals)],
                 K.ideal(1))
        return Cl(I)

    def has_global_minimal_model(self):
        r"""
        Return whether this elliptic curve has a global minimal model.

        OUTPUT:

        boolean; ``True`` iff a global minimal model exists, i.e. an
        integral model which is minimal at every prime.

        EXAMPLES::

            sage: K.<a> = NumberField(x^2 - 10)
            sage: E = EllipticCurve([0, 0, 0, 4536*a+14148, -163728*a-474336])
            sage: E.is_global_minimal_model()
            False
            sage: E.has_global_minimal_model()
            True
        """
        return self.global_minimality_class().is_one()

    def global_minimal_model(self, proof=None, semi_global=False):
        r"""
        Return a model of ``self`` that is integral, and minimal.

        .. NOTE::

            Over fields of class number greater than 1, a global
            minimal model may not exist.  If it does not, set the
            parameter ``semi_global`` to ``True`` to obtain a model
            minimal at all but one prime.

        INPUT:

        - ``proof`` -- whether to only use provably correct methods
          (default controlled by global proof module).  Note that the
          proof module is number_field, not elliptic_curves, since the
          functions that actually need the flag are in number fields.

        - ``semi_global``-- boolean (default: ``False``); if there is no
          global minimal mode, return a semi-global minimal model
          (minimal at all but one prime) instead, if True; raise an
          error if False.  No effect if a global minimal model exists.

        OUTPUT:

        A global integral and minimal model, or an integral model
        minimal at all but one prime of there is no global minimal
        model and the flag ``semi_global`` is True.

        EXAMPLES::

            sage: K.<a> = NumberField(x^2 - 38)
            sage: E = EllipticCurve([0,0,0, 21796814856932765568243810*a - 134364590724198567128296995, 121774567239345229314269094644186997594*a - 750668847495706904791115375024037711300])
            sage: E2 = E.global_minimal_model()
            sage: E2
            Elliptic Curve defined by
             y^2 + a*x*y + (a+1)*y = x^3 + (a+1)*x^2 + (4*a+15)*x + (4*a+21)
             over Number Field in a with defining polynomial x^2 - 38
            sage: E2.local_data()
            []

        See :issue:`11347`::

            sage: K.<g> = NumberField(x^2 - x - 1)
            sage: E = EllipticCurve(K, [0, 0, 0, -1/48, 161/864])
            sage: E2 = E.integral_model().global_minimal_model(); E2
            Elliptic Curve defined by y^2 + x*y + y = x^3 + x^2
             over Number Field in g with defining polynomial x^2 - x - 1
            sage: [(p.norm(), e) for p, e in E2.conductor().factor()]
            [(9, 1), (5, 1)]
            sage: [(p.norm(), e) for p, e in E2.discriminant().factor()]
            [(-5, 2), (9, 1)]

        See :issue:`14472`, this used not to work over a relative extension::

            sage: K1.<w> = NumberField(x^2 + x + 1)
            sage: m = polygen(K1)
            sage: K2.<v> = K1.extension(m^2 - w + 1)
            sage: E = EllipticCurve([0*v, -432])
            sage: E.global_minimal_model()
            Elliptic Curve defined by y^2 + y = x^3
             over Number Field in v with defining polynomial x^2 - w + 1 over its base field

        See :issue:`18662`: for fields of class number greater than 1,
        even when global minimal models did exist, their computation
        was not implemented.  Now it is::

            sage: K.<a> = NumberField(x^2 - 10)
            sage: K.class_number()
            2
            sage: E = EllipticCurve([0, 0, 0, -186408*a - 589491, 78055704*a + 246833838])
            sage: E.discriminant().norm()
            16375845905239507992576
            sage: E.discriminant().norm().factor()
            2^31 * 3^27
            sage: E.has_global_minimal_model()
            True
            sage: Emin = E.global_minimal_model(); Emin
            Elliptic Curve defined by
             y^2 + (a+1)*x*y + (a+1)*y = x^3 + (-a)*x^2 + (a-12)*x + (-2*a+2)
             over Number Field in a with defining polynomial x^2 - 10
            sage: Emin.discriminant().norm()
            3456
            sage: Emin.discriminant().norm().factor()
            2^7 * 3^3

        If there is no global minimal model, this method will raise an
        error unless you set the parameter ``semi_global`` to ``True``::

            sage: K.<a> = NumberField(x^2 - 10)
            sage: K.class_number()
            2
            sage: E = EllipticCurve([a, a, 0, 3*a+8, 4*a+3])
            sage: E.has_global_minimal_model()
            False
            sage: E.global_minimal_model()
            Traceback (most recent call last):
            ...
            ValueError: Elliptic Curve defined by
            y^2 + a*x*y = x^3 + a*x^2 + (3*a+8)*x + (4*a+3) over Number Field
            in a with defining polynomial x^2 - 10 has no global minimal model!
            For a semi-global minimal model use semi_global=True
            sage: E.global_minimal_model(semi_global=True)
            Elliptic Curve defined by
             y^2 + a*x*y = x^3 + a*x^2 + (3*a+8)*x + (4*a+3) over Number Field in a
             with defining polynomial x^2 - 10

        An example of a curve with everywhere good reduction but which
        has no model with unit discriminant::

            sage: K.<a> = NumberField(x^2 - x - 16)
            sage: K.class_number()
            2
            sage: E = EllipticCurve([0,0,0,-15221331*a - 53748576, -79617688290*a - 281140318368])
            sage: Emin = E.global_minimal_model(semi_global=True)
            sage: Emin.ainvs()
            (a, a - 1, a, 605*a - 2728, 15887*a - 71972)
            sage: Emin.discriminant()
            -17*a - 16
            sage: Emin.discriminant().norm()
            -4096
            sage: Emin.minimal_discriminant_ideal()
            Fractional ideal (1)
            sage: E.conductor()
            Fractional ideal (1)
        """
        if proof is None:
            import sage.structure.proof.proof
            # We use the "number_field" flag because the actual proof dependence is in PARI's number field functions.
            proof = sage.structure.proof.proof.get_flag(None, "number_field")

        if self.has_global_minimal_model() or semi_global:
            if self.base_ring().class_number() == 1:
                E = self.global_integral_model()
                for P in E.base_ring()(E.discriminant()).support():
                    E = E.local_data(P,proof, globally=True).minimal_model()
            else:
                from .kraus import semi_global_minimal_model
                E, P = semi_global_minimal_model(self)
            return E._scale_by_units()._reduce_model()

        raise ValueError("%s has no global minimal model!  For a semi-global minimal model use semi_global=True" % self)

    def reduction(self,place):
        r"""
        Return the reduction of the elliptic curve at a place of good reduction.

        INPUT:

        - ``place`` -- a prime ideal in the base field of the curve

        OUTPUT: an elliptic curve over a finite field, the residue field of the place

        EXAMPLES::

            sage: K.<i> = QuadraticField(-1)
            sage: EK = EllipticCurve([0, 0, 0, i, i+3])
            sage: v = K.fractional_ideal(2*i+3)
            sage: EK.reduction(v)
            Elliptic Curve defined by y^2  = x^3 + 5*x + 8
             over Residue field of Fractional ideal (2*i + 3)
            sage: EK.reduction(K.ideal(1+i))
            Traceback (most recent call last):
            ...
            ValueError: The curve must have good reduction at the place.
            sage: EK.reduction(K.ideal(2))
            Traceback (most recent call last):
            ...
            ValueError: The ideal must be prime.
            sage: K = QQ.extension(x^2 + x + 1, "a")
            sage: E = EllipticCurve([1024*K.0, 1024*K.0])
            sage: E.reduction(2*K)
            Elliptic Curve defined by y^2 + (abar+1)*y = x^3
             over Residue field in abar of Fractional ideal (2)
        """
        K = self.base_field()
        OK = K.ring_of_integers()
        try:
            place = K.ideal(place)
        except TypeError:
            raise TypeError("The parameter must be an ideal of the base field of the elliptic curve")
        if not place.is_prime():
            raise ValueError("The ideal must be prime.")
        disc = self.discriminant()
        if not K.ideal(disc).valuation(place) == 0:
            local_data = self.local_data(place)
            if local_data.has_good_reduction():
                Fv = OK.residue_field(place)
                return local_data.minimal_model().change_ring(Fv)
            raise ValueError("The curve must have good reduction at the place.")
        Fv = OK.residue_field(place)
        return self.change_ring(Fv)

    @cached_method
    def torsion_subgroup(self):
        r"""
        Return the torsion subgroup of this elliptic curve.

        OUTPUT: the :class:`EllipticCurveTorsionSubgroup` associated to this elliptic
        curve.

        EXAMPLES::

            sage: E = EllipticCurve('11a1')
            sage: x = polygen(ZZ, 'x')
            sage: K.<t> = NumberField(x^4 + x^3 + 11*x^2 + 41*x + 101)
            sage: EK = E.base_extend(K)
            sage: tor = EK.torsion_subgroup()  # long time (2s on sage.math, 2014)
            sage: tor  # long time
            Torsion Subgroup isomorphic to Z/5 + Z/5 associated to the Elliptic Curve
             defined by y^2 + y = x^3 + (-1)*x^2 + (-10)*x + (-20) over Number Field
             in t with defining polynomial x^4 + x^3 + 11*x^2 + 41*x + 101
            sage: tor.gens()  # long time
            ((16 : 60 : 1), (t : 1/11*t^3 + 6/11*t^2 + 19/11*t + 48/11 : 1))

        ::

            sage: E = EllipticCurve('15a1')
            sage: K.<t> = NumberField(x^2 + 2*x + 10)
            sage: EK = E.base_extend(K)
            sage: EK.torsion_subgroup()
            Torsion Subgroup isomorphic to Z/4 + Z/4 associated to the
             Elliptic Curve defined by y^2 + x*y + y = x^3 + x^2 + (-10)*x + (-10)
             over Number Field in t with defining polynomial x^2 + 2*x + 10

        ::

            sage: E = EllipticCurve('19a1')
            sage: K.<t> = NumberField(x^9-3*x^8-4*x^7+16*x^6-3*x^5-21*x^4+5*x^3+7*x^2-7*x+1)
            sage: EK = E.base_extend(K)
            sage: EK.torsion_subgroup()
            Torsion Subgroup isomorphic to Z/9 associated to the Elliptic Curve defined
             by y^2 + y = x^3 + x^2 + (-9)*x + (-15) over Number Field in t with defining
             polynomial x^9 - 3*x^8 - 4*x^7 + 16*x^6 - 3*x^5 - 21*x^4 + 5*x^3 + 7*x^2 - 7*x + 1

        ::

            sage: K.<i> = QuadraticField(-1)
            sage: EK = EllipticCurve([0, 0, 0, i, i+3])
            sage: EK.torsion_subgroup ()
            Torsion Subgroup isomorphic to Trivial group associated to the
             Elliptic Curve defined by y^2 = x^3 + i*x + (i+3)
             over Number Field in i with defining polynomial x^2 + 1 with i = 1*I

        .. SEEALSO::

            Use :meth:`~sage.schemes.elliptic_curves.ell_field.EllipticCurve_field.division_field`
            to determine the field of definition of the `\ell`-torsion subgroup.
        """
        from .ell_torsion import EllipticCurveTorsionSubgroup
        return EllipticCurveTorsionSubgroup(self)

    @cached_method
    def torsion_order(self):
        r"""
        Return the order of the torsion subgroup of this elliptic curve.

        OUTPUT: integer

        EXAMPLES::

            sage: E = EllipticCurve('11a1')
            sage: x = polygen(ZZ, 'x')
            sage: K.<t> = NumberField(x^4 + x^3 + 11*x^2 + 41*x + 101)
            sage: EK = E.base_extend(K)
            sage: EK.torsion_order()  # long time (2s on sage.math, 2014)
            25

        ::

            sage: E = EllipticCurve('15a1')
            sage: K.<t> = NumberField(x^2 + 2*x + 10)
            sage: EK = E.base_extend(K)
            sage: EK.torsion_order()
            16

        ::

            sage: E = EllipticCurve('19a1')
            sage: K.<t> = NumberField(x^9 - 3*x^8 - 4*x^7 + 16*x^6 - 3*x^5 - 21*x^4 + 5*x^3 + 7*x^2 - 7*x + 1)
            sage: EK = E.base_extend(K)
            sage: EK.torsion_order()
            9

        ::

            sage: K.<i> = QuadraticField(-1)
            sage: EK = EllipticCurve([0, 0, 0, i, i + 3])
            sage: EK.torsion_order()
            1
        """
        return self.torsion_subgroup().order()

    def torsion_points(self):
        r"""
        Return a list of the torsion points of this elliptic curve.

        OUTPUT: sorted list of the torsion points

        EXAMPLES::

            sage: E = EllipticCurve('11a1')
            sage: E.torsion_points()
            [(0 : 1 : 0), (5 : -6 : 1), (5 : 5 : 1), (16 : -61 : 1), (16 : 60 : 1)]
            sage: x = polygen(ZZ, 'x')
            sage: K.<t> = NumberField(x^4 + x^3 + 11*x^2 + 41*x + 101)
            sage: EK = E.base_extend(K)
            sage: EK.torsion_points()  # long time (1s on sage.math, 2014)
            [(0 : 1 : 0),
             (t : 1/11*t^3 + 6/11*t^2 + 19/11*t + 48/11 : 1),
             (1/11*t^3 - 5/11*t^2 + 19/11*t - 40/11 : -6/11*t^3 - 3/11*t^2 - 26/11*t - 321/11 : 1),
             (1/11*t^3 - 5/11*t^2 + 19/11*t - 40/11 : 6/11*t^3 + 3/11*t^2 + 26/11*t + 310/11 : 1),
             (t : -1/11*t^3 - 6/11*t^2 - 19/11*t - 59/11 : 1),
             (16 : 60 : 1),
             (-3/55*t^3 - 7/55*t^2 - 2/55*t - 133/55 : 6/55*t^3 + 3/55*t^2 + 25/11*t + 156/55 : 1),
             (14/121*t^3 - 15/121*t^2 + 90/121*t + 232/121 : 16/121*t^3 - 69/121*t^2 + 293/121*t - 46/121 : 1),
             (-26/121*t^3 + 20/121*t^2 - 219/121*t - 995/121 : -15/121*t^3 - 156/121*t^2 + 232/121*t - 2887/121 : 1),
             (10/121*t^3 + 49/121*t^2 + 168/121*t + 73/121 : -32/121*t^3 - 60/121*t^2 + 261/121*t + 686/121 : 1),
             (5 : 5 : 1),
             (-9/121*t^3 - 21/121*t^2 - 127/121*t - 377/121 : -7/121*t^3 + 24/121*t^2 + 197/121*t + 16/121 : 1),
             (3/55*t^3 + 7/55*t^2 + 2/55*t + 78/55 : 7/55*t^3 - 24/55*t^2 + 9/11*t + 17/55 : 1),
             (-5/121*t^3 + 36/121*t^2 - 84/121*t + 24/121 : -34/121*t^3 + 27/121*t^2 - 305/121*t - 829/121 : 1),
             (5/121*t^3 - 14/121*t^2 - 158/121*t - 453/121 : 49/121*t^3 + 129/121*t^2 + 315/121*t + 86/121 : 1),
             (5 : -6 : 1),
             (5/121*t^3 - 14/121*t^2 - 158/121*t - 453/121 : -49/121*t^3 - 129/121*t^2 - 315/121*t - 207/121 : 1),
             (-5/121*t^3 + 36/121*t^2 - 84/121*t + 24/121 : 34/121*t^3 - 27/121*t^2 + 305/121*t + 708/121 : 1),
             (3/55*t^3 + 7/55*t^2 + 2/55*t + 78/55 : -7/55*t^3 + 24/55*t^2 - 9/11*t - 72/55 : 1),
             (-9/121*t^3 - 21/121*t^2 - 127/121*t - 377/121 : 7/121*t^3 - 24/121*t^2 - 197/121*t - 137/121 : 1),
             (16 : -61 : 1),
             (10/121*t^3 + 49/121*t^2 + 168/121*t + 73/121 : 32/121*t^3 + 60/121*t^2 - 261/121*t - 807/121 : 1),
             (-26/121*t^3 + 20/121*t^2 - 219/121*t - 995/121 : 15/121*t^3 + 156/121*t^2 - 232/121*t + 2766/121 : 1),
             (14/121*t^3 - 15/121*t^2 + 90/121*t + 232/121 : -16/121*t^3 + 69/121*t^2 - 293/121*t - 75/121 : 1),
             (-3/55*t^3 - 7/55*t^2 - 2/55*t - 133/55 : -6/55*t^3 - 3/55*t^2 - 25/11*t - 211/55 : 1)]

        ::

            sage: E = EllipticCurve('15a1')
            sage: K.<t> = NumberField(x^2 + 2*x + 10)
            sage: EK = E.base_extend(K)
            sage: EK.torsion_points()
            [(0 : 1 : 0),
             (-7 : -5*t - 2 : 1),
             (-7 : 5*t + 8 : 1),
             (-13/4 : 9/8 : 1),
             (-2 : -2 : 1),
             (-2 : 3 : 1),
             (-t - 2 : -t - 7 : 1),
             (-t - 2 : 2*t + 8 : 1),
             (-1 : 0 : 1),
             (t : t - 5 : 1),
             (t : -2*t + 4 : 1),
             (1/2 : -5/4*t - 2 : 1),
             (1/2 : 5/4*t + 1/2 : 1),
             (3 : -2 : 1),
             (8 : -27 : 1),
             (8 : 18 : 1)]

        ::

            sage: K.<i> = QuadraticField(-1)
            sage: EK = EllipticCurve(K, [0,0,0,0,-1])
            sage: EK.torsion_points()
             [(0 : 1 : 0),
              (-2 : -3*i : 1),
              (-2 : 3*i : 1),
              (0 : -i : 1),
              (0 : i : 1),
              (1 : 0 : 1)]
        """
        T = self.torsion_subgroup()  # cached
        return sorted(T.points())           # these are also cached in T

    def rank_bounds(self, **kwds):
        r"""
        Return the lower and upper bounds using :meth:`~simon_two_descent`.
        The results of :meth:`~simon_two_descent` are cached.

        .. NOTE::

            The optional parameters control the Simon two descent algorithm;
            see the documentation of :meth:`~simon_two_descent` for more
            details.

        INPUT:

        - ``verbose`` -- 0, 1, 2, or 3 (default: 0); the verbosity level

        - ``lim1`` -- (default: 2) limit on trivial points on quartics

        - ``lim3`` -- (default: 4) limit on points on ELS quartics

        - ``limtriv`` -- (default: 2) limit on trivial points on elliptic curve

        - ``maxprob`` -- (default: 20)

        - ``limbigprime`` -- (default: 30) to distinguish between
          small and large prime numbers. Use probabilistic tests for
          large primes. If 0, do not use probabilistic tests.

        - ``known_points`` -- (default: ``None``) list of known points on
          the curve

        OUTPUT: lower and upper bounds for the rank of the Mordell-Weil group

        .. NOTE::

            For non-quadratic number fields, this code does
            return, but it takes a long time.

        EXAMPLES::

            sage: x = polygen(ZZ, 'x')
            sage: K.<a> = NumberField(x^2 + 23, 'a')
            sage: E = EllipticCurve(K, '37')
            sage: E == loads(dumps(E))
            True
            sage: E.rank_bounds()
            (2, 2)

        Here is a curve with two-torsion, again the bounds coincide::

            sage: Qrt5.<rt5> = NumberField(x^2 - 5)
            sage: E = EllipticCurve([0, 5-rt5, 0, rt5, 0])
            sage: E.rank_bounds()
            (1, 1)

        Finally an example with non-trivial 2-torsion in Sha. So the
        2-descent will not be able to determine the rank, but can only
        give bounds::

            sage: E = EllipticCurve("15a5")
            sage: K.<t> = NumberField(x^2 - 6)
            sage: EK = E.base_extend(K)
            sage: EK.rank_bounds(lim1=1, lim3=1, limtriv=1)
            (0, 2)

        IMPLEMENTATION:

        Uses Denis Simon's PARI/GP scripts from
        http://www.math.unicaen.fr/~simon/.
        """
        lower, upper, gens = self.simon_two_descent(**kwds)
        # this was corrected in trac 13593. upper is the dimension
        # of the 2-selmer group, so we can certainly remove the
        # 2-torsion of the Mordell-Weil group.
        upper -= self.two_torsion_rank()
        return lower, upper

    def rank(self, **kwds):
        r"""
        Return the rank of this elliptic curve, if it can be determined.

        .. NOTE::

            The optional parameters control the Simon two descent algorithm;
            see the documentation of :meth:`~simon_two_descent` for more
            details.

        INPUT:

        - ``verbose`` -- 0, 1, 2, or 3 (default: 0); the verbosity level

        - ``lim1`` -- (default: 2) limit on trivial points on quartics

        - ``lim3`` -- (default: 4) limit on points on ELS quartics

        - ``limtriv`` -- (default: 2) limit on trivial points on elliptic curve

        - ``maxprob`` -- (default: 20)

        - ``limbigprime`` -- (default: 30) to distinguish between
          small and large prime numbers. Use probabilistic tests for
          large primes. If 0, do not use probabilistic tests.

        - ``known_points`` -- (default: ``None``) list of known points on
          the curve

        OUTPUT:

        If the upper and lower bounds given by Simon two-descent are
        the same, then the rank has been uniquely identified and we
        return this. Otherwise, we raise a :exc:`ValueError` with an error
        message specifying the upper and lower bounds.

        .. NOTE::

            For non-quadratic number fields, this code does return, but it takes
            a long time.

        EXAMPLES::

            sage: x = polygen(ZZ, 'x')
            sage: K.<a> = NumberField(x^2 + 23, 'a')
            sage: E = EllipticCurve(K, '37')
            sage: E == loads(dumps(E))
            True
            sage: E.rank()
            2

        Here is a curve with two-torsion in the Tate-Shafarevich group,
        so here the bounds given by the algorithm do not uniquely
        determine the rank::

            sage: E = EllipticCurve("15a5")
            sage: K.<t> = NumberField(x^2 - 6)
            sage: EK = E.base_extend(K)
            sage: EK.rank(lim1=1, lim3=1, limtriv=1)
            Traceback (most recent call last):
            ...
            ValueError: There is insufficient data to determine the rank -
            2-descent gave lower bound 0 and upper bound 2

        IMPLEMENTATION:

        Uses Denis Simon's PARI/GP scripts from
        http://www.math.unicaen.fr/~simon/.
        """
        lower, upper = self.rank_bounds(**kwds)
        if lower == upper:
            return lower
        else:
            raise ValueError('There is insufficient data to determine the rank - 2-descent gave lower bound %s and upper bound %s' % (lower, upper))

    def gens(self, **kwds):
        r"""
        Return some points of infinite order on this elliptic curve.

        Contrary to what the name of this method suggests, the points
        it returns do not always generate a subgroup of full rank in
        the Mordell-Weil group, nor are they necessarily linearly
        independent.  Moreover, the number of points can be smaller or
        larger than what one could expect after calling :meth:`~rank`
        or :meth:`~rank_bounds`.

        .. NOTE::

            The optional parameters control the Simon two descent algorithm;
            see the documentation of :meth:`~simon_two_descent` for more
            details.

        INPUT:

        - ``verbose`` -- 0, 1, 2, or 3 (default: 0); the verbosity level

        - ``lim1`` -- (default: 2) limit on trivial points on quartics

        - ``lim3`` -- (default: 4) limit on points on ELS quartics

        - ``limtriv`` -- (default: 2) limit on trivial points on elliptic curve

        - ``maxprob`` -- (default: 20)

        - ``limbigprime`` -- (default: 30) to distinguish between
          small and large prime numbers. Use probabilistic tests for
          large primes. If 0, do not use probabilistic tests.

        - ``known_points`` -- (default: ``None``) list of known points on
          the curve

        OUTPUT: a set of points of infinite order given by the Simon two-descent

        .. NOTE::

            For non-quadratic number fields, this code does return, but it takes
            a long time.

        EXAMPLES::

            sage: x = polygen(ZZ, 'x')
            sage: K.<a> = NumberField(x^2 + 23, 'a')
            sage: E = EllipticCurve(K,[0,0,0,101,0])
            sage: E.gens()
            [(23831509/8669448*a - 2867471/8669448 : 76507317707/18049790736*a - 424166479633/18049790736 : 1),
             (-2031032029/969232392*a + 58813561/969232392 : -15575984630401/21336681877488*a + 451041199309/21336681877488 : 1),
             (-186948623/4656964 : 549438861195/10049728312*a : 1)]

        It can happen that no points are found if the height bounds
        used in the search are too small (see :issue:`10745`)::

            sage: K.<t> = NumberField(x^4 + x^2 - 7)
            sage: E = EllipticCurve(K, [1, 0, 5*t^2 + 16, 0, 0])
            sage: E.gens(lim1=1, lim3=1)
            []
            sage: E.rank()
            1
            sage: gg=E.gens(lim3=13); gg  # long time (about 4s)
            [(... : 1)]

        Check that the the point found has infinite order, and that it is on the curve::

            sage: P=gg[0]; P.order()  # long time
            +Infinity
            sage: E.defining_polynomial()(*P)  # long time
            0

        Here is a curve of rank 2::

            sage: K.<t> = NumberField(x^2 - 17)
            sage: E = EllipticCurve(K, [-4, 0])
            sage: E.gens()
            [(-1/2*t + 1/2 : -1/2*t + 1/2 : 1), (-t + 3 : -2*t + 10 : 1)]
            sage: E.rank()
            2

        Test that points of finite order are not included (see :issue:`13593`)::

            sage: E = EllipticCurve("17a3")
            sage: K.<t> = NumberField(x^2 + 3)
            sage: EK = E.base_extend(K)
            sage: EK.rank()
            0
            sage: EK.gens()
            []

        IMPLEMENTATION:

        For curves over quadratic fields which are base-changes from
        `\QQ`, we delegate the work to :meth:`gens_quadratic` where
        methods over `\QQ` suffice.  Otherwise, we use Denis Simon's
        PARI/GP scripts from http://www.math.unicaen.fr/~simon/.
        """
        try:
            return self.gens_quadratic(**kwds)
        except ValueError:
            self.simon_two_descent(**kwds)
            return self._known_points

    def period_lattice(self, embedding):
        r"""
        Return the period lattice of the elliptic curve for the given
        embedding of its base field with respect to the differential
        `dx/(2y + a_1x + a_3)`.

        INPUT:

        - ``embedding`` -- an embedding of the base number field into `\RR` or `\CC`

        .. NOTE::

            The precision of the embedding is ignored: we only use the
            given embedding to determine which embedding into ``QQbar``
            to use.  Once the lattice has been initialized, periods can
            be computed to arbitrary precision.

        EXAMPLES:

        First define a field with two real embeddings::

            sage: K.<a> = NumberField(x^3 - 2)
            sage: E = EllipticCurve([0,0,0,a,2])
            sage: embs = K.embeddings(CC); len(embs)
            3

        For each embedding we have a different period lattice::

            sage: E.period_lattice(embs[0])
            Period lattice associated to Elliptic Curve defined by y^2 = x^3 + a*x + 2
             over Number Field in a with defining polynomial x^3 - 2
             with respect to the embedding Ring morphism:
              From: Number Field in a with defining polynomial x^3 - 2
              To:   Algebraic Field
              Defn: a |--> -0.6299605249474365? - 1.091123635971722?*I

            sage: E.period_lattice(embs[1])
            Period lattice associated to Elliptic Curve defined by y^2 = x^3 + a*x + 2
             over Number Field in a with defining polynomial x^3 - 2
             with respect to the embedding Ring morphism:
              From: Number Field in a with defining polynomial x^3 - 2
              To:   Algebraic Field
              Defn: a |--> -0.6299605249474365? + 1.091123635971722?*I

            sage: E.period_lattice(embs[2])
            Period lattice associated to Elliptic Curve defined by y^2 = x^3 + a*x + 2
             over Number Field in a with defining polynomial x^3 - 2
             with respect to the embedding Ring morphism:
              From: Number Field in a with defining polynomial x^3 - 2
              To:   Algebraic Field
              Defn: a |--> 1.259921049894873?

        Although the original embeddings have only the default
        precision, we can obtain the basis with higher precision
        later::

            sage: L = E.period_lattice(embs[0])
            sage: L.basis()
            (1.86405007647981 - 0.903761485143226*I, -0.149344633143919 - 2.06619546272945*I)

            sage: L.basis(prec=100)
            (1.8640500764798108425920506200 - 0.90376148514322594749786960975*I,
             -0.14934463314391922099120107422 - 2.0661954627294548995621225062*I)
        """
        from sage.schemes.elliptic_curves.period_lattice import PeriodLattice_ell
        return PeriodLattice_ell(self,embedding)

    def real_components(self, embedding):
        """
        Return the number of real components with respect to a real embedding of the base field.

        EXAMPLES::

            sage: K.<a> = QuadraticField(5)
            sage: embs = K.real_embeddings()
            sage: E = EllipticCurve([0,1,1,a,a])
            sage: [e(E.discriminant()) > 0 for e in embs]
            [True, False]
            sage: [E.real_components(e) for e in embs]
            [2, 1]

        TESTS::

            sage: K.<a> = QuadraticField(-1)
            sage: e = K.complex_embeddings()[0]
            sage: E = EllipticCurve([0,1,1,a,a])
            sage: E.real_components(e)
            Traceback (most recent call last):
            ...
            ValueError: invalid embedding specified: should be real

            sage: E.real_components('banana')
            Traceback (most recent call last):
            ...
            ValueError: invalid embedding

            sage: K.<a> = QuadraticField(-1)
            sage: E = EllipticCurve([0,1,1,a,a])
            sage: K1.<a> = QuadraticField(5)
            sage: e = K1.real_embeddings()[0]
            sage: E.real_components(e)
            Traceback (most recent call last):
            ...
            ValueError: invalid embedding specified: should have domain ...
        """
        try:
            if embedding.domain() is not self.base_field():
                raise ValueError("invalid embedding specified: should have domain {}".format(self.base_field()))
            if not isinstance(embedding.codomain(), sage.rings.abc.RealField):
                raise ValueError("invalid embedding specified: should be real")
        except AttributeError:
                raise ValueError("invalid embedding")

        from sage.rings.number_field.number_field import refine_embedding
        from sage.rings.infinity import Infinity
        e = refine_embedding(embedding,Infinity)

        return 2 if e(self.discriminant()) > 0 else 1

    def height_function(self):
        """
        Return the canonical height function attached to ``self``.

        EXAMPLES::

            sage: x = polygen(ZZ, 'x')
            sage: K.<a> = NumberField(x^2 - 5)
            sage: E = EllipticCurve(K, '11a3')
            sage: E.height_function()
            EllipticCurveCanonicalHeight object associated to
             Elliptic Curve defined by y^2 + y = x^3 + (-1)*x^2
             over Number Field in a with defining polynomial x^2 - 5
        """
        if not hasattr(self, '_height_function'):
            from sage.schemes.elliptic_curves.height import EllipticCurveCanonicalHeight
            self._height_function = EllipticCurveCanonicalHeight(self)
        return self._height_function

    ##########################################################
    # Isogeny class
    ##########################################################
    def isogeny_class(self, reducible_primes=None, algorithm='Billerey', minimal_models=True):
        r"""
        Return the isogeny class of this elliptic curve.

        INPUT:

        - ``reducible_primes`` -- list of integers or ``None`` (default)); if
          not ``None`` then this should be a list of primes; in computing
          the isogeny class, only composites isogenies of these
          degrees will be used.

        - ``algorithm`` -- string (default: ``'Billerey'``); the algorithm
          to use to compute the reducible primes.  Ignored for CM
          curves or if ``reducible_primes`` is provided.  Values are
          ``'Billerey'`` (default), ``'Larson'``, and ``'heuristic'``.

        - ``minimal_models`` -- boolean (default: ``True``); if ``True``,
          all curves in the class will be minimal or semi-minimal
          models.  Over fields of larger degree it can be expensive to
          compute these so set to ``False``.

        OUTPUT:

        An instance of the class
        :class:`sage.schemes.elliptic_curves.isogeny_class.IsogenyClass_EC_NumberField`.
        From this object may be obtained a list of curves in the
        class, a matrix of the degrees of the isogenies between them,
        and the isogenies themselves.

        .. NOTE::

            If using the algorithm ``'heuristic'`` for non-CM curves, the
            result is not guaranteed to be the complete isogeny class,
            since only reducible primes up to the default bound in
            :meth:`reducible_primes_naive` (currently 1000) are
            tested.  However, no examples of non-CM elliptic curves
            with reducible primes greater than 100 have yet been
            computed so the output is likely to be correct.

        .. NOTE::

            By default, the curves in the isogeny class will all be
            minimal models if these exist (for example, when the class
            number is `1`); otherwise they will be minimal at all but
            one prime.  This behaviour can be switched off if desired,
            for example over fields where the computation of the class
            group would be too expensive.

        EXAMPLES::

            sage: K.<i> = QuadraticField(-1)
            sage: E = EllipticCurve(K, [0,0,0,0,1])
            sage: C = E.isogeny_class(); C
            Isogeny class of Elliptic Curve defined by y^2 = x^3 + 1
             over Number Field in i with defining polynomial x^2 + 1 with i = 1*I

        The curves in the class (sorted)::

            sage: [E1.ainvs() for E1 in C]
            [(0, 0, 0, 0, -27),
            (0, 0, 0, 0, 1),
            (i + 1, i, 0, 3, -i),
            (i + 1, i, 0, 33, 91*i)]

        The matrix of degrees of cyclic isogenies between curves::

            sage: C.matrix()
            [1 3 6 2]
            [3 1 2 6]
            [6 2 1 3]
            [2 6 3 1]

        The array of isogenies themselves is not filled out but only
        contains those used to construct the class, the other entries
        containing the integer 0.  This will be changed when the
        class :class:`EllipticCurveIsogeny` allowed composition.  In
        this case we used `2`-isogenies to go from 0 to 2 and from 1
        to 3, and `3`-isogenies to go from 0 to 1 and from 2 to 3::

            sage: isogs = C.isogenies()
            sage: [((i,j), isogs[i][j].degree())
            ....:  for i in range(4) for j in range(4) if isogs[i][j] != 0]
            [((0, 1), 3),
             ((0, 3), 2),
             ((1, 0), 3),
             ((1, 2), 2),
             ((2, 1), 2),
             ((2, 3), 3),
             ((3, 0), 2),
             ((3, 2), 3)]
            sage: [((i,j), isogs[i][j].x_rational_map())
            ....:  for i in range(4) for j in range(4) if isogs[i][j] != 0]
            [((0, 1), (1/9*x^3 - 12)/x^2),
             ((0, 3), (1/2*i*x^2 - 2*i*x + 15*i)/(x - 3)),
             ((1, 0), (x^3 + 4)/x^2),
             ((1, 2), (1/2*i*x^2 + i)/(x + 1)),
             ((2, 1), (-1/2*i*x^2 - 1/2*i)/(x - 1/2*i)),
             ((2, 3), (x^3 - 2*i*x^2 - 7*x + 4*i)/(x^2 - 2*i*x - 1)),
             ((3, 0), (-1/2*i*x^2 + 2*x - 5/2*i)/(x + 7/2*i)),
             ((3, 2), (1/9*x^3 + 2/3*i*x^2 - 13/3*x - 116/9*i)/(x^2 + 10*i*x - 25))]

        The isogeny class may be visualized by obtaining its graph and
        plotting it::

            sage: G = C.graph()
            sage: G.show(edge_labels=True) # long time

            sage: K.<i> = QuadraticField(-1)
            sage: E = EllipticCurve([1+i, -i, i, 1, 0])
            sage: C = E.isogeny_class(); C # long time
            Isogeny class of
             Elliptic Curve defined by y^2 + (i+1)*x*y + i*y = x^3 + (-i)*x^2 + x
              over Number Field in i with defining polynomial x^2 + 1 with i = 1*I
            sage: len(C) # long time
            6
            sage: C.matrix() # long time
            [ 1  3  9 18  6  2]
            [ 3  1  3  6  2  6]
            [ 9  3  1  2  6 18]
            [18  6  2  1  3  9]
            [ 6  2  6  3  1  3]
            [ 2  6 18  9  3  1]
            sage: [E1.ainvs() for E1 in C] # long time
            [(i + 1, i - 1, i, -i - 1, -i + 1),
             (i + 1, i - 1, i, 14*i + 4, 7*i + 14),
             (i + 1, i - 1, i, 59*i + 99, 372*i - 410),
             (i + 1, -i, i, -240*i - 399, 2869*i + 2627),
             (i + 1, -i, i, -5*i - 4, 2*i + 5),
             (i + 1, -i, i, 1, 0)]

        An example with CM by `\sqrt{-5}`::

            sage: pol = PolynomialRing(QQ,'x')([1,0,3,0,1])
            sage: K.<c> = NumberField(pol)
            sage: j = 1480640 + 565760*c^2
            sage: E = EllipticCurve(j=j)
            sage: E.has_cm()
            True
            sage: E.has_rational_cm()
            True
            sage: E.cm_discriminant()
            -20
            sage: C = E.isogeny_class()
            sage: len(C)
            2
            sage: C.matrix()
            [1 2]
            [2 1]
            sage: [E.ainvs() for E in C]
            [(0, 0, 0, 83490*c^2 - 147015, -64739840*c^2 - 84465260),
            (0, 0, 0, -161535*c^2 + 70785, -62264180*c^3 + 6229080*c)]
            sage: C.isogenies()[0][1]
            Isogeny of degree 2
             from Elliptic Curve defined by
                  y^2 = x^3 + (83490*c^2-147015)*x + (-64739840*c^2-84465260)
                  over Number Field in c with defining polynomial x^4 + 3*x^2 + 1
               to Elliptic Curve defined by
                  y^2 = x^3 + (-161535*c^2+70785)*x + (-62264180*c^3+6229080*c)
                  over Number Field in c with defining polynomial x^4 + 3*x^2 + 1

        An example with CM by `\sqrt{-23}` (class number `3`)::

            sage: pol = PolynomialRing(QQ,'x')([1,-3,5,-5,5,-3,1])
            sage: L.<a> = NumberField(pol)
            sage: js = hilbert_class_polynomial(-23).roots(L,multiplicities=False); len(js)
            3
            sage: E = EllipticCurve(j=js[0])
            sage: E.has_rational_cm()
            True
            sage: len(E.isogenies_prime_degree())  # long time
            3
            sage: C = E.isogeny_class(); len(C)  # long time
            6

        The reason for the isogeny class having size six while the
        class number is only `3` is that the class also contains three
        curves with CM by the order of discriminant `-92=4\cdot(-23)`,
        which also has class number `3`.  The curves in the class are
        sorted first by CM discriminant (then lexicographically using
        a-invariants)::

            sage: [F.cm_discriminant() for F in C]  # long time
            [-23, -23, -23, -92, -92, -92]

        `2` splits in the order with discriminant `-23`, into two
        primes of order `3` in the class group, each of which induces
        a `2`-isogeny to a curve with the same endomorphism ring; the
        third `2`-isogeny is to a curve with the smaller endomorphism
        ring::

            sage: [phi.codomain().cm_discriminant() for phi in E.isogenies_prime_degree()]  # long time
            [-92, -23, -23]

            sage: C.matrix()  # long time # random
            [1 2 2 4 4 2]
            [2 1 2 4 2 4]
            [2 2 1 2 4 4]
            [4 4 2 1 3 3]
            [4 2 4 3 1 3]
            [2 4 4 3 3 1]

        The graph of this isogeny class has a shape which does not
        occur over `\QQ`: a triangular prism.  Note that for curves
        without CM, the graph has an edge between two curves if and
        only if they are connected by an isogeny of prime degree, and
        this degree is uniquely determined by the two curves, but in
        the CM case this property does not hold, since for pairs of
        curves in the class with the same endomorphism ring `O`, the
        set of degrees of isogenies between them is the set of
        integers represented by a primitive integral binary quadratic
        form of discriminant `\text{disc}(O)`, and this form
        represents infinitely many primes.  In the matrix we give a
        small prime represented by the appropriate form.  In this
        example, the matrix is formed by four `3\times3` blocks.  The
        isogenies of degree `2` indicated by the upper left `3\times3`
        block of the matrix could be replaced by isogenies of any
        degree represented by the quadratic form `2x^2+xy+3y^2` of
        discriminant `-23`.  Similarly in the lower right block, the
        entries of `3` could be represented by any integers
        represented by the quadratic form `3x^2+2xy+8y^2` of
        discriminant `-92`.  In the top right block and lower left
        blocks, by contrast, the prime entries `2` are unique
        determined::

            sage: G = C.graph()  # long time
            sage: G.adjacency_matrix()  # long time # random
            [0 1 1 0 0 1]
            [1 0 1 0 1 0]
            [1 1 0 1 0 0]
            [0 0 1 0 1 1]
            [0 1 0 1 0 1]
            [1 0 0 1 1 0]
            sage: Graph(polytopes.simplex(2).prism().adjacency_matrix()).is_isomorphic(G) # long time
            True

        To display the graph without any edge labels::

            sage: G.show()  # not tested

        To display the graph with edge labels: by default, for curves
        with rational CM, the labels are the coefficients of the
        associated quadratic forms::

            sage: G.show(edge_labels=True)  # not tested

        For an alternative view, first relabel the edges using only 2
        labels to distinguish between isogenies between curves with
        the same endomorphism ring and isogenies between curves with
        different endomorphism rings, then use a 3-dimensional plot
        which can be rotated::

            sage: for i, j, l in G.edge_iterator():  # long time
            ....:     G.set_edge_label(i, j, l.count(','))
            sage: G.show3d(color_by_label=True)  # long time

        A class number `6` example.  First we set up the fields: ``pol``
        defines the same field as ``pol26`` but is simpler::

            sage: pol26 = hilbert_class_polynomial(-4*26)
            sage: pol = x^6 - x^5 + 2*x^4 + x^3 - 2*x^2 - x - 1
            sage: K.<a> = NumberField(pol)
            sage: L.<b> = K.extension(x^2 + 26)

        Only `2` of the `j`-invariants with discriminant -104 are in
        `K`, though all are in `L`::

            sage: len(pol26.roots(K))
            2
            sage: len(pol26.roots(L))
            6

        We create an elliptic curve defined over `K` with one of the
        `j`-invariants in `K`::

            sage: j1 = pol26.roots(K)[0][0]
            sage: E = EllipticCurve(j=j1)
            sage: E.has_cm()
            True
            sage: E.has_rational_cm()
            False
            sage: E.has_rational_cm(L)
            True

        Over `K` the isogeny class has size `4`, with `2` curves for
        each of the `2` `K`-rational `j`-invariants::

            sage: C = E.isogeny_class(); len(C) # long time (~11s)
            4
            sage: C.matrix()                    # long time
            [ 1 13  2 26]
            [13  1 26  2]
            [ 2 26  1 13]
            [26  2 13  1]
            sage: len(Set([EE.j_invariant() for EE in C.curves]))  # long time
            2

        Over `L`, the isogeny class grows to size `6` (the class
        number)::

            sage: EL = E.change_ring(L)
            sage: CL = EL.isogeny_class(minimal_models=False) # long time
            sage: len(CL) # long time
            6
            sage: s1 = Set([EE.j_invariant() for EE in CL.curves]) # long time
            sage: s2 = Set(pol26.roots(L, multiplicities=False)) # long time
            sage: s1 == s2 # long time
            True

        In each position in the matrix of degrees, we see primes (or
        `1`).  In fact the set of degrees of cyclic isogenies from
        curve `i` to curve `j` is infinite, and is the set of all
        integers represented by one of the primitive binary quadratic
        forms of discriminant `-104`, from which we have selected a
        small prime::

            sage: CL.matrix() # long time # random (see :issue:`19229`)
            [1 2 3 3 5 5]
            [2 1 5 5 3 3]
            [3 5 1 3 2 5]
            [3 5 3 1 5 2]
            [5 3 2 5 1 3]
            [5 3 5 2 3 1]

        To see the array of binary quadratic forms::

            sage: CL.qf_matrix()  # long time # random (see :issue:`19229`)
            [[[1], [2, 0, 13], [3, -2, 9], [3, -2, 9], [5, -4, 6], [5, -4, 6]],
             [[2, 0, 13], [1], [5, -4, 6], [5, -4, 6], [3, -2, 9], [3, -2, 9]],
             [[3, -2, 9], [5, -4, 6], [1], [3, -2, 9], [2, 0, 13], [5, -4, 6]],
             [[3, -2, 9], [5, -4, 6], [3, -2, 9], [1], [5, -4, 6], [2, 0, 13]],
             [[5, -4, 6], [3, -2, 9], [2, 0, 13], [5, -4, 6], [1], [3, -2, 9]],
             [[5, -4, 6], [3, -2, 9], [5, -4, 6], [2, 0, 13], [3, -2, 9], [1]]]

        As in the non-CM case, the isogeny class may be visualized by
        obtaining its graph and plotting it.  Since there are more
        edges than in the non-CM case, it may be preferable to omit
        the edge_labels::

            sage: G = C.graph()
            sage: G.show(edge_labels=False) # long time

        It is possible to display a 3-dimensional plot, with colours
        to represent the different edge labels, in a form which can be
        rotated!::

            sage: G.show3d(color_by_label=True) # long time

        Over larger number fields several options make computations
        tractable.  Here we use algorithm 'heuristic' which avoids a
        rigorous computation of the reducible primes, only testing
        those less than 1000, and setting ``minimal_models`` to
        ``False`` avoid having to compute the class group of `K`.  To
        obtain minimal models set proof.number_field(False); the class
        group computation takes an additional 10s::

            sage: K.<z> = CyclotomicField(53)
            sage: E = EllipticCurve(K,[0,6,0,2,0])
            sage: C = E.isogeny_class(algorithm='heuristic', minimal_models=False); C # long time (10s)
            Isogeny class of Elliptic Curve defined by y^2 = x^3 + 6*x^2 + 2*x
             over Cyclotomic Field of order 53 and degree 52
            sage: C.curves # long time
            [Elliptic Curve defined by y^2 = x^3 + 6*x^2 + (-8)*x + (-48)
              over Cyclotomic Field of order 53 and degree 52,
             Elliptic Curve defined by y^2 = x^3 + 6*x^2 + 2*x
              over Cyclotomic Field of order 53 and degree 52]

        TESTS:

        An example which failed until fixed at :issue:`19229`::

            sage: K.<a> = NumberField(x^2 - x + 1)
            sage: E = EllipticCurve([a+1, 1, 1, 0, 0])
            sage: C = E.isogeny_class(); len(C) # long time
            4

        Check that :issue:`36780` is fixed::

            sage: L5.<r5> = NumberField(x^2-5)
            sage: F = EllipticCurve(L5,[0,-4325477943600 *r5-4195572876000])
            sage: F.isogeny_class().matrix()
            [ 1 25 75  3  5 15]
            [25  1  3 75  5 15]
            [75  3  1 25 15  5]
            [ 3 75 25  1 15  5]
            [ 5  5 15 15  1  3]
            [15 15  5  5  3  1]
        """
        try:
            return self._isoclass
        except AttributeError:
            from sage.schemes.elliptic_curves.isogeny_class import IsogenyClass_EC_NumberField
            self._isoclass = IsogenyClass_EC_NumberField(self, reducible_primes=reducible_primes, algorithm=algorithm, minimal_models=minimal_models)
            return self._isoclass

    def isogenies_prime_degree(self, l=None, algorithm='Billerey', minimal_models=True):
        r"""
        Return a list of `\ell`-isogenies from ``self``, where `\ell` is a
        prime.

        INPUT:

        - ``l`` -- either ``None`` or a prime or a list of primes

        - ``algorithm`` -- string (default: ``'Billerey'``); the algorithm
          to use to compute the reducible primes when ``l`` is ``None``.
          Ignored for CM curves or if ``l`` is provided.  Values are
          'Billerey' (default), 'Larson', and 'heuristic'.

        - ``minimal_models`` -- boolean (default: ``True``); if ``True``,
          all curves computed will be minimal or semi-minimal models.
          Over fields of larger degree it can be expensive to compute
          these so set to ``False``.

        OUTPUT: list; `\ell`-isogenies for the given `\ell` or if `\ell` is
        ``None``, all isogenies of prime degree (see below for the CM case)

        .. NOTE::

            Over `\QQ`, the codomains of the isogenies returned are
            standard minimal models.  Over other number fields they are
            global minimal models if these exist, otherwise models
            which are minimal at all but one prime.

        .. NOTE::

            For curves with rational CM, isogenies of primes degree
            exist for infinitely many primes `\ell`, though there are
            only finitely many isogenous curves up to isomorphism.  The
            list returned only includes one isogeny of prime degree for
            each codomain.

        EXAMPLES::

            sage: K.<i> = QuadraticField(-1)
            sage: E = EllipticCurve(K, [0,0,0,0,1])
            sage: isogs = E.isogenies_prime_degree()
            sage: [phi.degree() for phi in isogs]
            [2, 3]

            sage: pol = PolynomialRing(QQ,'x')([1,-3,5,-5,5,-3,1])
            sage: L.<a> = NumberField(pol)
            sage: js = hilbert_class_polynomial(-23).roots(L, multiplicities=False); len(js)
            3
            sage: E = EllipticCurve(j=js[0])
            sage: len(E.isogenies_prime_degree())  # long time
            3

        Set ``minimal_models`` to False to avoid computing minimal
        models of the isogenous curves, since that can be
        time-consuming since it requires computation of the class
        group::

            sage: proof.number_field(False)
            sage: K.<z> = CyclotomicField(53)
            sage: E = EllipticCurve(K, [0,6,0,2,0])
            sage: E.isogenies_prime_degree(2, minimal_models=False)
            [Isogeny of degree 2
              from Elliptic Curve defined by y^2 = x^3 + 6*x^2 + 2*x
                   over Cyclotomic Field of order 53 and degree 52
                to Elliptic Curve defined by y^2 = x^3 + 6*x^2 + (-8)*x + (-48)
                   over Cyclotomic Field of order 53 and degree 52]
            sage: E.isogenies_prime_degree(2, minimal_models=True) # not tested (10s)
            [Isogeny of degree 2
              from Elliptic Curve defined by y^2 = x^3 + 6*x^2 + 2*x
                   over Cyclotomic Field of order 53 and degree 52
                to Elliptic Curve defined by y^2 = x^3 + (-20)*x + (-16)
                   over Cyclotomic Field of order 53 and degree 52]

        TESTS::

            sage: E.isogenies_prime_degree(4)
            Traceback (most recent call last):
            ...
            ValueError: 4 is not prime.
        """
        from .isogeny_small_degree import isogenies_prime_degree

        if l is not None and not isinstance(l, list):
            try:
                l = ZZ(l)
            except TypeError:
                raise ValueError("%s is not a prime integer" % l)
            try:
                if l.is_prime(proof=False):
                    return isogenies_prime_degree(self, l, minimal_models=minimal_models)
                else:
                    raise ValueError("%s is not prime." % l)
            except AttributeError:
                raise ValueError("%s is not prime." % l)

        if l is None:
            from .isogeny_class import possible_isogeny_degrees
            L = possible_isogeny_degrees(self)
            return self.isogenies_prime_degree(L, minimal_models=minimal_models)

        isogs = sum([self.isogenies_prime_degree(p, minimal_models=minimal_models) for p in l],
                    [])

        if self.has_rational_cm():
            # eliminate any endomorphisms and repeated codomains
            isogs = [phi for phi in isogs if not self.is_isomorphic(phi.codomain())]
            codoms = [phi.codomain() for phi in isogs]
            isogs = [phi for i, phi in enumerate(isogs)
                     if not any(E.is_isomorphic(codoms[i])
                                for E in codoms[:i])]
        return isogs

    def is_isogenous(self, other, proof=True, maxnorm=100):
        """
        Return whether or not ``self`` is isogenous to ``other``.

        INPUT:

        - ``other`` -- another elliptic curve

        - ``proof`` -- boolean (default: ``True``); if ``False``, the function
          will return ``True`` whenever the two curves have the same
          conductor and are isogenous modulo `p` for all primes `p` of
          norm up to ``maxnorm``.  If ``True``, the function returns
          False when the previous condition does not hold, and if it
          does hold we compute the complete isogeny class to see if
          the curves are indeed isogenous.

        - ``maxnorm`` -- integer (default: 100); the maximum norm of
          primes `p` for which isogeny modulo `p` will be checked

        OUTPUT:

        boolean; ``True`` if there is an isogeny from curve ``self`` to
        curve ``other``.

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
            sage: F = NumberField(x^2 - 2, 's'); F
            Number Field in s with defining polynomial x^2 - 2
            sage: E1 = EllipticCurve(F, [7,8])
            sage: E2 = EllipticCurve(F, [0,5,0,1,0])
            sage: E3 = EllipticCurve(F, [0,-10,0,21,0])
            sage: E1.is_isogenous(E2)
            False
            sage: E1.is_isogenous(E1)
            True
            sage: E2.is_isogenous(E2)
            True
            sage: E2.is_isogenous(E1)
            False
            sage: E2.is_isogenous(E3)
            True

        ::

            sage: x = polygen(QQ, 'x')
            sage: F = NumberField(x^2 - 2, 's'); F
            Number Field in s with defining polynomial x^2 - 2
            sage: E = EllipticCurve('14a1')
            sage: EE = EllipticCurve('14a2')
            sage: E1 = E.change_ring(F)
            sage: E2 = EE.change_ring(F)
            sage: E1.is_isogenous(E2)
            True

        ::

            sage: x = polygen(QQ, 'x')
            sage: F = NumberField(x^2 - 2, 's'); F
            Number Field in s with defining polynomial x^2 - 2
            sage: k.<a> = NumberField(x^3 + 7)
            sage: E = EllipticCurve(F, [7,8])
            sage: EE = EllipticCurve(k, [2, 2])
            sage: E.is_isogenous(EE)
            Traceback (most recent call last):
            ...
            ValueError: Second argument must be defined over the same number field.

        Some examples from Cremona's 1981 tables::

            sage: K.<i> = QuadraticField(-1)
            sage: E1 = EllipticCurve([i + 1, 0, 1, -240*i - 400, -2869*i - 2627])
            sage: E1.conductor()
            Fractional ideal (4*i + 7)
            sage: E2 = EllipticCurve([1+i,0,1,0,0])
            sage: E2.conductor()
            Fractional ideal (4*i + 7)
            sage: E1.is_isogenous(E2) # long time
            True
            sage: E1.is_isogenous(E2, proof=False) # faster  (~170ms)
            True

        In this case E1 and E2 are in fact 9-isogenous, as may be
        deduced from the following::

            sage: E3 = EllipticCurve([i + 1, 0, 1, -5*i - 5, -2*i - 5])
            sage: E3.is_isogenous(E1)
            True
            sage: E3.is_isogenous(E2)
            True
            sage: E1.isogeny_degree(E2) # long time
            9

        TESTS:

        Check that :issue:`15890` is fixed::

            sage: K.<s> = QuadraticField(229)
            sage: c4 = 2173 - 235*(1 - s)/2
            sage: c6 = -124369 + 15988*(1 - s)/2
            sage: c4c = 2173 - 235*(1 + s)/2
            sage: c6c = -124369 + 15988*(1 + s)/2
            sage: E = EllipticCurve_from_c4c6(c4, c6)
            sage: Ec = EllipticCurve_from_c4c6(c4c, c6c)
            sage: E.is_isogenous(Ec)
            True

        Check that :issue:`17295` is fixed::

            sage: k.<s> = QuadraticField(2)
            sage: K.<b> = k.extension(x^2 - 3)
            sage: E = EllipticCurve(k, [-3*s*(4 + 5*s), 2*s*(2 + 14*s + 11*s^2)])
            sage: Ec = EllipticCurve(k, [3*s*(4 - 5*s), -2*s*(2 - 14*s + 11*s^2)])
            sage: EK = E.base_extend(K)
            sage: EcK = Ec.base_extend(K)
            sage: EK.is_isogenous(EcK)      # long time (about 3.5 s)
            True
        """
        if not isinstance(other, EllipticCurve_generic):
            raise ValueError("Second argument is not an Elliptic Curve.")
        if self.is_isomorphic(other):
            return True
        K = self.base_field()
        if K != other.base_field():
            raise ValueError("Second argument must be defined over the same number field.")

        E1 = self.integral_model()
        E2 = other.integral_model()
        N = E1.conductor()
        if N != E2.conductor():
            return False

        PI = K.primes_of_degree_one_iter()
        while True:
            P = next(PI)
            if P.absolute_norm() > maxnorm:
                break
            if not P.divides(N):
                if E1.reduction(P).cardinality() != E2.reduction(P).cardinality():
                    return False

        if not proof:
            return True

        #  We first try the easiest cases: primes for which X_0(l) has genus 0:

        for l in [2, 3, 5, 7, 13]:
            if any(E2.is_isomorphic(f.codomain())
                   for f in E1.isogenies_prime_degree(l)):
                return True

        #  Next we try the primes for which X_0^+(l) has genus 0 for
        #  which isogeny-finding is faster than in general:

        from .isogeny_small_degree import hyperelliptic_primes
        for l in hyperelliptic_primes:
            if any(E2.is_isomorphic(f.codomain())
                   for f in E1.isogenies_prime_degree(l)):
                return True

        # Next we try looking modulo some more primes:

        while True:
            if P.absolute_norm() > 10 * maxnorm:
                break
            if not P.divides(N):
                OP = K.residue_field(P)
                if E1.change_ring(OP).cardinality() != E2.change_ring(OP).cardinality():
                    return False
            P = next(PI)

        # Finally we compute the full isogeny class of E1 and check if
        # E2 is isomorphic to any curve in the class:

        return any(E2.is_isomorphic(E3) for E3 in E1.isogeny_class().curves)

        raise NotImplementedError("Curves appear to be isogenous (same conductor, isogenous modulo all primes of norm up to %s), but no isogeny has been constructed." % (10*maxnorm))

    def isogeny_degree(self, other):
        """
        Return the minimal degree of an isogeny between ``self`` and
        ``other``, or 0 if no isogeny exists.

        INPUT:

        - ``other`` -- another elliptic curve

        OUTPUT: integer; the degree of an isogeny from ``self`` to ``other``, or 0

        EXAMPLES::

            sage: x = QQ['x'].0
            sage: F = NumberField(x^2 - 2, 's'); F
            Number Field in s with defining polynomial x^2 - 2
            sage: E = EllipticCurve('14a1')
            sage: EE = EllipticCurve('14a2')
            sage: E1 = E.change_ring(F)
            sage: E2 = EE.change_ring(F)
            sage: E1.isogeny_degree(E2)  # long time
            2
            sage: E2.isogeny_degree(E2)
            1
            sage: E5 = EllipticCurve('14a5').change_ring(F)
            sage: E1.isogeny_degree(E5)  # long time
            6

            sage: E = EllipticCurve('11a1')
            sage: [E2.label() for E2 in cremona_curves([11..20]) if E.isogeny_degree(E2)]
            ['11a1', '11a2', '11a3']

            sage: K.<i> = QuadraticField(-1)
            sage: E = EllipticCurve([1+i, -i, i, 1, 0])
            sage: C = E.isogeny_class() # long time
            sage: [E.isogeny_degree(F) for F in C] # long time
            [2, 6, 18, 9, 3, 1]
        """
        # First deal with some easy cases:
        if self.conductor() != other.conductor():
            return ZZ.zero()

        if self.is_isomorphic(other):
            return ZZ.one()

        C = self.isogeny_class()
        i = C.index(self)  # may not be 0 since curves are sorted
        try:
            return C.matrix()[i][C.index(other)]
        except ValueError:
            return ZZ.zero()

    def reducible_primes(self, algorithm='Billerey', max_l=None,
                         num_l=None, verbose=False):
        r"""
        Return a finite set of primes `\ell` for which `E` has a
        K-rational `\ell`-isogeny.

        For curves without CM the list returned is exactly the finite
        set of primes `\ell` for which the mod-`\ell` Galois
        representation is reducible.  For curves with CM this set is
        infinite; we return a (not necessarily minimal) finite list
        of primes `\ell` such that every curve isogenous to this curve
        can be obtained by a finite sequence of isogenies of degree one
        of the primes in the list.

        INPUT:

        - ``algorithm`` -- string; only relevant for non-CM curves
          Either 'Billerey", to use the methods of [Bil2011]_,
          'Larson' to use Larson's implementation using Galois
          representations, or 'heuristic' (see below).

        - ``max_l`` -- integer or ``None``; only relevant for non-CM
          curves and algorithms 'Billerey' and 'heuristic.  Controls
          the maximum prime used in either algorithm.  If ``None``,
          use the default for that algorithm.

        - ``num_l`` -- integer or ``None``; only relevant for non-CM
          curves and algorithm 'Billerey'.  Controls the maximum
          number of primes used in the algorithm.  If ``None``, use
          the default for that algorithm.

        .. NOTE::

            The 'heuristic' algorithm only checks primes up to the
            bound ``max_l``.  This is faster but not guaranteed to be
            complete.  Both the Billerey and Larson algorithms are
            rigorous.

        EXAMPLES::

            sage: x = polygen(ZZ, 'x')
            sage: K = NumberField(x**2 - 29, 'a'); a = K.gen()
            sage: E = EllipticCurve([1, 0, ((5 + a)/2)**2, 0, 0])
            sage: rho = E.galois_representation()
            sage: rho.reducible_primes() # long time
            [3, 5]
            sage: E.reducible_primes() # long time
            [3, 5]
            sage: K = NumberField(x**2 + 1, 'a')
            sage: E = EllipticCurve_from_j(K(1728)) # CM over K
            sage: rho = E.galois_representation()
            sage: rho.reducible_primes() # CM curves always return [0]
            [0]
            sage: E.reducible_primes()
            [2, 5]
            sage: E = EllipticCurve_from_j(K(0)) # CM but NOT over K
            sage: rho = E.galois_representation()
            sage: rho.reducible_primes() # long time
            [2, 3]
            sage: E.reducible_primes()
            [2, 3]
            sage: E = EllipticCurve_from_j(K(2268945/128)).global_minimal_model() # c.f. [Sut2012]
            sage: rho = E.galois_representation()
            sage: rho.isogeny_bound() # ..but there is no 7-isogeny, long time
            [7]
            sage: rho.reducible_primes() # long time
            []
            sage: E.reducible_primes() # long time
            []
        """
        from sage.schemes.elliptic_curves.isogeny_class import possible_isogeny_degrees
        return possible_isogeny_degrees(self, max_l=max_l, num_l=num_l, exact=True, verbose=verbose)

    def lll_reduce(self, points, height_matrix=None, precision=None):
        r"""
        Return an LLL-reduced basis from a given basis, with transform
        matrix.

        INPUT:

        - ``points`` -- list of points on this elliptic
          curve, which should be independent

        - ``height_matrix`` -- the height-pairing matrix of
          the points, or ``None``. If ``None``, it will be computed.

        - ``precision`` -- number of bits of precision of intermediate
          computations (default: ``None``, for default RealField
          precision; ignored if ``height_matrix`` is supplied)

        OUTPUT: a tuple (newpoints, U) where U is a unimodular integer
        matrix, new_points is the transform of points by U, such that
        new_points has LLL-reduced height pairing matrix

        .. NOTE::

            If the input points are not independent, the output
            depends on the undocumented behaviour of PARI's
            :pari:`qflllgram` function when applied to a Gram matrix which
            is not positive definite.

        EXAMPLES:

        Some examples over `\QQ`::

            sage: E = EllipticCurve([0, 1, 1, -2, 42])
            sage: Pi = E.gens(); Pi
            [(-4 : 1 : 1), (-3 : 5 : 1), (-11/4 : 43/8 : 1), (-2 : 6 : 1)]
            sage: Qi, U = E.lll_reduce(Pi)
            sage: all(sum(U[i,j]*Pi[i] for i in range(4)) == Qi[j] for j in range(4))
            True
            sage: sorted(Qi)
            [(-4 : 1 : 1), (-3 : 5 : 1), (-2 : 6 : 1), (0 : 6 : 1)]
            sage: U.det()
            1
            sage: E.regulator_of_points(Pi)
            4.59088036960573
            sage: E.regulator_of_points(Qi)
            4.59088036960574

        ::

            sage: E = EllipticCurve([1,0,1,-120039822036992245303534619191166796374,504224992484910670010801799168082726759443756222911415116])
            sage: xi = [2005024558054813068,
            ....:       -4690836759490453344,
            ....:       4700156326649806635,
            ....:       6785546256295273860,
            ....:       6823803569166584943,
            ....:       7788809602110240789,
            ....:       27385442304350994620556,
            ....:       54284682060285253719/4,
            ....:       -94200235260395075139/25,
            ....:       -3463661055331841724647/576,
            ....:       -6684065934033506970637/676,
            ....:       -956077386192640344198/2209,
            ....:       -27067471797013364392578/2809,
            ....:       -25538866857137199063309/3721,
            ....:       -1026325011760259051894331/108241,
            ....:       9351361230729481250627334/1366561,
            ....:       10100878635879432897339615/1423249,
            ....:       11499655868211022625340735/17522596,
            ....:       110352253665081002517811734/21353641,
            ....:       414280096426033094143668538257/285204544,
            ....:       36101712290699828042930087436/4098432361,
            ....:       45442463408503524215460183165/5424617104,
            ....:       983886013344700707678587482584/141566320009,
            ....:       1124614335716851053281176544216033/152487126016]
            sage: points = [E.lift_x(x) for x in xi]
            sage: newpoints, U = E.lll_reduce(points)  # long time (35s on sage.math, 2011)
            sage: [P[0] for P in newpoints]            # long time
            [6823803569166584943, 5949539878899294213, 2005024558054813068, 5864879778877955778, 23955263915878682727/4, 5922188321411938518, 5286988283823825378, 11465667352242779838, -11451575907286171572, 3502708072571012181, 1500143935183238709184/225, 27180522378120223419/4, -5811874164190604461581/625, 26807786527159569093, 7041412654828066743, 475656155255883588, 265757454726766017891/49, 7272142121019825303, 50628679173833693415/4, 6951643522366348968, 6842515151518070703, 111593750389650846885/16, 2607467890531740394315/9, -1829928525835506297]

        An example to show the explicit use of the height pairing matrix::

            sage: E = EllipticCurve([0, 1, 1, -2, 42])
            sage: Pi = E.gens()
            sage: H = E.height_pairing_matrix(Pi,3)
            sage: E.lll_reduce(Pi,height_matrix=H)
            (
                                                                      [1 0 0 1]
                                                                      [0 1 0 1]
                                                                      [0 0 0 1]
            [(-4 : 1 : 1), (-3 : 5 : 1), (-2 : 6 : 1), (1 : -7 : 1)], [0 0 1 1]
            )

        Some examples over number fields (see :issue:`9411`)::

            sage: K.<a> = QuadraticField(-23, 'a')
            sage: E = EllipticCurve(K, [0,0,1,-1,0])
            sage: P = E(-2, -(a+1)/2)
            sage: Q = E(0,-1)
            sage: E.lll_reduce([P,Q])
            (
                                                    [ 0 -1]
            [(0 : -1 : 1), (-2 : 1/2*a - 1/2 : 1)], [ 1  0]
            )

        ::

            sage: K.<a> = QuadraticField(-5)
            sage: E = EllipticCurve(K, [0,a])
            sage: points = [E.point([-211/841*a - 6044/841,-209584/24389*a + 53634/24389]),
            ....:           E.point([-17/18*a - 1/9, -109/108*a - 277/108])]
            sage: E.lll_reduce(points)
            (
            [(-a + 4 : -3*a + 7 : 1), (-17/18*a - 1/9 : -109/108*a - 277/108 : 1)],
            <BLANKLINE>
            [1 0]
            [1 1]
            )
        """
        r = len(points)
        if height_matrix is None:
            height_matrix = self.height_pairing_matrix(points, precision)
        U = height_matrix.__pari__().lllgram().sage()
        new_points = [sum([U[j, i]*points[j] for j in range(r)])
                      for i in range(r)]
        return new_points, U

    def galois_representation(self):
        r"""
        The compatible family of the Galois representation
        attached to this elliptic curve.

        Given an elliptic curve `E` over a number field `K`
        and a rational prime number `p`, the `p^n`-torsion
        `E[p^n]` points of `E` is a representation of the
        absolute Galois group of `K`. As `n` varies
        we obtain the Tate module `T_p E` which is a
        a representation of `G_K` on a free `\ZZ_p`-module
        of rank `2`. As `p` varies the representations
        are compatible.

        EXAMPLES::

            sage: x = polygen(ZZ, 'x')
            sage: K = NumberField(x**2 + 1, 'a')
            sage: E = EllipticCurve('11a1').change_ring(K)
            sage: rho = E.galois_representation()
            sage: rho
            Compatible family of Galois representations associated to the
             Elliptic Curve defined by y^2 + y = x^3 + (-1)*x^2 + (-10)*x + (-20)
              over Number Field in a with defining polynomial x^2 + 1
            sage: rho.is_surjective(3)
            True
            sage: rho.is_surjective(5)  # long time (4s on sage.math, 2014)
            False
            sage: rho.non_surjective()
            [5]
        """
        from .gal_reps_number_field import GaloisRepresentation
        return GaloisRepresentation(self)

    @cached_method
    def cm_discriminant(self):
        """
        Return the CM discriminant of the `j`-invariant of this curve, or 0.

        OUTPUT:

        An integer `D` which is either `0` if this curve `E` does not
        have Complex Multiplication) (CM), or an imaginary quadratic
        discriminant if `j(E)` is the `j`-invariant of the order with
        discriminant `D`.

        .. NOTE::

            If `E` has CM but the discriminant `D` is not a square in
            the base field `K` then the extra endomorphisms will not be
            defined over `K`.  See also :meth:`has_rational_cm`.

        EXAMPLES::

            sage: EllipticCurve(j=0).cm_discriminant()
            -3
            sage: EllipticCurve(j=1).cm_discriminant()
            Traceback (most recent call last):
            ...
            ValueError: Elliptic Curve defined by y^2 + x*y = x^3 + 36*x + 3455
            over Rational Field does not have CM
            sage: EllipticCurve(j=1728).cm_discriminant()
            -4
            sage: EllipticCurve(j=8000).cm_discriminant()
            -8
            sage: K.<a> = QuadraticField(5)
            sage: EllipticCurve(j=282880*a + 632000).cm_discriminant()
            -20
            sage: x = polygen(ZZ, 'x')
            sage: K.<a> = NumberField(x^3 - 2)
            sage: EllipticCurve(j=31710790944000*a^2 + 39953093016000*a + 50337742902000).cm_discriminant()
            -108
        """
        from sage.schemes.elliptic_curves.cm import is_cm_j_invariant
        flag, df = is_cm_j_invariant(self.j_invariant())
        if flag:
            d, f = df
            return d*f**2
        else:  # no CM
            return ZZ.zero()

    @cached_method
    def has_cm(self):
        """
        Return whether or not this curve has a CM `j`-invariant.

        OUTPUT:

        ``True`` if this curve has CM over the algebraic closure
        of the base field, otherwise ``False``.  See also
        :meth:`cm_discriminant` and :meth:`has_rational_cm`.

        .. NOTE::

            Even if `E` has CM in this sense (that its `j`-invariant is
            a CM `j`-invariant), if the associated negative
            discriminant `D` is not a square in the base field `K`, the
            extra endomorphisms will not be defined over `K`.  See also
            the method :meth:`has_rational_cm` which tests whether `E`
            has extra endomorphisms defined over `K` or a given
            extension of `K`.

        EXAMPLES::

            sage: EllipticCurve(j=0).has_cm()
            True
            sage: EllipticCurve(j=1).has_cm()
            False
            sage: EllipticCurve(j=1728).has_cm()
            True
            sage: EllipticCurve(j=8000).has_cm()
            True
            sage: K.<a> = QuadraticField(5)
            sage: EllipticCurve(j=282880*a + 632000).has_cm()
            True
            sage: x = polygen(ZZ, 'x')
            sage: K.<a> = NumberField(x^3 - 2)
            sage: EllipticCurve(j=31710790944000*a^2 + 39953093016000*a + 50337742902000).has_cm()
            True
        """
        return not self.cm_discriminant().is_zero()

    @cached_method
    def has_rational_cm(self, field=None):
        r"""
        Return whether or not this curve has CM defined over its
        base field or a given extension.

        INPUT:

        - ``field`` -- a field, which should be an extension of the
          base field of the curve.  If ``field`` is ``None`` (the
          default), it is taken to be the base field of the curve.

        OUTPUT:

        ``True`` if the ring of endomorphisms of this curve over
        the given field is larger than `\ZZ`; otherwise ``False``.
        See also :meth:`cm_discriminant` and :meth:`has_cm`.

        .. NOTE::

            If `E` has CM but the discriminant `D` is not a square in
            the given field `K` then the extra endomorphisms will not
            be defined over `K`, and this function will return
            ``False``.  See also :meth:`has_cm`.  To obtain the CM
            discriminant, use :meth:`cm_discriminant`.

        EXAMPLES::

            sage: E = EllipticCurve(j=0)
            sage: E.has_cm()
            True
            sage: E.has_rational_cm()
            False
            sage: D = E.cm_discriminant(); D
            -3
            sage: E.has_rational_cm(QuadraticField(D))
            True

            sage: E = EllipticCurve(j=1728)
            sage: E.has_cm()
            True
            sage: E.has_rational_cm()
            False
            sage: D = E.cm_discriminant(); D
            -4
            sage: E.has_rational_cm(QuadraticField(D))
            True

        Higher degree examples::

            sage: K.<a> = QuadraticField(5)
            sage: E = EllipticCurve(j=282880*a + 632000)
            sage: E.has_cm()
            True
            sage: E.has_rational_cm()
            False
            sage: E.cm_discriminant()
            -20
            sage: E.has_rational_cm(K.extension(x^2 + 5, 'b'))
            True

        An error is raised if a field is given which is not an extension of the base field::

            sage: E.has_rational_cm(QuadraticField(-20))
            Traceback (most recent call last):
            ...
            ValueError: Error in has_rational_cm: Number Field in a
            with defining polynomial x^2 + 20 with a = 4.472135954999579?*I
            is not an extension field of Number Field in a
            with defining polynomial x^2 - 5 with a = 2.236067977499790?

            sage: x = polygen(ZZ, 'x')
            sage: K.<a> = NumberField(x^3 - 2)
            sage: E = EllipticCurve(j=31710790944000*a^2 + 39953093016000*a + 50337742902000)
            sage: E.has_cm()
            True
            sage: E.has_rational_cm()
            False
            sage: D = E.cm_discriminant(); D
            -108
            sage: E.has_rational_cm(K.extension(x^2 + 108,'b'))
            True
        """
        D = self.cm_discriminant()
        if D.is_zero():
            return False
        if field is None:
            return self.base_field()(D).is_square()
        if self.base_field().embeddings(field):
            D = field(D)
            return D.is_square()
        raise ValueError("Error in has_rational_cm: %s is not an extension field of %s"
                         % (field,self.base_field()))

    @cached_method
    def is_Q_curve(self, maxp=100, certificate=False, verbose=False):
        r"""
        Return ``True`` if this is a `\QQ`-curve, with optional certificate.

        INPUT:

        - ``maxp`` -- integer (default: 100); bound on primes used for
          checking necessary local conditions.  The result will not
          depend on this, but using a larger value may return
          ``False`` faster.

        - ``certificate`` -- boolean (default: ``False``); if ``True`` then
          a second value is returned giving a certificate for the
          `\QQ`-curve property

        OUTPUT:

        If ``certificate`` is ``False``: either ``True`` (if `E` is a
        `\QQ`-curve), or ``False``.

        If ``certificate`` is ``True``: a tuple consisting of a boolean
        flag as before and a certificate, defined as follows:

        - when the flag is ``True``, so `E` is a `\QQ`-curve:

          - either {'CM':`D`} where `D` is a negative discriminant, when
            `E` has potential CM with discriminant `D`;

          - otherwise {'CM': `0`, 'core_poly': `f`, 'rho': `\rho`,
            'r': `r`, 'N': `N`}, when `E` is a non-CM `\QQ`-curve,
            where the core polynomial `f` is an irreducible monic
            polynomial over `QQ` of degree `2^\rho`, all of whose
            roots are `j`-invariants of curves isogenous to `E`, the
            core level `N` is a square-free integer with `r` prime
            factors which is the LCM of the degrees of the isogenies
            between these conjugates.  For example, if there exists
            a curve `E'` isogenous to `E` with `j(E')=j\in\QQ`, then
            the certificate is {'CM':0, 'r':0, 'rho':0, 'core_poly':
            x-j, 'N':1}.

        - when the flag is ``False``, so `E` is not a `\QQ`-curve, the
          certificate is a prime `p` such that the reductions of `E`
          at the primes dividing `p` are inconsistent with the
          property of being a `\QQ`-curve.  See the documentation for
          :meth:`sage.src.schemes.elliptic_curves.Qcurves.is_Q_curve`
          for details.

        ALGORITHM:

        See the documentation for
        :meth:`sage.src.schemes.elliptic_curves.Qcurves.is_Q_curve`, and
        [CrNa2020]_ for details.

        EXAMPLES:

        A non-CM curve over `\QQ` and a CM curve over `\QQ` are both
        trivially `\QQ`-curves::

            sage: E = EllipticCurve([1,2,3,4,5])
            sage: flag, cert = E.is_Q_curve(certificate=True)
            sage: flag
            True
            sage: cert
            {'CM': 0, 'N': 1, 'core_poly': x, 'r': 0, 'rho': 0}

            sage: E = EllipticCurve(j=8000)
            sage: flag, cert = E.is_Q_curve(certificate=True)
            sage: flag
            True
            sage: cert
            {'CM': -8}

        A non-`\QQ`-curve over a quartic field.  The local data at bad
        primes above `3` is inconsistent::

            sage: R.<x> = PolynomialRing(QQ)
            sage: K.<a> = NumberField(R([3, 0, -5, 0, 1]))
            sage: E = EllipticCurve([K([-3,-4,1,1]), K([4,-1,-1,0]), K([-2,0,1,0]), K([-621,778,138,-178]), K([9509,2046,-24728,10380])])
            sage: E.is_Q_curve(certificate=True, verbose=True)
            Checking whether Elliptic Curve defined by y^2 + (a^3+a^2-4*a-3)*x*y + (a^2-2)*y = x^3 + (-a^2-a+4)*x^2 + (-178*a^3+138*a^2+778*a-621)*x + (10380*a^3-24728*a^2+2046*a+9509) over Number Field in a with defining polynomial x^4 - 5*x^2 + 3 is a Q-curve
            No: inconsistency at the 2 primes dividing 3
            - potentially multiplicative: [True, False]
            (False, 3)

        A non-`\QQ`-curve over a quadratic field.  The local data at bad
        primes is consistent, but the local test at good primes above `13`
        is not::

            sage: K.<a> = NumberField(R([-10, 0, 1]))
            sage: E = EllipticCurve([K([0,1]), K([-1,-1]), K([0,0]), K([-236,40]), K([-1840,464])])
            sage: E.is_Q_curve(certificate=True, verbose=True)
            Checking whether Elliptic Curve defined by y^2 + a*x*y = x^3 + (-a-1)*x^2 + (40*a-236)*x + (464*a-1840) over Number Field in a with defining polynomial x^2 - 10 is a Q-curve
            Applying local tests at good primes above p<=100
            No: inconsistency at the 2 ordinary primes dividing 13
            - Frobenius discriminants mod squares: [-1, -3]
            No: local test at p=13 failed
            (False, 13)

        A quadratic `\QQ`-curve with CM discriminant `-15`
        (so the `j`-invariant is not in `\QQ`)::

            sage: R.<x> = PolynomialRing(QQ)
            sage: K.<a> = NumberField(R([-1, -1, 1]))
            sage: E = EllipticCurve([K([1,0]), K([-1,0]), K([0,1]), K([0,-2]), K([0,1])])
            sage: E.is_Q_curve(certificate=True, verbose=True)
            Checking whether Elliptic Curve defined by y^2 + x*y + a*y = x^3 + (-1)*x^2 + (-2*a)*x + a over Number Field in a with defining polynomial x^2 - x - 1 is a Q-curve
            Yes: E is CM (discriminant -15)
            (True, {'CM': -15})

        An example over `\QQ(\sqrt{2},\sqrt{3})`.  The `j`-invariant
        is in `\QQ(\sqrt{6})`, so computations will be done over that
        field, and in fact there is an isogenous curve with rational
        `j`, so we have a so-called rational `\QQ`-curve::

            sage: K.<a> = NumberField(R([1, 0, -4, 0, 1]))
            sage: E = EllipticCurve([K([-2,-4,1,1]), K([0,1,0,0]), K([0,1,0,0]), K([-4780,9170,1265,-2463]), K([163923,-316598,-43876,84852])])
            sage: flag, cert = E.is_Q_curve(certificate=True) # long time
            sage: flag # long time
            True
            sage: cert # long time
            {'CM': 0, 'N': 1, 'core_degs': [1], 'core_poly': x - 85184/3, 'r': 0, 'rho': 0}

        Over the same field, a so-called strict `\QQ`-curve which is
        not isogenous to one with rational `j`, but whose core field
        is quadratic. In fact the isogeny class over `K` consists of
        `6` curves, four with conjugate quartic `j`-invariants and `2`
        with quadratic conjugate `j`-invariants in `\QQ(\sqrt{3})`
        (but which are not base-changes from the quadratic subfield)::

            sage: E = EllipticCurve([K([0,-3,0,1]), K([1,4,0,-1]), K([0,0,0,0]), K([-2,-16,0,4]), K([-19,-32,4,8])])
            sage: flag, cert = E.is_Q_curve(certificate=True) # long time
            sage: flag # long time
            True
            sage: cert # long time
            {'CM': 0,
             'N': 2,
             'core_degs': [1, 2],
             'core_poly': x^2 - 840064*x + 1593413632,
             'r': 1,
             'rho': 1}
        """
        from sage.schemes.elliptic_curves.Qcurves import is_Q_curve as isQ
        return isQ(self, maxp, certificate, verbose)

    def saturation(self, points, verbose=False,
                   max_prime=0, one_prime=0, odd_primes_only=False,
                   lower_ht_bound=None, reg=None, debug=False):
        r"""
        Given a list of rational points on `E` over `K`, compute the
        saturation in `E(K)` of the subgroup they generate.

        INPUT:

        - ``points (list)`` -- list of points on E.  Points of finite
          order are ignored; the remaining points should be independent,
          or an error is raised.

        - ``verbose`` -- boolean (default: ``False``); if ``True``, give
          verbose output

        - ``max_prime`` -- integer (default: 0); saturation is performed
          for all primes up to ``max_prime``. If ``max_prime`` is 0,
          perform saturation at *all* primes, i.e., compute the true
          saturation.

        - ``odd_primes_only`` -- boolean (default: ``False``); only do
          saturation at odd primes

        - ``one_prime`` -- integer (default: 0); if nonzero, only do
          saturation at this prime

        The following two inputs are optional, and may be provided to speed
        up the computation.

        - ``lower_ht_bound`` -- real (default: ``None``); lower bound of
          the regulator `E(K)`, if known

        - ``reg`` -- real (default: ``None``); regulator of the span of
          points, if known

        - ``debug`` -- integer (default: 0); used for debugging and
          testing

        OUTPUT:

        - ``saturation`` -- list; points that form a basis for the
          saturation

        - ``index`` -- integer; the index of the group generated by the
          input points in their saturation

        - ``regulator`` -- real with default precision, or ``None``;
          regulator of saturated points

        EXAMPLES::

            sage: K.<i> = QuadraticField(-1)
            sage: E = EllipticCurve('389a1')
            sage: EK = E.change_ring(K)
            sage: P = EK(-1,1); Q = EK(0,-1)

            sage: EK.saturation([2*P], max_prime=2)
            ([(-1 : 1 : 1)], 2, 0.686667083305587)
            sage: EK.saturation([12*P], max_prime=2)
            ([(26/361 : -5720/6859 : 1)], 4, 6.18000374975028)
            sage: EK.saturation([12*P], lower_ht_bound=0.1)
            ([(-1 : 1 : 1)], 12, 0.686667083305587)
            sage: EK.saturation([2*P, Q], max_prime=2)
            ([(-1 : 1 : 1), (0 : -1 : 1)], 2, 0.152460177943144)
            sage: EK.saturation([P + Q, P - Q], lower_ht_bound=.1, debug=2)
            ([(-1 : 1 : 1), (1 : 0 : 1)], 2, 0.152460177943144)
            sage: EK.saturation([P + Q, 17*Q], lower_ht_bound=0.1)  # long time
            ([(4 : 8 : 1), (0 : -1 : 1)], 17, 0.152460177943143)

            sage: R = EK(i-2,-i-3)
            sage: EK.saturation([P + R, P + Q, Q + R], lower_ht_bound=0.1)
            ([(841/1369*i - 171/1369 : 41334/50653*i - 74525/50653 : 1),
              (4 : 8 : 1),
              (-1/25*i + 18/25 : -69/125*i - 58/125 : 1)],
             2,
             0.103174443217351)
            sage: EK.saturation([26*Q], lower_ht_bound=0.1)
            ([(0 : -1 : 1)], 26, 0.327000773651605)

        Another number field::

            sage: E = EllipticCurve('389a1')
            sage: K.<a> = NumberField(x^3 - x + 1)
            sage: EK = E.change_ring(K)
            sage: P = EK(-1,1); Q = EK(0,-1)
            sage: EK.saturation([P + Q, P - Q], lower_ht_bound=0.1)
            ([(-1 : 1 : 1), (1 : 0 : 1)], 2, 0.152460177943144)
            sage: EK.saturation([3*P, P + 5*Q], lower_ht_bound=0.1)
            ([(-185/2209 : -119510/103823 : 1), (80041/34225 : -26714961/6331625 : 1)],
             15,
             0.152460177943144)

        A different curve::

            sage: K.<a> = QuadraticField(3)
            sage: E = EllipticCurve('37a1')
            sage: EK = E.change_ring(K)
            sage: P = EK(0,0); Q = EK(2-a, 2*a-4)
            sage: EK.saturation([3*P - Q, 3*P + Q], lower_ht_bound=.01)
            ([(0 : 0 : 1), (1/2*a : -1/4*a - 1/4 : 1)], 6, 0.0317814530725985)

        The points must be linearly independent::

            sage: EK.saturation([2*P, 3*Q, P-Q])
            Traceback (most recent call last):
            ...
            ValueError: points not linearly independent in saturation()

        Degenerate case::

            sage: EK.saturation([])
            ([], 1, 1.00000000000000)

        ALGORITHM:

        For rank 1 subgroups, simply do trial division up to the maximal
        prime divisor. For higher rank subgroups, perform trial division
        on all linear combinations for small primes, and look for
        projections `E(K) \rightarrow \oplus E(k) \otimes \mathbf{F}_p` which
        are either full rank or provide `p`-divisible linear combinations,
        where the `k` here are residue fields of `K`.

        TESTS::

            sage: K.<i> = QuadraticField(-1)
            sage: E = EllipticCurve('389a1')
            sage: EK = E.change_ring(K)
            sage: P = EK(-1,1); Q = EK(0,-1)

            sage: EK.saturation([P+Q, P-Q], lower_ht_bound=.1, debug=2)
            ([(-1 : 1 : 1), (1 : 0 : 1)], 2, 0.152460177943144)
            sage: EK.saturation([5*P + 6*Q, 5*P - 3*Q], lower_ht_bound=.1)
            ([(-3/4 : -15/8 : 1), (159965/16129 : -67536260/2048383 : 1)],
             45,
             0.152460177943144)
            sage: EK.saturation([5*P + 6*Q, 5*P - 3*Q], lower_ht_bound=.1, debug=2)
            ([(-3/4 : -15/8 : 1), (159965/16129 : -67536260/2048383 : 1)],
             45,
             0.152460177943144)

        See :issue:`27387`::

            sage: K.<a> = NumberField(x^2 - x - 26)
            sage: E = EllipticCurve([a, 1-a, 0, 93-16*a, 3150-560*a])
            sage: P = E([65-35*a/3, (959*a-5377)/9])
            sage: E.saturation([P], one_prime=2)
            ([(-1/4*a + 3/4 : 59/8*a - 317/8 : 1)], 2, 0.344624259712631)
        """
        full_saturation = (max_prime == 0) and (one_prime == 0)
        Plist = (self(P) for P in points)
        Plist = [P for P in Plist if P.has_infinite_order()]
        n = len(Plist)
        index = ZZ.one()

        if n == 0:
            return Plist, index, RealField()(1)

        # compute the list of primes p at which p-saturation is
        # required.

        heights = [P.height() for P in Plist]
        if reg is None:
            reg = self.regulator_of_points(Plist)
        if reg / min(heights) < 1e-6:
            raise ValueError("points not linearly independent in saturation()")
        sat_reg = reg

        from sage.rings.fast_arith import prime_range
        if full_saturation:
            if lower_ht_bound is None:
                # TODO (robertwb): verify this for rank > 1
                if verbose:
                    print("Computing lower height bound..")
                lower_ht_bound = self.height_function().min(.1, 5) ** n
                if verbose:
                    print("..done: %s" % lower_ht_bound)
            index_bound = (reg/lower_ht_bound).sqrt()
            if index_bound < 2:
                if verbose:
                    print("Saturation index bound < 2, points are saturated already.")
                return Plist, index, RealField()(1)
            else:
                if verbose:
                    print("p-saturating for primes p < {}".format(index_bound.ceil()))
            prime_list = prime_range(index_bound.ceil())
        else:
            if one_prime:
                prime_list = [one_prime]
            else:
                prime_list = prime_range(max_prime+1)
        if odd_primes_only and 2 in prime_list:
            prime_list.remove(2)

        # Now saturate at each prime in prime_list.  The dict
        # lin_combs keeps the values of linear combinations of the
        # points, indexed by coefficient tuples, for efficiency; it is
        # reset whenever the point list changes.

        from sage.schemes.elliptic_curves.saturation import EllipticCurveSaturator
        saturator = EllipticCurveSaturator(self, verbose)
        for p in prime_list:
            if full_saturation and (p > index_bound):
                break
            if verbose:
                print("Saturating at p=%s" % p)
            newPlist, expo = saturator.full_p_saturation(Plist, p)
            if expo:
                if verbose:
                    print(" --gaining index %s^%s" % (p, expo))
                pe = p**expo
                index *= pe
                if full_saturation:
                    index_bound /= pe
                sat_reg /= pe**2
                Plist = newPlist
            else:
                if verbose:
                    print(" --already %s-saturated" % p)

        return Plist, index, sat_reg

    def gens_quadratic(self, **kwds):
        r"""
        Return generators for the Mordell-Weil group modulo torsion, for a
        curve which is a base change from `\QQ` to a quadratic field.

        EXAMPLES::

            sage: E = EllipticCurve([1,2,3,40,50])
            sage: E.conductor()
            2123582
            sage: E.gens()
            [(5 : 17 : 1)]
            sage: K.<i> = QuadraticField(-1)
            sage: EK = E.change_ring(K)
            sage: EK.gens_quadratic()
            [(5 : 17 : 1), (-13 : 48*i + 5 : 1)]

            sage: E.change_ring(QuadraticField(3, 'a')).gens_quadratic()
            [(5 : 17 : 1), (-1 : 2*a - 1 : 1), (11/4 : 33/4*a - 23/8 : 1)]

            sage: K.<a> = QuadraticField(-7)
            sage: E = EllipticCurve([0,0,0,197,0])
            sage: E.conductor()
            2483776
            sage: E.gens()
            [(47995604297578081/7389879786648100 : -25038161802544048018837479/635266655830129794121000 : 1)]
            sage: K.<a> = QuadraticField(7)
            sage: E.change_ring(K).gens_quadratic()
            [(-1209642055/59583566*a + 1639995844/29791783 : -377240626321899/1720892553212*a + 138577803462855/245841793316 : 1),
             (1/28 : 393/392*a : 1),
             (-61*a + 162 : 1098*a - 2916 : 1)]

            sage: E = EllipticCurve([1, a])
            sage: E.gens_quadratic()
            Traceback (most recent call last):
            ...
            ValueError: gens_quadratic() requires the elliptic curve to be a base change from Q
        """
        if not kwds:
            try:
                return list(self.__gens)
            except AttributeError:
                pass

        K = self.base_ring()
        if K.absolute_degree() != 2:
            raise ValueError("gens_quadratic() requires the base field to be quadratic")

        EE = self.descend_to(QQ)
        if not EE:
            raise ValueError("gens_quadratic() requires the elliptic curve to be a base change from Q")

        # In all cases there are exactly two distinct curves /Q whose
        # base-change to K is the original.  NB These need not be
        # quadratic twists of each other!  For example, '32a1' and
        # '32a2' are not quadratic twists of each other (each is its
        # own twist by -1) but they become isomorphic over
        # Q(sqrt(-1)).

        EQ1 = EE[0]
        EQ2 = EE[1]
        iso1 = EQ1.change_ring(K).isomorphism_to(self)
        iso2 = EQ2.change_ring(K).isomorphism_to(self)
        gens1 = [iso1(P) for P in EQ1.gens(**kwds)]
        gens2 = [iso2(P) for P in EQ2.gens(**kwds)]
        gens = self.saturation(gens1 + gens2, max_prime=2)[0]
        self.__gens = gens
        return gens

    def rational_points(self, **kwds):
        r"""
        Find rational points on the elliptic curve, all arguments are passed
        on to :meth:`sage.schemes.generic.algebraic_scheme.rational_points`.

        EXAMPLES::

            sage: E = EllipticCurve('37a')
            sage: E.rational_points(bound=8) # long time
            [(-1 : -1 : 1),
             (-1 : 0 : 1),
             (0 : -1 : 1),
             (0 : 0 : 1),
             (0 : 1 : 0),
             (1/4 : -5/8 : 1),
             (1/4 : -3/8 : 1),
             (1 : -1 : 1),
             (1 : 0 : 1),
             (2 : -3 : 1),
             (2 : 2 : 1)]

        Check that :issue:`26677` is fixed::

            sage: E = EllipticCurve("11a1")
            sage: E.rational_points(bound=5)
            [(0 : 1 : 0), (5 : 5 : 1)]
            sage: E.rational_points(bound=6) # long time
            [(0 : 1 : 0), (5 : -6 : 1), (5 : 5 : 1)]

        An example over a number field::

            sage: E = EllipticCurve([1,0])
            sage: pts = E.rational_points(bound=2, F=QuadraticField(-1))
            sage: pts
            [(-a : 0 : 1), (0 : 0 : 1), (0 : 1 : 0), (a : 0 : 1)]
            sage: pts[0] + pts[1]
            (a : 0 : 1)
        """
        from sage.schemes.curves.constructor import Curve
        # we change E to be a plain curve to allow the generic rational
        # points code to reduce mod any prime, whereas an EllipticCurve
        # can only be base changed to good primes.
        E = self
        if 'F' in kwds:
            E = E.change_ring(kwds['F'])

        return [E(pt) for pt in Curve(self).rational_points(**kwds)]
