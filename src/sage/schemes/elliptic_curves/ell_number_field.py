# -*- coding: utf-8 -*-
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

EXAMPLE::

    sage: K.<i> = NumberField(x^2+1)
    sage: E = EllipticCurve([0,4+i])
    sage: E.discriminant()
    -3456*i - 6480
    sage: P= E([i,2])
    sage: P+P
    (-2*i + 9/16 : -9/4*i - 101/64 : 1)

::

    sage: E.has_good_reduction(2+i)
    True
    sage: E.local_data(4+i)
    Local data at Fractional ideal (i + 4):
    Reduction type: bad additive
    Local minimal model: Elliptic Curve defined by y^2 = x^3 + (i+4) over Number Field in i with defining polynomial x^2 + 1
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
    [Isogeny of degree 3 from Elliptic Curve defined by y^2 = x^3 + (i+4) over Number Field in i with defining polynomial x^2 + 1 to Elliptic Curve defined by y^2 = x^3 + (-27*i-108) over Number Field in i with defining polynomial x^2 + 1]

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

#*****************************************************************************
#       Copyright (C) 2007 Robert Bradshaw <robertwb@math.washington.edu>
#                          William Stein   <wstein@gmail.com>
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
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from ell_field import EllipticCurve_field
import ell_point
import sage.matrix.all as matrix
from sage.rings.ring import Ring
from sage.rings.arith import gcd, prime_divisors
from sage.misc.misc import prod
import sage.databases.cremona
import ell_torsion
from ell_generic import is_EllipticCurve

from sage.libs.pari.all import pari, PariError
from gp_simon import simon_two_descent
from constructor import EllipticCurve
from sage.rings.all import PolynomialRing, ZZ, RealField
from sage.misc.misc import verbose, forall
from sage.rings.integer import Integer
from sage.rings.infinity import Infinity # just for verbose output
from sage.rings.arith import valuation

import gal_reps_number_field

class EllipticCurve_number_field(EllipticCurve_field):
    r"""
    Elliptic curve over a number field.

    EXAMPLES::

        sage: K.<i>=NumberField(x^2+1)
        sage: EllipticCurve([i, i - 1, i + 1, 24*i + 15, 14*i + 35])
        Elliptic Curve defined by y^2 + i*x*y + (i+1)*y = x^3 + (i-1)*x^2 + (24*i+15)*x + (14*i+35) over Number Field in i with defining polynomial x^2 + 1
    """
    def __init__(self, x, y=None):
        r"""
        Allow some ways to create an elliptic curve over a number
        field in addition to the generic ones.

        INPUT:

        - ``x``, ``y`` -- see examples.

        EXAMPLES:

        A curve from the database of curves over `\QQ`, but over a larger field:

            sage: K.<i>=NumberField(x^2+1)
            sage: EllipticCurve(K,'389a1')
            Elliptic Curve defined by y^2 + y = x^3 + x^2 + (-2)*x over Number Field in i with defining polynomial x^2 + 1

        Making the field of definition explicitly larger::

            sage: EllipticCurve(K,[0,-1,1,0,0])
            Elliptic Curve defined by y^2 + y = x^3 + (-1)*x^2 over Number Field in i with defining polynomial x^2 + 1

        """
        if y is None:
            if isinstance(x, list):
                ainvs = x
                field = ainvs[0].parent()
        else:
            if isinstance(y, str):
                field = x
                X = sage.databases.cremona.CremonaDatabase()[y]
                ainvs = list(X.a_invariants())
            else:
                field = x
                ainvs = y
        if not (isinstance(field, Ring) and isinstance(ainvs,list)):
            raise TypeError

        EllipticCurve_field.__init__(self, [field(x) for x in ainvs])
        self._point = ell_point.EllipticCurvePoint_number_field

    def simon_two_descent(self, verbose=0, lim1=5, lim3=50, limtriv=10, maxprob=20, limbigprime=30):
        r"""
        Computes lower and upper bounds on the rank of the Mordell-Weil group,
        and a list of independent points. Used internally by the :meth:`~rank`,
        :meth:`~rank_bounds` and :meth:`~gens` methods.

        INPUT:

        - ``verbose`` -- 0, 1, 2, or 3 (default: 0), the verbosity level

        - ``lim1``    -- (default: 5) limit on trivial points on quartics

        - ``lim3``  -- (default: 50) limit on points on ELS quartics

        - ``limtriv`` -- (default: 10) limit on trivial points on elliptic curve

        - ``maxprob`` -- (default: 20)

        - ``limbigprime`` -- (default: 30) to distinguish between
          small and large prime numbers. Use probabilistic tests for
          large primes. If 0, don't use probabilistic tests.

        OUTPUT:

        ``(lower, upper, list)`` where ``lower`` is a lower bound on
        the rank, ``upper`` is an upper bound (the 2-Selmer rank) and
        ``list`` is a list of independent points on the Weierstrass
        model.  The length of ``list`` is equal to either ``lower``,
        or ``lower-1``, since when ``lower`` is less than ``upper``
        and of different parity, the value of ``lower`` is increased by
        1.

        .. note::

           For non-quadratic number fields, this code does return, but
           it takes a long time.

        ALGORITHM:

        Uses Denis Simon's PARI/GP scripts from
        http://www.math.unicaen.fr/~simon/.

        EXAMPLES::

            sage: K.<a> = NumberField(x^2 + 23, 'a')
            sage: E = EllipticCurve(K, '37')
            sage: E == loads(dumps(E))
            True
            sage: E.simon_two_descent()
            (2, 2, [(-1 : 0 : 1), (1/2*a - 5/2 : -1/2*a - 13/2 : 1)])
            sage: E.simon_two_descent(lim1=3, lim3=20, limtriv=5, maxprob=7, limbigprime=10)
            (2, 2, [(-1 : 0 : 1), (-1/8*a + 5/8 : -3/16*a - 9/16 : 1)])

        ::

            sage: K.<a> = NumberField(x^2 + 7, 'a')
            sage: E = EllipticCurve(K, [0,0,0,1,a]); E
            Elliptic Curve defined by y^2  = x^3 + x + a over Number Field in a with defining polynomial x^2 + 7

            sage: v = E.simon_two_descent(verbose=1); v
            courbe elliptique : Y^2 = x^3 + Mod(1, y^2 + 7)*x + Mod(y, y^2 + 7)
            points triviaux sur la courbe = [[1, 1, 0], [Mod(1/2*y + 3/2, y^2 + 7), Mod(-y - 2, y^2 + 7), 1]]
            #S(E/K)[2]    = 2
            #E(K)/2E(K)   = 2
            #III(E/K)[2]  = 1
            rang(E/K)     = 1
            listpointsmwr = [[Mod(1/2*y + 3/2, y^2 + 7), Mod(-y - 2, y^2 + 7), 1]]
            (1, 1, [(1/2*a + 3/2 : -a - 2 : 1)])

            sage: v = E.simon_two_descent(verbose=2)  # random output
            K = bnfinit(y^2 + 7);
            a = Mod(y,K.pol);
            bnfellrank(K, [0,0,0,1,a]);
            courbe elliptique : Y^2 = x^3 + Mod(1, y^2 + 7)*x + Mod(y, y^2 + 7)
            A = 0
            B = Mod(1, y^2 + 7)
            C = Mod(y, y^2 + 7)
            LS2gen = [Mod(Mod(-5, y^2 + 7)*x^2 + Mod(-3*y, y^2 + 7)*x + Mod(8, y^2 + 7), x^3 + Mod(1, y^2 + 7)*x + Mod(y, y^2 + 7)), Mod(Mod(1, y^2 + 7)*x^2 + Mod(1/2*y - 1/2, y^2 + 7)*x - 1, x^3 + Mod(1, y^2 + 7)*x + Mod(y, y^2 + 7))]
            #LS2gen = 2
            Recherche de points triviaux sur la courbe
            points triviaux sur la courbe = [[1, 1, 0], [Mod(1/2*y + 3/2, y^2 + 7), Mod(-y - 2, y^2 + 7), 1]]
            zc = Mod(Mod(-5, y^2 + 7)*x^2 + Mod(-3*y, y^2 + 7)*x + Mod(8, y^2 + 7), x^3 + Mod(1, y^2 + 7)*x + Mod(y, y^2 + 7))
            symbole de Hilbert (Mod(2, y^2 + 7),Mod(-5, y^2 + 7)) = -1
            zc = Mod(Mod(1, y^2 + 7)*x^2 + Mod(1/2*y - 1/2, y^2 + 7)*x + Mod(-1, y^2 + 7), x^3 + Mod(1, y^2 + 7)*x + Mod(y, y^2 + 7))
            symbole de Hilbert (Mod(-2*y + 2, y^2 + 7),Mod(1, y^2 + 7)) = 0
            sol de Legendre = [1, 0, 1]~
            zc*z1^2 = Mod(Mod(2*y - 2, y^2 + 7)*x + Mod(2*y + 10, y^2 + 7), x^3 + Mod(1, y^2 + 7)*x + Mod(y, y^2 + 7))
            quartique : (-1/2*y + 1/2)*Y^2 = x^4 + (-3*y - 15)*x^2 + (-8*y - 16)*x + (-11/2*y - 15/2)
            reduite: Y^2 = (-1/2*y + 1/2)*x^4 - 4*x^3 + (-3*y + 3)*x^2 + (2*y - 2)*x + (1/2*y + 3/2)
            non ELS en [2, [0, 1]~, 1, 1, [1, 1]~]
            zc = Mod(Mod(1, y^2 + 7)*x^2 + Mod(1/2*y + 1/2, y^2 + 7)*x + Mod(-1, y^2 + 7), x^3 + Mod(1, y^2 + 7)*x + Mod(y, y^2 + 7))
            vient du point trivial [Mod(1/2*y + 3/2, y^2 + 7), Mod(-y - 2, y^2 + 7), 1]
            m1 = 1
            m2 = 1
            #S(E/K)[2]    = 2
            #E(K)/2E(K)   = 2
            #III(E/K)[2]  = 1
            rang(E/K)     = 1
            listpointsmwr = [[Mod(1/2*y + 3/2, y^2 + 7), Mod(-y - 2, y^2 + 7), 1]]
            v =  [1, 1, [[Mod(1/2*y + 3/2, y^2 + 7), Mod(-y - 2, y^2 + 7)]]]
            sage: v
            (1, 1, [(1/2*a + 3/2 : -a - 2 : 1)])


        A curve with 2-torsion::

            sage: K.<a> = NumberField(x^2 + 7)
            sage: E = EllipticCurve(K, '15a')
            sage: E.simon_two_descent()  # long time (3s on sage.math, 2013), points can vary
            (1, 3, [...])

        A failure in the PARI/GP script ell.gp (VERSION 25/03/2009) is reported::

            sage: K = CyclotomicField(43).subfields(3)[0][0]
            sage: E = EllipticCurve(K, '37')
            sage: E.simon_two_descent()  # long time (4s on sage.math, 2013)
            Traceback (most recent call last):
            ...
            RuntimeError:
              ***   at top-level: ans=bnfellrank(K,[0,0,1,
              ***                     ^--------------------
              ***   in function bnfellrank: ...eqtheta,rnfeq,bbnf];rang=
              ***   bnfell2descent_gen(b
              ***   ^--------------------
              ***   in function bnfell2descent_gen: ...riv,r=nfsqrt(nf,norm(zc))
              ***   [1];if(DEBUGLEVEL_el
              ***   ^--------------------
              ***   array index (1) out of allowed range [none].
            An error occurred while running Simon's 2-descent program

        """

        try:
            result = self._simon_two_descent_data[lim1,lim3,limtriv,maxprob,limbigprime]
            if verbose == 0:
                return result
        except AttributeError:
            self._simon_two_descent_data = {}
        except KeyError:
            pass

        t = simon_two_descent(self,
                              verbose=verbose, lim1=lim1, lim3=lim3, limtriv=limtriv,
                              maxprob=maxprob, limbigprime=limbigprime)
        prob_rank = Integer(t[0])
        two_selmer_rank = Integer(t[1])
        prob_gens = [self(P) for P in t[2]]
        self._simon_two_descent_data[lim1,lim3,limtriv,maxprob,limbigprime] = (prob_rank, two_selmer_rank, prob_gens)
        return prob_rank, two_selmer_rank, prob_gens

    def height_pairing_matrix(self, points=None, precision=None):
        r"""
        Returns the height pairing matrix of the given points.

        INPUT:

        - points - either a list of points, which must be on this
          curve, or (default) None, in which case self.gens() will be
          used.

        - precision - number of bits of precision of result
          (default: None, for default RealField precision)

        EXAMPLES::

            sage: E = EllipticCurve([0, 0, 1, -1, 0])
            sage: E.height_pairing_matrix()
            [0.0511114082399688]

        For rank 0 curves, the result is a valid 0x0 matrix::

            sage: EllipticCurve('11a').height_pairing_matrix()
            []
            sage: E=EllipticCurve('5077a1')
            sage: E.height_pairing_matrix([E.lift_x(x) for x in [-2,-7/4,1]], precision=100)
            [  1.3685725053539301120518194471  -1.3095767070865761992624519454 -0.63486715783715592064475542573]
            [ -1.3095767070865761992624519454   2.7173593928122930896610589220   1.0998184305667292139777571432]
            [-0.63486715783715592064475542573   1.0998184305667292139777571432  0.66820516565192793503314205089]

            sage: E = EllipticCurve('389a1')
            sage: E = EllipticCurve('389a1')
            sage: P,Q = E.point([-1,1,1]),E.point([0,-1,1])
            sage: E.height_pairing_matrix([P,Q])
            [0.686667083305587 0.268478098806726]
            [0.268478098806726 0.327000773651605]

        Over a number field::

            sage: x = polygen(QQ)
            sage: K.<t> = NumberField(x^2+47)
            sage: EK = E.base_extend(K)
            sage: EK.height_pairing_matrix([EK(P),EK(Q)])
            [0.686667083305587 0.268478098806726]
            [0.268478098806726 0.327000773651605]

        ::

            sage: K.<i> = QuadraticField(-1)
            sage: E = EllipticCurve([0,0,0,i,i])
            sage: P = E(-9+4*i,-18-25*i)
            sage: Q = E(i,-i)
            sage: E.height_pairing_matrix([P,Q])
            [  2.16941934493768 -0.870059380421505]
            [-0.870059380421505  0.424585837470709]
            sage: E.regulator_of_points([P,Q])
            0.164101403936070
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
        M = matrix.MatrixSpace(RR, r)
        mat = M()
        for j in range(r):
            mat[j,j] = points[j].height(precision=precision)
        for j in range(r):
            for k in range(j+1,r):
                mat[j,k]=((points[j]+points[k]).height(precision=precision) - mat[j,j] - mat[k,k])/2
                mat[k,j]=mat[j,k]
        return mat

    def regulator_of_points(self, points=[], precision=None):
        """
        Returns the regulator of the given points on this curve.

        INPUT:

        - ``points`` -(default: empty list)  a list of points on this curve

        - ``precision`` - int or None (default: None): the precision
          in bits of the result (default real precision if None)

        EXAMPLES::

            sage: E = EllipticCurve('37a1')
            sage: P = E(0,0)
            sage: Q = E(1,0)
            sage: E.regulator_of_points([P,Q])
            0.000000000000000
            sage: 2*P==Q
            True

        ::

            sage: E = EllipticCurve('5077a1')
            sage: points = [E.lift_x(x) for x in [-2,-7/4,1]]
            sage: E.regulator_of_points(points)
            0.417143558758384
            sage: E.regulator_of_points(points,precision=100)
            0.41714355875838396981711954462

        ::

            sage: E = EllipticCurve('389a')
            sage: E.regulator_of_points()
            1.00000000000000
            sage: points = [P,Q] = [E(-1,1),E(0,-1)]
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
            sage: E = EllipticCurve(K,[1,1])
            sage: P = E(0,1)
            sage: P.height()
            0.476223106404866
            sage: E.regulator_of_points([P])
            0.476223106404866

        ::

            sage: E = EllipticCurve('11a1')
            sage: x = polygen(QQ)
            sage: K.<t> = NumberField(x^2+47)
            sage: EK = E.base_extend(K)
            sage: T = EK(5,5)
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
            sage: K.<t> = NumberField(x^2+47)
            sage: EK = E.base_extend(K)
            sage: EK.regulator_of_points([EK(P),EK(Q)])
            0.152460177943144

        ::

            sage: K.<i> = QuadraticField(-1)
            sage: E = EllipticCurve([0,0,0,i,i])
            sage: P = E(-9+4*i,-18-25*i)
            sage: Q = E(i,-i)
            sage: E.height_pairing_matrix([P,Q])
            [  2.16941934493768 -0.870059380421505]
            [-0.870059380421505  0.424585837470709]
            sage: E.regulator_of_points([P,Q])
            0.164101403936070

        """
        if points is None:
            points = []
        mat = self.height_pairing_matrix(points=points, precision=precision)
        return mat.det(algorithm="hessenberg")


    def is_local_integral_model(self,*P):
        r"""
        Tests if self is integral at the prime ideal `P`, or at all the
        primes if `P` is a list or tuple.

        INPUT:

        - ``*P`` -- a prime ideal, or a list or tuple of primes.

        EXAMPLES::

            sage: K.<i> = NumberField(x^2+1)
            sage: P1,P2 = K.primes_above(5)
            sage: E = EllipticCurve([i/5,i/5,i/5,i/5,i/5])
            sage: E.is_local_integral_model(P1,P2)
            False
            sage: Emin = E.local_integral_model(P1,P2)
            sage: Emin.is_local_integral_model(P1,P2)
            True
        """
        if len(P)==1: P=P[0]
        if isinstance(P,(tuple,list)):
            return forall(P, lambda x : self.is_local_integral_model(x))[0]
        return forall(self.ainvs(), lambda x : x.valuation(P) >= 0)[0]

    def local_integral_model(self,*P):
        r"""
        Return a model of self which is integral at the prime ideal
        `P`.

        .. note::

           The integrality at other primes is not affected, even if
           `P` is non-principal.

        INPUT:

        - ``*P`` -- a prime ideal, or a list or tuple of primes.

        EXAMPLES::

            sage: K.<i> = NumberField(x^2+1)
            sage: P1,P2 = K.primes_above(5)
            sage: E = EllipticCurve([i/5,i/5,i/5,i/5,i/5])
            sage: E.local_integral_model((P1,P2))
            Elliptic Curve defined by y^2 + (-i)*x*y + (-25*i)*y = x^3 + 5*i*x^2 + 125*i*x + 3125*i over Number Field in i with defining polynomial x^2 + 1
        """
        if len(P)==1: P=P[0]
        if isinstance(P,(tuple,list)):
            E=self
            for Pi in P: E=E.local_integral_model(Pi)
            return E
        ai = self.a_invariants()
        e  = min([(ai[i].valuation(P)/[1,2,3,4,6][i]) for i in range(5)]).floor()
        pi = self.base_field().uniformizer(P, 'negative')
        return EllipticCurve([ai[i]/pi**(e*[1,2,3,4,6][i]) for i in range(5)])

    def is_global_integral_model(self):
        r"""
        Return true iff self is integral at all primes.

        EXAMPLES::

            sage: K.<i> = NumberField(x^2+1)
            sage: E = EllipticCurve([i/5,i/5,i/5,i/5,i/5])
            sage: P1,P2 = K.primes_above(5)
            sage: Emin = E.global_integral_model()
            sage: Emin.is_global_integral_model()
            True
        """
        return forall(self.a_invariants(), lambda x : x.is_integral())[0]

    def global_integral_model(self):
        r"""
        Return a model of self which is integral at all primes.

        EXAMPLES::

            sage: K.<i> = NumberField(x^2+1)
            sage: E = EllipticCurve([i/5,i/5,i/5,i/5,i/5])
            sage: P1,P2 = K.primes_above(5)
            sage: E.global_integral_model()
            Elliptic Curve defined by y^2 + (-i)*x*y + (-25*i)*y = x^3 + 5*i*x^2 + 125*i*x + 3125*i over Number Field in i with defining polynomial x^2 + 1

        :trac:`7935`::

            sage: K.<a> = NumberField(x^2-38)
            sage: E = EllipticCurve([a,1/2])
            sage: E.global_integral_model()
            Elliptic Curve defined by y^2 = x^3 + 1444*a*x + 27436 over Number Field in a with defining polynomial x^2 - 38

        :trac:`9266`::

            sage: K.<s> = NumberField(x^2-5)
            sage: w = (1+s)/2
            sage: E = EllipticCurve(K,[2,w])
            sage: E.global_integral_model()
            Elliptic Curve defined by y^2 = x^3 + 2*x + (1/2*s+1/2) over Number Field in s with defining polynomial x^2 - 5

        :trac:`12151`::

            sage: K.<v> = NumberField(x^2 + 161*x - 150)
            sage: E = EllipticCurve([25105/216*v - 3839/36, 634768555/7776*v - 98002625/1296, 634768555/7776*v - 98002625/1296, 0, 0])
            sage: E.global_integral_model()
            Elliptic Curve defined by y^2 + (-502639783*v+465618899)*x*y + (-6603604211463489399460860*v+6117229527723443603191500)*y = x^3 + (1526887622075335620*v-1414427901517840500)*x^2 over Number Field in v with defining polynomial x^2 + 161*x - 150

        :trac:`14476`::

            sage: R.<t> = QQ[]
            sage: K.<g> = NumberField(t^4 - t^3 - 3*t^2 - t + 1)
            sage: E = EllipticCurve([ -43/625*g^3 + 14/625*g^2 - 4/625*g + 706/625, -4862/78125*g^3 - 4074/78125*g^2 - 711/78125*g + 10304/78125,  -4862/78125*g^3 - 4074/78125*g^2 - 711/78125*g + 10304/78125, 0,0])
            sage: E.global_integral_model()
            Elliptic Curve defined by y^2 + (15*g^3-48*g-42)*x*y + (-111510*g^3-162162*g^2-44145*g+37638)*y = x^3 + (-954*g^3-1134*g^2+81*g+576)*x^2 over Number Field in g with defining polynomial t^4 - t^3 - 3*t^2 - t + 1

        """
        K = self.base_field()
        ai = self.a_invariants()
        Ps = [[ ff[0] for ff in a.denominator_ideal().factor() ] for a in ai if not a.is_integral() ]
        Ps = sage.misc.misc.union(sage.misc.flatten.flatten(Ps))
        for P in Ps:
            pi = K.uniformizer(P,'positive')
            e  = min([(ai[i].valuation(P)/[1,2,3,4,6][i]) for i in range(5)]).floor()
            if e < 0 :
                ai = [ai[i]/pi**(e*[1,2,3,4,6][i]) for i in range(5)]
            if all(a.is_integral() for a in ai):
                break
        for z in ai:
            assert z.is_integral(), "bug in global_integral_model: %s" % list(ai)
        return EllipticCurve(list(ai))

    integral_model = global_integral_model

    def _reduce_model(self):
        r"""

        Transforms the elliptic curve to a model in which `a_1`,
        `a_2`, `a_3` are reduced modulo 2, 3, 2 respectively.

        .. note::

           This only works on integral models, i.e. it requires that
           `a_1`, `a_2` and `a_3` lie in the ring of integers of the base
           field.

        EXAMPLES::

            sage: K.<a>=NumberField(x^2-38)
            sage: E=EllipticCurve([a, -5*a + 19, -39*a + 237, 368258520200522046806318224*a - 2270097978636731786720858047, 8456608930180227786550494643437985949781*a - 52130038506835491453281450568107193773505])
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
            (a1, a2, a3, a4, a6) = [ZK(a) for a in self.a_invariants()]
        except TypeError:
            import sys
            raise TypeError, "_reduce_model() requires an integral model.", sys.exc_info()[2]

        # N.B. Must define s, r, t in the right order.
        if ZK.absolute_degree() == 1:
            s = ((-a1)/2).round('up')
            r = ((-a2 + s*a1 +s*s)/3).round()
            t = ((-a3 - r*a1)/2).round('up')
        else:
            pariK = K._pari_()
            s = K(pariK.nfeltdiveuc(-a1, 2))
            r = K(pariK.nfeltdiveuc(-a2 + s*a1 + s*s, 3))
            t = K(pariK.nfeltdiveuc(-a3 - r*a1, 2))

        return self.rst_transform(r, s, t)

    def local_data(self, P=None, proof = None, algorithm="pari"):
        r"""
        Local data for this elliptic curve at the prime `P`.

        INPUT:

        - ``P`` -- either None or a prime ideal of the base field of self.

        - ``proof`` -- whether to only use provably correct methods
          (default controlled by global proof module).  Note that the
          proof module is number_field, not elliptic_curves, since the
          functions that actually need the flag are in number fields.

        - ``algorithm`` (string, default: "pari") -- Ignored unless the
          base field is `\QQ`.  If "pari", use the PARI C-library
          ``ellglobalred`` implementation of Tate's algorithm over
          `\QQ`. If "generic", use the general number field
          implementation.

        OUTPUT:

        If `P` is specified, returns the ``EllipticCurveLocalData``
        object associated to the prime `P` for this curve.  Otherwise,
        returns a list of such objects, one for each prime `P` in the
        support of the discriminant of this model.

        .. note::

           The model is not required to be integral on input.

           For principal `P`, a generator is used as a uniformizer,
           and integrality or minimality at other primes is not
           affected.  For non-principal `P`, the minimal model
           returned will preserve integrality at other primes, but not
           minimality.

        EXAMPLES::

            sage: K.<i> = NumberField(x^2+1)
            sage: E = EllipticCurve([1 + i, 0, 1, 0, 0])
            sage: E.local_data()
            [Local data at Fractional ideal (i - 2):
            Reduction type: bad non-split multiplicative
            Local minimal model: Elliptic Curve defined by y^2 + (i+1)*x*y + y = x^3 over Number Field in i with defining polynomial x^2 + 1
            Minimal discriminant valuation: 1
            Conductor exponent: 1
            Kodaira Symbol: I1
            Tamagawa Number: 1,
            Local data at Fractional ideal (-3*i - 2):
            Reduction type: bad split multiplicative
            Local minimal model: Elliptic Curve defined by y^2 + (i+1)*x*y + y = x^3 over Number Field in i with defining polynomial x^2 + 1
            Minimal discriminant valuation: 2
            Conductor exponent: 1
            Kodaira Symbol: I2
            Tamagawa Number: 2]
            sage: E.local_data(K.ideal(3))
            Local data at Fractional ideal (3):
            Reduction type: good
            Local minimal model: Elliptic Curve defined by y^2 + (i+1)*x*y + y = x^3 over Number Field in i with defining polynomial x^2 + 1
            Minimal discriminant valuation: 0
            Conductor exponent: 0
            Kodaira Symbol: I0
            Tamagawa Number: 1

        An example raised in :trac:`3897`::

            sage: E = EllipticCurve([1,1])
            sage: E.local_data(3)
            Local data at Principal ideal (3) of Integer Ring:
            Reduction type: good
            Local minimal model: Elliptic Curve defined by y^2 = x^3 + x + 1 over Rational Field
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

        return self._get_local_data(P,proof,algorithm)

    def _get_local_data(self, P, proof, algorithm="pari"):
        r"""
        Internal function to create data for this elliptic curve at the prime `P`.

        This function handles the caching of local data.  It is called
        by local_data() which is the user interface and which parses
        the input parameters `P` and proof.

        INPUT:

        - ``P`` -- either None or a prime ideal of the base field of self.

        - ``proof`` -- whether to only use provably correct methods
          (default controlled by global proof module).  Note that the
          proof module is number_field, not elliptic_curves, since the
          functions that actually need the flag are in number fields.

        - ``algorithm`` (string, default: "pari") -- Ignored unless the
          base field is `\QQ`.  If "pari", use the PARI C-library
          ``ellglobalred`` implementation of Tate's algorithm over
          `\QQ`. If "generic", use the general number field
          implementation.

        EXAMPLES::

            sage: K.<i> = NumberField(x^2+1)
            sage: E = EllipticCurve(K,[0,1,0,-160,308])
            sage: p = K.ideal(i+1)
            sage: E._get_local_data(p, False)
            Local data at Fractional ideal (i + 1):
            Reduction type: good
            Local minimal model: Elliptic Curve defined by y^2 + x*y + y = x^3 + x^2 + (-10)*x + (-10) over Number Field in i with defining polynomial x^2 + 1
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
            return self._local_data[P, proof, algorithm]
        except AttributeError:
            self._local_data = {}
        except KeyError:
            pass
        from sage.schemes.elliptic_curves.ell_local_data import EllipticCurveLocalData
        self._local_data[P, proof, algorithm] = EllipticCurveLocalData(self, P, proof, algorithm)
        return self._local_data[P, proof, algorithm]

    def local_minimal_model(self, P, proof = None, algorithm="pari"):
        r"""
        Returns a model which is integral at all primes and minimal at `P`.

        INPUT:

        - ``P`` -- either None or a prime ideal of the base field of self.

        - ``proof`` -- whether to only use provably correct methods
          (default controlled by global proof module).  Note that the
          proof module is number_field, not elliptic_curves, since the
          functions that actually need the flag are in number fields.

        - ``algorithm`` (string, default: "pari") -- Ignored unless the
          base field is `\QQ`.  If "pari", use the PARI C-library
          ``ellglobalred`` implementation of Tate's algorithm over
          `\QQ`. If "generic", use the general number field
          implementation.

        OUTPUT:

        A model of the curve which is minimal (and integral) at `P`.

        .. note::

           The model is not required to be integral on input.

           For principal `P`, a generator is used as a uniformizer,
           and integrality or minimality at other primes is not
           affected.  For non-principal `P`, the minimal model
           returned will preserve integrality at other primes, but not
           minimality.

        EXAMPLES::

            sage: K.<a>=NumberField(x^2-5)
            sage: E=EllipticCurve([20, 225, 750, 1250*a + 6250, 62500*a + 15625])
            sage: P=K.ideal(a)
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
        Return True if this elliptic curve has good reduction at the prime `P`.

        INPUT:

        - ``P`` -- a prime ideal of the base field of self, or a field
          element generating such an ideal.

        OUTPUT:

        (bool) -- True if the curve has good reduction at `P`, else False.

        .. note::

           This requires determining a local integral minimal model;
           we do not just check that the discriminant of the current
           model has valuation zero.

        EXAMPLES::

            sage: E=EllipticCurve('14a1')
            sage: [(p,E.has_good_reduction(p)) for p in prime_range(15)]
            [(2, False), (3, True), (5, True), (7, False), (11, True), (13, True)]

            sage: K.<a>=NumberField(x^3-2)
            sage: P17a, P17b = [P for P,e in K.factor(17)]
            sage: E = EllipticCurve([0,0,0,0,2*a+1])
            sage: [(p,E.has_good_reduction(p)) for p in [P17a,P17b]]
            [(Fractional ideal (4*a^2 - 2*a + 1), True),
            (Fractional ideal (2*a + 1), False)]
        """
        return self.local_data(P).has_good_reduction()

    def has_bad_reduction(self, P):
        r"""
        Return True if this elliptic curve has bad reduction at the prime `P`.

        INPUT:

        - ``P`` -- a prime ideal of the base field of self, or a field
          element generating such an ideal.

        OUTPUT:

        (bool) True if the curve has bad reduction at `P`, else False.

        .. note::

           This requires determining a local integral minimal model;
           we do not just check that the discriminant of the current
           model has valuation zero.

        EXAMPLES::

            sage: E=EllipticCurve('14a1')
            sage: [(p,E.has_bad_reduction(p)) for p in prime_range(15)]
            [(2, True), (3, False), (5, False), (7, True), (11, False), (13, False)]

            sage: K.<a>=NumberField(x^3-2)
            sage: P17a, P17b = [P for P,e in K.factor(17)]
            sage: E = EllipticCurve([0,0,0,0,2*a+1])
            sage: [(p,E.has_bad_reduction(p)) for p in [P17a,P17b]]
            [(Fractional ideal (4*a^2 - 2*a + 1), False),
            (Fractional ideal (2*a + 1), True)]
        """
        return self.local_data(P).has_bad_reduction()

    def has_multiplicative_reduction(self, P):
        r"""
        Return True if this elliptic curve has (bad) multiplicative reduction at the prime `P`.

        .. note::

           See also ``has_split_multiplicative_reduction()`` and
           ``has_nonsplit_multiplicative_reduction()``.

        INPUT:

        - ``P`` -- a prime ideal of the base field of self, or a field
                 element generating such an ideal.

        OUTPUT:

        (bool) True if the curve has multiplicative reduction at `P`,
        else False.

        EXAMPLES::

            sage: E=EllipticCurve('14a1')
            sage: [(p,E.has_multiplicative_reduction(p)) for p in prime_range(15)]
            [(2, True), (3, False), (5, False), (7, True), (11, False), (13, False)]

            sage: K.<a>=NumberField(x^3-2)
            sage: P17a, P17b = [P for P,e in K.factor(17)]
            sage: E = EllipticCurve([0,0,0,0,2*a+1])
            sage: [(p,E.has_multiplicative_reduction(p)) for p in [P17a,P17b]]
            [(Fractional ideal (4*a^2 - 2*a + 1), False), (Fractional ideal (2*a + 1), False)]
        """
        return self.local_data(P).has_multiplicative_reduction()

    def has_split_multiplicative_reduction(self, P):
        r"""
        Return True if this elliptic curve has (bad) split multiplicative reduction at the prime `P`.

        INPUT:

        - ``P`` -- a prime ideal of the base field of self, or a field
          element generating such an ideal.

        OUTPUT:

        (bool) True if the curve has split multiplicative reduction at
        `P`, else False.

        EXAMPLES::

            sage: E=EllipticCurve('14a1')
            sage: [(p,E.has_split_multiplicative_reduction(p)) for p in prime_range(15)]
            [(2, False), (3, False), (5, False), (7, True), (11, False), (13, False)]

            sage: K.<a>=NumberField(x^3-2)
            sage: P17a, P17b = [P for P,e in K.factor(17)]
            sage: E = EllipticCurve([0,0,0,0,2*a+1])
            sage: [(p,E.has_split_multiplicative_reduction(p)) for p in [P17a,P17b]]
            [(Fractional ideal (4*a^2 - 2*a + 1), False), (Fractional ideal (2*a + 1), False)]
        """
        return self.local_data(P).has_split_multiplicative_reduction()

    def has_nonsplit_multiplicative_reduction(self, P):
        r"""
        Return True if this elliptic curve has (bad) non-split multiplicative reduction at the prime `P`.

        INPUT:

        - ``P`` -- a prime ideal of the base field of self, or a field
          element generating such an ideal.

        OUTPUT:

        (bool) True if the curve has non-split multiplicative
        reduction at `P`, else False.

        EXAMPLES::

            sage: E=EllipticCurve('14a1')
            sage: [(p,E.has_nonsplit_multiplicative_reduction(p)) for p in prime_range(15)]
            [(2, True), (3, False), (5, False), (7, False), (11, False), (13, False)]

            sage: K.<a>=NumberField(x^3-2)
            sage: P17a, P17b = [P for P,e in K.factor(17)]
            sage: E = EllipticCurve([0,0,0,0,2*a+1])
            sage: [(p,E.has_nonsplit_multiplicative_reduction(p)) for p in [P17a,P17b]]
            [(Fractional ideal (4*a^2 - 2*a + 1), False), (Fractional ideal (2*a + 1), False)]
        """
        return self.local_data(P).has_nonsplit_multiplicative_reduction()

    def has_additive_reduction(self, P):
        r"""
        Return True if this elliptic curve has (bad) additive reduction at the prime `P`.

        INPUT:

        - ``P`` -- a prime ideal of the base field of self, or a field
          element generating such an ideal.

        OUTPUT:

        (bool) True if the curve has additive reduction at `P`, else False.

        EXAMPLES::

            sage: E=EllipticCurve('27a1')
            sage: [(p,E.has_additive_reduction(p)) for p in prime_range(15)]
            [(2, False), (3, True), (5, False), (7, False), (11, False), (13, False)]

            sage: K.<a>=NumberField(x^3-2)
            sage: P17a, P17b = [P for P,e in K.factor(17)]
            sage: E = EllipticCurve([0,0,0,0,2*a+1])
            sage: [(p,E.has_additive_reduction(p)) for p in [P17a,P17b]]
            [(Fractional ideal (4*a^2 - 2*a + 1), False),
            (Fractional ideal (2*a + 1), True)]
        """
        return self.local_data(P).has_additive_reduction()

    def tamagawa_number(self, P, proof = None):
        r"""
        Returns the Tamagawa number of this elliptic curve at the prime `P`.

        INPUT:

        - ``P`` -- either None or a prime ideal of the base field of self.

        - ``proof`` -- whether to only use provably correct methods
          (default controlled by global proof module).  Note that the
          proof module is number_field, not elliptic_curves, since the
          functions that actually need the flag are in number fields.

        OUTPUT:

        (positive integer) The Tamagawa number of the curve at `P`.

        EXAMPLES::

            sage: K.<a>=NumberField(x^2-5)
            sage: E=EllipticCurve([20, 225, 750, 625*a + 6875, 31250*a + 46875])
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
            sage: K.<a>=NumberField(x^2+3)
            sage: eK = e.base_extend(K)
            sage: eK.tamagawa_numbers()
            [4, 6, 1]
        """
        return [self.tamagawa_number(p) for p in prime_divisors(self.conductor())]

    def tamagawa_exponent(self, P, proof = None):
        r"""
        Returns the Tamagawa index of this elliptic curve at the prime `P`.

        INPUT:

        - ``P`` -- either None or a prime ideal of the base field of self.

        - ``proof`` -- whether to only use provably correct methods
          (default controlled by global proof module).  Note that the
          proof module is number_field, not elliptic_curves, since the
          functions that actually need the flag are in number fields.

        OUTPUT:

        (positive integer) The Tamagawa index of the curve at P.

        EXAMPLES::

            sage: K.<a>=NumberField(x^2-5)
            sage: E=EllipticCurve([20, 225, 750, 625*a + 6875, 31250*a + 46875])
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

    def tamagawa_product_bsd(self):
        r"""
        Given an elliptic curve `E` over a number field `K`, this function returns the
        integer `C(E/K)` that appears in the Birch and Swinnerton-Dyer conjecture accounting
        for the local information at finite places. If the model is a global minimal model then `C(E/K)` is
        simply the product of the Tamagawa numbers `c_v` where `v` runs over all prime ideals of `K`. Otherwise, if the model has to be changed at a place `v` a correction factor appears.
        The definition is such that `C(E/K)` times the periods at the infinite places is invariant
        under change of the Weierstrass model. See [Ta2] and [Do] for details.

        .. note::

            This definition is slightly different from the definition of ``tamagawa_product``
            for curves defined over `\QQ`. Over the rational number it is always defined to be the product
            of the Tamagawa numbers, so the two definitions only agree when the model is global minimal.

        OUTPUT:

        A rational number

        EXAMPLES::

            sage: K.<i> = NumberField(x^2+1)
            sage: E = EllipticCurve([0,2+i])
            sage: E.tamagawa_product_bsd()
            1

            sage: E = EllipticCurve([(2*i+1)^2,i*(2*i+1)^7])
            sage: E.tamagawa_product_bsd()
            4

        An example where the Neron model changes over K::

            sage: K.<t> = NumberField(x^5-10*x^3+5*x^2+10*x+1)
            sage: E = EllipticCurve(K,'75a1')
            sage: E.tamagawa_product_bsd()
            5
            sage: da = E.local_data()
            sage: [dav.tamagawa_number() for dav in da]
            [1, 1]

        An example over `\mathbb{Q}` (:trac:`9413`)::

            sage: E = EllipticCurve('30a')
            sage: E.tamagawa_product_bsd()
            6

        REFERENCES:

        - [Ta2] Tate, John, On the conjectures of Birch and Swinnerton-Dyer and a geometric analog. Seminaire Bourbaki, Vol. 9, Exp. No. 306.

        - [Do] Dokchitser, Tim and Vladimir, On the Birch-Swinnerton-Dyer quotients modulo squares, Annals of Math., 2010.

        """
        da = self.local_data()
        pr = 1
        for dav in da:
            pp = dav.prime()
            cv = dav.tamagawa_number()
            # uu is the quotient of the Neron differential at pp divided by
            # the differential associated to this particular equation E
            uu = self.isomorphism_to(dav.minimal_model()).u
            if self.base_field().absolute_degree() == 1:
                p = pp.gens_reduced()[0]
                f = 1
                v = valuation(ZZ(uu),p)
            else:
                p = pp.smallest_integer()
                f = pp.residue_class_degree()
                v = valuation(uu,pp)
            uu_abs_val = p**(f*v)
            pr *= cv * uu_abs_val
        return pr

    def kodaira_symbol(self, P, proof = None):
        r"""
        Returns the Kodaira Symbol of this elliptic curve at the prime `P`.

        INPUT:

        - ``P`` -- either None or a prime ideal of the base field of self.

        - ``proof`` -- whether to only use provably correct methods
          (default controlled by global proof module).  Note that the
          proof module is number_field, not elliptic_curves, since the
          functions that actually need the flag are in number fields.

        OUTPUT:

        The Kodaira Symbol of the curve at P, represented as a string.

        EXAMPLES::

            sage: K.<a>=NumberField(x^2-5)
            sage: E=EllipticCurve([20, 225, 750, 625*a + 6875, 31250*a + 46875])
            sage: bad_primes = E.discriminant().support(); bad_primes
            [Fractional ideal (-a), Fractional ideal (7/2*a - 81/2), Fractional ideal (-a - 52), Fractional ideal (2)]
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
        Returns the conductor of this elliptic curve as a fractional
        ideal of the base field.

        OUTPUT:

        (fractional ideal) The conductor of the curve.

        EXAMPLES::

            sage: K.<i>=NumberField(x^2+1)
            sage: EllipticCurve([i, i - 1, i + 1, 24*i + 15, 14*i + 35]).conductor()
            Fractional ideal (21*i - 3)
            sage: K.<a>=NumberField(x^2-x+3)
            sage: EllipticCurve([1 + a , -1 + a , 1 + a , -11 + a , 5 -9*a  ]).conductor()
            Fractional ideal (6*a)

        A not so well known curve with everywhere good reduction::

            sage: K.<a>=NumberField(x^2-38)
            sage: E=EllipticCurve([0,0,0, 21796814856932765568243810*a - 134364590724198567128296995, 121774567239345229314269094644186997594*a - 750668847495706904791115375024037711300])
            sage: E.conductor()
            Fractional ideal (1)

        An example which used to fail (see :trac:`5307`)::

            sage: K.<w>=NumberField(x^2+x+6)
            sage: E=EllipticCurve([w,-1,0,-w-6,0])
            sage: E.conductor()
            Fractional ideal (86304, w + 5898)

        An example raised in :trac:`11346`::

            sage: K.<g> = NumberField(x^2 - x - 1)
            sage: E1 = EllipticCurve(K,[0,0,0,-1/48,-161/864])
            sage: [(p.smallest_integer(),e) for p,e in E1.conductor().factor()]
            [(2, 4), (3, 1), (5, 1)]
        """
        try:
            return self._conductor
        except AttributeError:
            pass

        # Note: for number fields other than QQ we could initialize
        # N=K.ideal(1) or N=OK.ideal(1), which are the same, but for
        # K==QQ it has to be ZZ.ideal(1).
        OK = self.base_ring().ring_of_integers()
        self._conductor = prod([d.prime()**(d.conductor_valuation()) \
                                for d in self.local_data()],\
                               OK.ideal(1))
        return self._conductor

    def global_minimal_model(self, proof = None):
        r"""
        Returns a model of self that is integral, minimal at all primes.

        .. note::

           This is only implemented for class number 1.  In general,
           such a model may or may not exist.

        INPUT:

        - ``proof`` -- whether to only use provably correct methods
          (default controlled by global proof module).  Note that the
          proof module is number_field, not elliptic_curves, since the
          functions that actually need the flag are in number fields.

        OUTPUT:

        A global integral and minimal model.

        EXAMPLES::

            sage: K.<a> = NumberField(x^2-38)
            sage: E = EllipticCurve([0,0,0, 21796814856932765568243810*a - 134364590724198567128296995, 121774567239345229314269094644186997594*a - 750668847495706904791115375024037711300])

            sage: E2 = E.global_minimal_model()
            sage: E2 # random (the global minimal model is not unique)
            Elliptic Curve defined by y^2 + a*x*y + (a+1)*y = x^3 + (a+1)*x^2 + (368258520200522046806318444*a-2270097978636731786720859345)*x + (8456608930173478039472018047583706316424*a-52130038506793883217874390501829588391299) over Number Field in a with defining polynomial x^2 - 38

            sage: E2.local_data()
            []

        See :trac:`11347`::

            sage: K.<g> = NumberField(x^2 - x - 1)
            sage: E = EllipticCurve(K,[0,0,0,-1/48,161/864]).integral_model().global_minimal_model(); E
            Elliptic Curve defined by y^2 + x*y + y = x^3 + x^2 over Number Field in g with defining polynomial x^2 - x - 1
            sage: [(p.norm(), e) for p, e in E.conductor().factor()]
            [(9, 1), (5, 1)]
            sage: [(p.norm(), e) for p, e in E.discriminant().factor()]
            [(-5, 2), (9, 1)]

        See :trac:`14472`, this used not to work over a relative extension::

            sage: K1.<w> = NumberField(x^2+x+1)
            sage: m = polygen(K1)
            sage: K2.<v> = K1.extension(m^2-w+1)
            sage: E = EllipticCurve([0*v,-432])
            sage: E.global_minimal_model()
            Elliptic Curve defined by y^2 + (v+w+1)*y = x^3 + ((6*w-10)*v+16*w+20) over Number Field in v with defining polynomial x^2 - w + 1 over its base field
        """
        if proof is None:
            import sage.structure.proof.proof
            # We use the "number_field" flag because the actual proof dependence is in PARI's number field functions.
            proof = sage.structure.proof.proof.get_flag(None, "number_field")
        K = self.base_ring()
        if K.class_number() != 1:
            raise ValueError, "global minimal models only exist in general for class number 1"

        E = self.global_integral_model()
        primes = E.base_ring()(E.discriminant()).support()
        for P in primes:
            E = E.local_data(P,proof).minimal_model()
        return E._reduce_model()

    def reduction(self,place):
       r"""
       Return the reduction of the elliptic curve at a place of good reduction.

       INPUT:

       - ``place`` -- a prime ideal in the base field of the curve

       OUTPUT:

       An elliptic curve over a finite field, the residue field of the place.

       EXAMPLES::

           sage: K.<i> = QuadraticField(-1)
           sage: EK = EllipticCurve([0,0,0,i,i+3])
           sage: v = K.fractional_ideal(2*i+3)
           sage: EK.reduction(v)
           Elliptic Curve defined by y^2  = x^3 + 5*x + 8 over Residue field of Fractional ideal (2*i + 3)
           sage: EK.reduction(K.ideal(1+i))
           Traceback (most recent call last):
           ...
           ValueError: The curve must have good reduction at the place.
           sage: EK.reduction(K.ideal(2))
           Traceback (most recent call last):
           ...
           ValueError: The ideal must be prime.
           sage: K=QQ.extension(x^2+x+1,"a")
           sage: E=EllipticCurve([1024*K.0,1024*K.0])
           sage: E.reduction(2*K)
           Elliptic Curve defined by y^2 + (abar+1)*y = x^3 over Residue field in abar of Fractional ideal (2)
       """
       K = self.base_field()
       OK = K.ring_of_integers()
       try:
           place = K.ideal(place)
       except TypeError:
           raise TypeError, "The parameter must be an ideal of the base field of the elliptic curve"
       if not place.is_prime():
           raise ValueError, "The ideal must be prime."
       disc = self.discriminant()
       if not K.ideal(disc).valuation(place) == 0:
           local_data=self.local_data(place)
           if local_data.has_good_reduction():
               Fv = OK.residue_field(place)
               return local_data.minimal_model().change_ring(Fv)
           raise ValueError, "The curve must have good reduction at the place."
       Fv = OK.residue_field(place)
       return self.change_ring(Fv)

    def _torsion_bound(self,number_of_places = 20):
        r"""
        An upper bound on the order of the torsion subgroup.

        INPUT:

        - ``number_of_places`` (positive integer, default = 20) -- the
          number of places that will be used to find the bound.

        OUTPUT:

        (integer) An upper bound on the torsion order.

        ALGORITHM:

        An upper bound on the order of the torsion.group of the
        elliptic curve is obtained by counting points modulo several
        primes of good reduction.  Note that the upper bound returned
        by this function is a multiple of the order of the torsion
        group, and in general will be greater than the order.

        EXAMPLES::

            sage: CDB=CremonaDatabase()
            sage: [E._torsion_bound() for E in CDB.iter([14])]
            [6, 6, 6, 6, 6, 6]
            sage: [E.torsion_order() for E in CDB.iter([14])]
            [6, 6, 2, 6, 2, 6]
        """
        E = self
        bound = 0
        k = 0
        K = E.base_field()
        OK = K.ring_of_integers()
        disc = E.discriminant()
        p = Integer(1)
        # runs through primes, decomposes them into prime ideals
        while k < number_of_places :
            p = p.next_prime()
            f = K.primes_above(p)
            # runs through prime ideals above p
            for qq in f:
                fqq = qq.residue_class_degree()
                charqq = qq.smallest_integer()
                # take only places with small residue field (so that the
                # number of points will be small)
                if fqq == 1 or charqq**fqq < 3*number_of_places:
                    # check if the model is integral at the place
                    if min([K.ideal(a).valuation(qq) for a in E.a_invariants()]) >= 0:
                        eqq = qq.ramification_index()
                        # check if the formal group at the place is torsion-free
                        # if so the torsion injects into the reduction
                        if eqq < charqq - 1 and disc.valuation(qq) == 0:
                            Etilda = E.reduction(qq)
                            Npp = Etilda.cardinality()
                            bound = gcd(bound,Npp)
                            if bound == 1:
                                return bound
                            k += 1
        return bound

    def torsion_subgroup(self):
        r"""
        Returns the torsion subgroup of this elliptic curve.

        OUTPUT:

        (``EllipticCurveTorsionSubgroup``) The
        ``EllipticCurveTorsionSubgroup`` associated to this elliptic
        curve.

        EXAMPLES::

            sage: E = EllipticCurve('11a1')
            sage: K.<t>=NumberField(x^4 + x^3 + 11*x^2 + 41*x + 101)
            sage: EK = E.base_extend(K)
            sage: tor = EK.torsion_subgroup()  # long time (3s on sage.math, 2013)
            sage: tor  # long time
            Torsion Subgroup isomorphic to Z/5 + Z/5 associated to the Elliptic Curve defined by y^2 + y = x^3 + (-1)*x^2 + (-10)*x + (-20) over Number Field in t with defining polynomial x^4 + x^3 + 11*x^2 + 41*x + 101
            sage: tor.gens()  # long time
            ((16 : 60 : 1), (t : 1/11*t^3 + 6/11*t^2 + 19/11*t + 48/11 : 1))

        ::

            sage: E = EllipticCurve('15a1')
            sage: K.<t>=NumberField(x^2 + 2*x + 10)
            sage: EK=E.base_extend(K)
            sage: EK.torsion_subgroup()
            Torsion Subgroup isomorphic to Z/4 + Z/4 associated to the Elliptic Curve defined by y^2 + x*y + y = x^3 + x^2 + (-10)*x + (-10) over Number Field in t with defining polynomial x^2 + 2*x + 10

        ::

            sage: E = EllipticCurve('19a1')
            sage: K.<t>=NumberField(x^9-3*x^8-4*x^7+16*x^6-3*x^5-21*x^4+5*x^3+7*x^2-7*x+1)
            sage: EK=E.base_extend(K)
            sage: EK.torsion_subgroup()
            Torsion Subgroup isomorphic to Z/9 associated to the Elliptic Curve defined by y^2 + y = x^3 + x^2 + (-9)*x + (-15) over Number Field in t with defining polynomial x^9 - 3*x^8 - 4*x^7 + 16*x^6 - 3*x^5 - 21*x^4 + 5*x^3 + 7*x^2 - 7*x + 1

        ::

            sage: K.<i> = QuadraticField(-1)
            sage: EK = EllipticCurve([0,0,0,i,i+3])
            sage: EK.torsion_subgroup ()
            Torsion Subgroup isomorphic to Trivial group associated to the Elliptic Curve defined by y^2  = x^3 + i*x + (i+3) over Number Field in i with defining polynomial x^2 + 1
        """
        try:
            return self.__torsion_subgroup
        except AttributeError:
            self.__torsion_subgroup = ell_torsion.EllipticCurveTorsionSubgroup(self)
            return self.__torsion_subgroup

    def torsion_order(self):
        r"""
        Returns the order of the torsion subgroup of this elliptic curve.

        OUTPUT:

        (integer) the order of the torsion subgroup of this elliptic curve.

        EXAMPLES::

            sage: E = EllipticCurve('11a1')
            sage: K.<t> = NumberField(x^4 + x^3 + 11*x^2 + 41*x + 101)
            sage: EK = E.base_extend(K)
            sage: EK.torsion_order()
            25

        ::

            sage: E = EllipticCurve('15a1')
            sage: K.<t> = NumberField(x^2 + 2*x + 10)
            sage: EK = E.base_extend(K)
            sage: EK.torsion_order()
            16

        ::

            sage: E = EllipticCurve('19a1')
            sage: K.<t> = NumberField(x^9-3*x^8-4*x^7+16*x^6-3*x^5-21*x^4+5*x^3+7*x^2-7*x+1)
            sage: EK = E.base_extend(K)
            sage: EK.torsion_order()
            9

        ::

            sage: K.<i> = QuadraticField(-1)
            sage: EK = EllipticCurve([0,0,0,i,i+3])
            sage: EK.torsion_order()
            1
         """
        try:
            return self.__torsion_order
        except AttributeError:
            self.__torsion_order = self.torsion_subgroup().order()
            return self.__torsion_order

    def torsion_points(self):
        r"""
        Returns a list of the torsion points of this elliptic curve.

        OUTPUT:

        (list) A sorted list of the torsion points.

        EXAMPLES::

            sage: E = EllipticCurve('11a1')
            sage: E.torsion_points()
            [(0 : 1 : 0), (5 : -6 : 1), (5 : 5 : 1), (16 : -61 : 1), (16 : 60 : 1)]
            sage: K.<t> = NumberField(x^4 + x^3 + 11*x^2 + 41*x + 101)
            sage: EK = E.base_extend(K)
            sage: EK.torsion_points()  # long time (3s on sage.math, 2013)
            [(16 : 60 : 1),
             (5 : 5 : 1),
             (5 : -6 : 1),
             (16 : -61 : 1),
             (t : 1/11*t^3 + 6/11*t^2 + 19/11*t + 48/11 : 1),
             (-3/55*t^3 - 7/55*t^2 - 2/55*t - 133/55 : 6/55*t^3 + 3/55*t^2 + 25/11*t + 156/55 : 1),
             (-9/121*t^3 - 21/121*t^2 - 127/121*t - 377/121 : -7/121*t^3 + 24/121*t^2 + 197/121*t + 16/121 : 1),
             (5/121*t^3 - 14/121*t^2 - 158/121*t - 453/121 : -49/121*t^3 - 129/121*t^2 - 315/121*t - 207/121 : 1),
             (10/121*t^3 + 49/121*t^2 + 168/121*t + 73/121 : 32/121*t^3 + 60/121*t^2 - 261/121*t - 807/121 : 1),
             (1/11*t^3 - 5/11*t^2 + 19/11*t - 40/11 : -6/11*t^3 - 3/11*t^2 - 26/11*t - 321/11 : 1),
             (14/121*t^3 - 15/121*t^2 + 90/121*t + 232/121 : 16/121*t^3 - 69/121*t^2 + 293/121*t - 46/121 : 1),
             (3/55*t^3 + 7/55*t^2 + 2/55*t + 78/55 : 7/55*t^3 - 24/55*t^2 + 9/11*t + 17/55 : 1),
             (-5/121*t^3 + 36/121*t^2 - 84/121*t + 24/121 : 34/121*t^3 - 27/121*t^2 + 305/121*t + 708/121 : 1),
             (-26/121*t^3 + 20/121*t^2 - 219/121*t - 995/121 : 15/121*t^3 + 156/121*t^2 - 232/121*t + 2766/121 : 1),
             (1/11*t^3 - 5/11*t^2 + 19/11*t - 40/11 : 6/11*t^3 + 3/11*t^2 + 26/11*t + 310/11 : 1),
             (-26/121*t^3 + 20/121*t^2 - 219/121*t - 995/121 : -15/121*t^3 - 156/121*t^2 + 232/121*t - 2887/121 : 1),
             (-5/121*t^3 + 36/121*t^2 - 84/121*t + 24/121 : -34/121*t^3 + 27/121*t^2 - 305/121*t - 829/121 : 1),
             (3/55*t^3 + 7/55*t^2 + 2/55*t + 78/55 : -7/55*t^3 + 24/55*t^2 - 9/11*t - 72/55 : 1),
             (14/121*t^3 - 15/121*t^2 + 90/121*t + 232/121 : -16/121*t^3 + 69/121*t^2 - 293/121*t - 75/121 : 1),
             (t : -1/11*t^3 - 6/11*t^2 - 19/11*t - 59/11 : 1),
             (10/121*t^3 + 49/121*t^2 + 168/121*t + 73/121 : -32/121*t^3 - 60/121*t^2 + 261/121*t + 686/121 : 1),
             (5/121*t^3 - 14/121*t^2 - 158/121*t - 453/121 : 49/121*t^3 + 129/121*t^2 + 315/121*t + 86/121 : 1),
             (-9/121*t^3 - 21/121*t^2 - 127/121*t - 377/121 : 7/121*t^3 - 24/121*t^2 - 197/121*t - 137/121 : 1),
             (-3/55*t^3 - 7/55*t^2 - 2/55*t - 133/55 : -6/55*t^3 - 3/55*t^2 - 25/11*t - 211/55 : 1),
             (0 : 1 : 0)]

        ::

            sage: E = EllipticCurve('15a1')
            sage: K.<t> = NumberField(x^2 + 2*x + 10)
            sage: EK = E.base_extend(K)
            sage: EK.torsion_points()
            [(-7 : -5*t - 2 : 1),
             (-7 : 5*t + 8 : 1),
             (-13/4 : 9/8 : 1),
             (-2 : -2 : 1),
             (-2 : 3 : 1),
             (-t - 2 : -t - 7 : 1),
             (-t - 2 : 2*t + 8 : 1),
             (-1 : 0 : 1),
             (t : t - 5 : 1),
             (t : -2*t + 4 : 1),
             (0 : 1 : 0),
             (1/2 : -5/4*t - 2 : 1),
             (1/2 : 5/4*t + 1/2 : 1),
             (3 : -2 : 1),
             (8 : -27 : 1),
             (8 : 18 : 1)]

        ::

            sage: K.<i> = QuadraticField(-1)
            sage: EK = EllipticCurve(K,[0,0,0,0,-1])
            sage: EK.torsion_points ()
             [(-2 : -3*i : 1), (-2 : 3*i : 1), (0 : -i : 1), (0 : i : 1), (0 : 1 : 0), (1 : 0 : 1)]
         """
        T = self.torsion_subgroup() # make sure it is cached
        return sorted(T.points())           # these are also cached in T

    def rank_bounds(self,verbose=0, lim1=5, lim3=50, limtriv=10, maxprob=20, limbigprime=30):
        r"""
        Returns the lower and upper bounds using :meth:`~simon_two_descent`.
        The results of :meth:`~simon_two_descent` are cached.

        .. NOTE::

            The optional parameters control the Simon two descent algorithm;
            see the documentation of :meth:`~simon_two_descent` for more
            details.

        INPUT:

        - ``verbose`` -- 0, 1, 2, or 3 (default: 0), the verbosity level

        - ``lim1``    -- (default: 5) limit on trivial points on quartics

        - ``lim3``  -- (default: 50) limit on points on ELS quartics

        - ``limtriv`` -- (default: 10) limit on trivial points on elliptic curve

        - ``maxprob`` -- (default: 20)

        - ``limbigprime`` -- (default: 30) to distinguish between
          small and large prime numbers. Use probabilistic tests for
          large primes. If 0, don't use probabilistic tests.

        OUTPUT:

        lower and upper bounds


        .. NOTE::

           For non-quadratic number fields, this code does
           return, but it takes a long time.

        EXAMPLES::

            sage: K.<a> = NumberField(x^2 + 23, 'a')
            sage: E = EllipticCurve(K, '37')
            sage: E == loads(dumps(E))
            True
            sage: E.rank_bounds()
            (2, 2)

        Here is a curve with two-torsion, so here the algorithm gives
        bounds on the rank::

            sage: Qrt5.<rt5>=NumberField(x^2-5)
            sage: E=EllipticCurve([0,5-rt5,0,rt5,0])
            sage: E.rank_bounds()
            (1, 2)

        IMPLEMENTATION:

        Uses Denis Simon's PARI/GP scripts from
        http://www.math.unicaen.fr/~simon/.

        """

        lower, upper, gens = self.simon_two_descent(verbose=verbose,lim1=lim1,lim3=lim3,limtriv=limtriv,maxprob=maxprob,limbigprime=limbigprime)
        return lower,upper

    def rank(self,verbose=0, lim1=5, lim3=50, limtriv=10, maxprob=20, limbigprime=30):
        r"""
        Return the rank of this elliptic curve, if it can be determined.

        .. NOTE::

            The optional parameters control the Simon two descent algorithm;
            see the documentation of :meth:`~simon_two_descent` for more
            details.

        INPUT:

        - ``verbose`` -- 0, 1, 2, or 3 (default: 0), the verbosity level

        - ``lim1``    -- (default: 5) limit on trivial points on quartics

        - ``lim3``  -- (default: 50) limit on points on ELS quartics

        - ``limtriv`` -- (default: 10) limit on trivial points on elliptic curve

        - ``maxprob`` -- (default: 20)

        - ``limbigprime`` -- (default: 30) to distinguish between
          small and large prime numbers. Use probabilistic tests for
          large primes. If 0, don't use probabilistic tests.

        OUTPUT:

        If the upper and lower bounds given by Simon two-descent are
        the same, then the rank has been uniquely identified and we
        return this. Otherwise, we raise a ValueError with an error
        message specifying the upper and lower bounds.

        .. NOTE::

           For non-quadratic number fields, this code does return, but it takes
           a long time.

        EXAMPLES::

            sage: K.<a> = NumberField(x^2 + 23, 'a')
            sage: E = EllipticCurve(K, '37')
            sage: E == loads(dumps(E))
            True
            sage: E.rank()
            2

        Here is a curve with two-torsion, so here the bounds given by the
        algorithm do not uniquely determine the rank::

            sage: Qrt5.<rt5>=NumberField(x^2-5)
            sage: E=EllipticCurve([0,5-rt5,0,rt5,0])
            sage: E.rank()
            Traceback (most recent call last):
            ...
            ValueError: There is insufficient data to determine the rank -
            2-descent gave lower bound 1 and upper bound 2

        IMPLEMENTATION:

        Uses Denis Simon's PARI/GP scripts from
        http://www.math.unicaen.fr/~simon/.

        """

        lower,upper = self.rank_bounds(verbose=verbose,lim1=lim1,lim3=lim3,limtriv=limtriv,maxprob=maxprob,limbigprime=limbigprime)
        if lower == upper:
            return lower
        else:
            raise ValueError, 'There is insufficient data to determine the rank - 2-descent gave lower bound %s and upper bound %s' % (lower, upper)

    def gens(self,verbose=0, lim1=5, lim3=50, limtriv=10, maxprob=20, limbigprime=30):
        r"""
        Returns some generators of this elliptic curve. Check :meth:`~rank` or
        :meth:`~rank_bounds` to verify the number of generators.

        .. NOTE::

            The optional parameters control the Simon two descent algorithm;
            see the documentation of :meth:`~simon_two_descent` for more
            details.

        INPUT:

        - ``verbose`` -- 0, 1, 2, or 3 (default: 0), the verbosity level

        - ``lim1``    -- (default: 5) limit on trivial points on quartics

        - ``lim3``  -- (default: 50) limit on points on ELS quartics

        - ``limtriv`` -- (default: 10) limit on trivial points on elliptic curve

        - ``maxprob`` -- (default: 20)

        - ``limbigprime`` -- (default: 30) to distinguish between
          small and large prime numbers. Use probabilistic tests for
          large primes. If 0, don't use probabilistic tests.

        OUTPUT:

        The linearly independent elements given by the Simon two-descent.

        .. NOTE::

           For non-quadratic number fields, this code does return, but it takes
           a long time.

        EXAMPLES::

            sage: K.<a> = NumberField(x^2 + 23, 'a')
            sage: E = EllipticCurve(K, '37')
            sage: E == loads(dumps(E))
            True
            sage: E.gens()
            [(-1 : 0 : 1), (1/2*a - 5/2 : -1/2*a - 13/2 : 1)]


        Here is a curve with two-torsion, so here the algorithm does not
        uniquely determine the rank::

            sage: Qrt5.<rt5>=NumberField(x^2-5)
            sage: E=EllipticCurve([0,5-rt5,0,rt5,0])
            sage: E.gens()
            [(3/2*rt5 + 5/2 : -9/2*rt5 - 15/2 : 1), (-1/2*rt5 + 3/2 : 3/2*rt5 - 9/2 : 1), (0 : 0 : 1)]

        IMPLEMENTATION:

        Uses Denis Simon's PARI/GP scripts from
        http://www.math.unicaen.fr/~simon/.
        """

        lower,upper,gens = self.simon_two_descent(verbose=verbose,lim1=lim1,lim3=lim3,limtriv=limtriv,maxprob=maxprob,limbigprime=limbigprime)
        return gens


    def period_lattice(self, embedding):
        r"""
        Returns the period lattice of the elliptic curve for the given
        embedding of its base field with respect to the differential
        `dx/(2y + a_1x + a_3)`.

        INPUT:

        - ``embedding`` - an embedding of the base number field into `\RR` or `\CC`.

        .. note::

           The precision of the embedding is ignored: we only use the
           given embedding to determine which embedding into ``QQbar``
           to use.  Once the lattice has been initialized, periods can
           be computed to arbitrary precision.


        EXAMPLES:

        First define a field with two real embeddings::

            sage: K.<a> = NumberField(x^3-2)
            sage: E=EllipticCurve([0,0,0,a,2])
            sage: embs=K.embeddings(CC); len(embs)
            3

        For each embedding we have a different period lattice::

            sage: E.period_lattice(embs[0])
            Period lattice associated to Elliptic Curve defined by y^2 = x^3 + a*x + 2 over Number Field in a with defining polynomial x^3 - 2 with respect to the embedding Ring morphism:
            From: Number Field in a with defining polynomial x^3 - 2
            To:   Algebraic Field
            Defn: a |--> -0.6299605249474365? - 1.091123635971722?*I

            sage: E.period_lattice(embs[1])
            Period lattice associated to Elliptic Curve defined by y^2 = x^3 + a*x + 2 over Number Field in a with defining polynomial x^3 - 2 with respect to the embedding Ring morphism:
            From: Number Field in a with defining polynomial x^3 - 2
            To:   Algebraic Field
            Defn: a |--> -0.6299605249474365? + 1.091123635971722?*I

            sage: E.period_lattice(embs[2])
            Period lattice associated to Elliptic Curve defined by y^2 = x^3 + a*x + 2 over Number Field in a with defining polynomial x^3 - 2 with respect to the embedding Ring morphism:
            From: Number Field in a with defining polynomial x^3 - 2
            To:   Algebraic Field
            Defn: a |--> 1.259921049894873?

        Although the original embeddings have only the default
        precision, we can obtain the basis with higher precision
        later::

            sage: L=E.period_lattice(embs[0])
            sage: L.basis()
            (1.86405007647981 - 0.903761485143226*I, -0.149344633143919 - 2.06619546272945*I)

            sage: L.basis(prec=100)
            (1.8640500764798108425920506200 - 0.90376148514322594749786960975*I, -0.14934463314391922099120107422 - 2.0661954627294548995621225062*I)
        """
        from sage.schemes.elliptic_curves.period_lattice import PeriodLattice_ell
        return PeriodLattice_ell(self,embedding)

    def is_isogenous(self, other, proof=True, maxnorm=100):
        """
        Returns whether or not self is isogenous to other.

        INPUT:

        - ``other`` -- another elliptic curve.

        - ``proof`` (default True) -- If ``False``, the function will
          return ``True`` whenever the two curves have the same
          conductor and are isogenous modulo `p` for all primes `p` of
          norm up to ``maxp``.  If ``True``, the function returns
          False when the previous condition does not hold, and if it
          does hold we attempt to see if the curves are indeed
          isogenous.  However, this has not been fully implemented
          (see examples below), so we may not be able to determine
          whether or not the curves are isogenous..

        - ``maxnorm`` (integer, default 100) -- The maximum norm of
          primes `p` for which isogeny modulo `p` will be checked.

        OUTPUT:

        (bool) True if there is an isogeny from curve ``self`` to
        curve ``other``.

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
            sage: F = NumberField(x^2 -2, 's'); F
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
            sage: F = NumberField(x^2 -2, 's'); F
            Number Field in s with defining polynomial x^2 - 2
            sage: E = EllipticCurve('14a1')
            sage: EE = EllipticCurve('14a2')
            sage: E1 = E.change_ring(F)
            sage: E2 = EE.change_ring(F)
            sage: E1.is_isogenous(E2)
            True

        ::

            sage: x = polygen(QQ, 'x')
            sage: F = NumberField(x^2 -2, 's'); F
            Number Field in s with defining polynomial x^2 - 2
            sage: k.<a> = NumberField(x^3+7)
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
            sage: E1.is_isogenous(E2)
            Traceback (most recent call last):
            ...
            NotImplementedError: Curves appear to be isogenous (same conductor, isogenous modulo all primes of norm up to 1000), but no isogeny has been constructed.
            sage: E1.is_isogenous(E2, proof=False)
            True

        In this case E1 and E2 are in fact 9-isogenous, as may be
        deduced from the following::

           sage: E3 = EllipticCurve([i + 1, 0, 1, -5*i - 5, -2*i - 5])
           sage: E3.is_isogenous(E1)
           True
           sage: E3.is_isogenous(E2)
           True
           sage: E1.isogeny_degree(E2)
           9

        """
        if not is_EllipticCurve(other):
            raise ValueError, "Second argument is not an Elliptic Curve."
        if self.is_isomorphic(other):
            return True
        K = self.base_field()
        if K != other.base_field():
            raise ValueError, "Second argument must be defined over the same number field."

        E1 = self.integral_model()
        E2 = other.integral_model()
        N = E1.conductor()
        if N != E2.conductor():
            return False

        PI = K.primes_of_degree_one_iter()
        while True:
            P = PI.next()
            if P.norm() > maxnorm: break
            if not P.divides(N):
                OP = K.residue_field(P)
                if E1.change_ring(OP).cardinality() != E2.change_ring(OP).cardinality():
                    return False

        if not proof:
            return True

        # We have not yet implemented isogenies of all possible
        # degrees, and do not know a bound on the possible degrees
        # over general number fields.  But here we do at least try
        # some easy cases:

        for l in [2,3,5,7,13]:
            if any([E2.is_isomorphic(f.codomain()) for f in E1.isogenies_prime_degree(l)]):
                return True

        # Next we try looking modulo some more primes:

        while True:
            if P.norm() > 10*maxnorm: break
            if not P.divides(N):
                OP = K.residue_field(P)
                if E1.change_ring(OP).cardinality() != E2.change_ring(OP).cardinality():
                    return False
            P = PI.next()

        # At this point is is highly likely that the curves are
        # isogenous, but we have not proved it.

        raise NotImplementedError, "Curves appear to be isogenous (same conductor, isogenous modulo all primes of norm up to %s), but no isogeny has been constructed." % (10*maxnorm)

    def isogeny_degree(self, other):
        """
        Returns the minimal degree of an isogeny between self and
        other, or 0 if no isogeny exists.

        INPUT:

        - ``other`` -- another elliptic curve.

        OUTPUT:

        (int) The degree of an isogeny from ``self`` to ``other``, or 0.

        .. warning::

           Not all isogenies over number fields are yet implemented.
           Currently the code only works if there is a chain of
           isogenies from ``self`` to ``other`` of degrees 2, 3, 5, 7
           and 13.

        EXAMPLES::

            sage: x = QQ['x'].0
            sage: F = NumberField(x^2 -2, 's'); F
            Number Field in s with defining polynomial x^2 - 2
            sage: E = EllipticCurve('14a1')
            sage: EE = EllipticCurve('14a2')
            sage: E1 = E.change_ring(F)
            sage: E2 = EE.change_ring(F)
            sage: E1.isogeny_degree(E2)
            2
            sage: E2.isogeny_degree(E2)
            1
            sage: E5 = EllipticCurve('14a5').change_ring(F)
            sage: E1.isogeny_degree(E5)
            6
        """
        if self.conductor() != other.conductor():
            return Integer(0)

        if self.is_isomorphic(other):
            return Integer(1)

        from sage.sets.set import Set

        curves = [self]
        degrees = [Integer(1)]
        l_list = [l for l in Set([ZZ(f.degree()) for f in self.isogenies_prime_degree([2,3,5,7,13])])]

        newcurves = []
        newdegs = []
        k = 0
        while k<len(curves):
            newcurves.extend([f.codomain() for f in curves[k].isogenies_prime_degree(l_list)])
            newdegs.extend([degrees[k]*f.degree() for f in curves[k].isogenies_prime_degree(l_list)])
            newisogpairs = dict(zip(newcurves, newdegs))
            i = 0
            while i<len(curves):
                j = 0
                while j<len(newcurves):
                    if curves[i].is_isomorphic(newcurves[j]):
                        newdegs.remove(newisogpairs[newcurves[j]])
                        newcurves.remove(newcurves[j])
                    else:
                        j = j+1
                i = i+1

            m = 0
            newisogpairs = dict(zip(newcurves, newdegs))
            while m<len(newcurves):
                if other.is_isomorphic(newcurves[m]):
                    return newisogpairs[newcurves[m]]
                m = m+1

            curves.extend(newcurves)
            degrees.extend(newdegs)
            k = k+1

        raise NotImplementedError("Not all isogenies implemented over general number fields.")

    def lll_reduce(self, points, height_matrix=None, precision=None):
        """
        Returns an LLL-reduced basis from a given basis, with transform
        matrix.

        INPUT:

        - ``points`` - a list of points on this elliptic
          curve, which should be independent.

        - ``height_matrix`` - the height-pairing matrix of
          the points, or ``None``. If ``None``, it will be computed.

        - ``precision`` - number of bits of precision of intermediate
          computations (default: ``None``, for default RealField
          precision; ignored if ``height_matrix`` is supplied)

        OUTPUT: A tuple (newpoints, U) where U is a unimodular integer
        matrix, new_points is the transform of points by U, such that
        new_points has LLL-reduced height pairing matrix

        .. note::

            If the input points are not independent, the output
            depends on the undocumented behaviour of PARI's
            ``qflllgram()`` function when applied to a gram matrix which
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
            sage: xi = [2005024558054813068,\
            -4690836759490453344,\
            4700156326649806635,\
            6785546256295273860,\
            6823803569166584943,\
            7788809602110240789,\
            27385442304350994620556,\
            54284682060285253719/4,\
            -94200235260395075139/25,\
            -3463661055331841724647/576,\
            -6684065934033506970637/676,\
            -956077386192640344198/2209,\
            -27067471797013364392578/2809,\
            -25538866857137199063309/3721,\
            -1026325011760259051894331/108241,\
            9351361230729481250627334/1366561,\
            10100878635879432897339615/1423249,\
            11499655868211022625340735/17522596,\
            110352253665081002517811734/21353641,\
            414280096426033094143668538257/285204544,\
            36101712290699828042930087436/4098432361,\
            45442463408503524215460183165/5424617104,\
            983886013344700707678587482584/141566320009,\
            1124614335716851053281176544216033/152487126016]
            sage: points = [E.lift_x(x) for x in xi]
            sage: newpoints, U = E.lll_reduce(points)  # long time (35s on sage.math, 2011)
            sage: [P[0] for P in newpoints]            # long time
            [6823803569166584943, 5949539878899294213, 2005024558054813068, 5864879778877955778, 23955263915878682727/4, 5922188321411938518, 5286988283823825378, 175620639884534615751/25, -11451575907286171572, 3502708072571012181, 1500143935183238709184/225, 27180522378120223419/4, -5811874164190604461581/625, 26807786527159569093, 7404442636649562303, 475656155255883588, 265757454726766017891/49, 7272142121019825303, 50628679173833693415/4, 6951643522366348968, 6842515151518070703, 111593750389650846885/16, 2607467890531740394315/9, -1829928525835506297]

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

        Some examples over number fields (see :trac:`9411`)::

            sage: K.<a> = QuadraticField(-23, 'a')
            sage: E = EllipticCurve(K, '37')
            sage: E.lll_reduce(E.gens())
            (
                                                    [1 1]
            [(-1 : 0 : 1), (-2 : 1/2*a - 1/2 : 1)], [0 1]
            )

        ::

            sage: K.<a> = QuadraticField(-5)
            sage: E = EllipticCurve(K,[0,a])
            sage: points = [E.point([-211/841*a - 6044/841,-209584/24389*a + 53634/24389]),E.point([-17/18*a - 1/9, -109/108*a - 277/108]) ]
            sage: E.lll_reduce(points)
            (
            [(-a + 4 : -3*a + 7 : 1), (-17/18*a - 1/9 : 109/108*a + 277/108 : 1)],
            [ 1  0]
            [ 1 -1]
            )
        """
        r = len(points)
        if height_matrix is None:
            height_matrix = self.height_pairing_matrix(points, precision)
        U = pari(height_matrix).lllgram().python()
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

            sage: K = NumberField(x**2 + 1, 'a')
            sage: E = EllipticCurve('11a1').change_ring(K)
            sage: rho = E.galois_representation()
            sage: rho
            Compatible family of Galois representations associated to the Elliptic Curve defined by y^2 + y = x^3 + (-1)*x^2 + (-10)*x + (-20) over Number Field in a with defining polynomial x^2 + 1
            sage: rho.is_surjective(3)
            True
            sage: rho.is_surjective(5)  # long time (9s on sage.math, 2013)
            False
            sage: rho.non_surjective()
            [5]
        """
        return gal_reps_number_field.GaloisRepresentation(self)
