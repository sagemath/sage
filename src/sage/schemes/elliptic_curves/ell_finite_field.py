"""
Elliptic curves over finite fields

AUTHORS:

- William Stein (2005): Initial version

- Robert Bradshaw et al....

- John Cremona (2008-02): Point counting and group structure for
  non-prime fields, Frobenius endomorphism and order, elliptic logs

- Mariah Lenox (2011-03): Added set_order method
"""

#*****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
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

from sage.misc.randstate import current_randstate

from sage.schemes.plane_curves.projective_curve import Hasse_bounds
from ell_field import EllipticCurve_field
from constructor import EllipticCurve, EllipticCurve_from_j
from sage.schemes.hyperelliptic_curves.hyperelliptic_finite_field import HyperellipticCurve_finite_field
import sage.rings.ring as ring
from sage.rings.all import Integer, ZZ, PolynomialRing, GF, polygen
from sage.rings.finite_rings.element_base import is_FiniteFieldElement
import sage.groups.generic as generic
import ell_point
from sage.rings.arith import gcd, lcm
from sage.structure.sequence import Sequence

import sage.plot.all as plot

import sage.libs.pari
pari = sage.libs.pari.all.pari

class EllipticCurve_finite_field(EllipticCurve_field, HyperellipticCurve_finite_field):
    """
    Elliptic curve over a finite field.

    EXAMPLES::

        sage: EllipticCurve(GF(101),[2,3])
        Elliptic Curve defined by y^2  = x^3 + 2*x + 3 over Finite Field of size 101

        sage: F=GF(101^2, 'a')
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
        Return a list of all the points on the curve, for prime fields only
        (see points() for the general case)

        EXAMPLES::

            sage: S=EllipticCurve(GF(97),[2,3])._points_via_group_structure()
            sage: len(S)
            100

        See trac #4687, where the following example did not work::

            sage: E=EllipticCurve(GF(2),[0, 0, 1, 1, 1])
            sage: E.points()
            [(0 : 1 : 0)]

        ::

            sage: E=EllipticCurve(GF(2),[0, 0, 1, 0, 1])
            sage: E.points()
            [(0 : 1 : 0), (1 : 0 : 1), (1 : 1 : 1)]

        ::

            sage: E=EllipticCurve(GF(4,'a'),[0, 0, 1, 0, 1])
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
        for m in range(1,ni[0]): H0.append(H0[-1]+pts[0])
        if ngens == 1:    # cyclic group
            return H0

        # else noncyclic group
        H1=[self(0)]
        for m in range(1,ni[1]): H1.append(H1[-1]+pts[1])
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

            sage: K = GF(p**2,'a')
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
        except AttributeError: pass

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
        Returns the cardinality of this elliptic curve over the base field or extensions.

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

        if n<1:
            raise ValueError("n must be a positive integer")

        if n==1:
            return self.cardinality()
        else:
            return [ self.cardinality(extension_degree=i) for  i in range(1,n+1)]

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
            sage: P = E.random_element(); P
            (16740 : 12486 : 1)
            sage: type(P)
            <class 'sage.schemes.elliptic_curves.ell_point.EllipticCurvePoint_finite_field'>
            sage: P in E
            True

        ::

            sage: k.<a> = GF(7^5)
            sage: E = EllipticCurve(k,[2,4])
            sage: P = E.random_element(); P
            (5*a^4 + 3*a^3 + 2*a^2 + a + 4 : 2*a^4 + 3*a^3 + 4*a^2 + a + 5 : 1)
            sage: type(P)
            <class 'sage.schemes.elliptic_curves.ell_point.EllipticCurvePoint_finite_field'>
            sage: P in E
            True

        ::

            sage: k.<a> = GF(2^5)
            sage: E = EllipticCurve(k,[a^2,a,1,a+1,1])
            sage: P = E.random_element(); P
            (a^4 + a : a^4 + a^3 + a^2 : 1)
            sage: type(P)
            <class 'sage.schemes.elliptic_curves.ell_point.EllipticCurvePoint_finite_field'>
            sage: P in E
            True

        Ensure that the entire point set is reachable::

            sage: E = EllipticCurve(GF(11), [2,1])
            sage: len(set(E.random_element() for _ in range(100)))
            16
            sage: E.cardinality()
            16

        TESTS:

        See trac #8311::

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

        .. note::

           This computes the curve cardinality, which may be
           time-consuming.

        EXAMPLES::

            sage: E=EllipticCurve(GF(101),[2,3])
            sage: E.trace_of_frobenius()
            6
            sage: E=EllipticCurve(GF(11^5,'a'),[2,5])
            sage: E.trace_of_frobenius()
            802

        The following shows that the issue from trac #2849 is fixed::

            sage: E=EllipticCurve(GF(3^5,'a'),[-1,-1])
            sage: E.trace_of_frobenius()
            -27
        """
        return 1 + self.base_field().order() - self.cardinality()

    def _cardinality_with_j_invariant_1728(self):
        r"""
        Special function to compute cardinality when j=1728.

        EXAMPLES: An example with q=p=1 (mod 4)

        ::

            sage: F=GF(10009)
            sage: [EllipticCurve(F,[0,0,0,11^i,0])._cardinality_with_j_invariant_1728() for i in range(4)]
            [10016, 10210, 10004, 9810]

        An example with q=p=3 (mod 4)

        ::

            sage: F=GF(10007)
            sage: [EllipticCurve(F,[0,0,0,5^i,0])._cardinality_with_j_invariant_1728() for i in range(4)]
            [10008, 10008, 10008, 10008]

        An example with `q=p^2`, p=3 (mod 4)

        ::

            sage: F.<a>=GF(10007^2,'a')
            sage: [EllipticCurve(F,[0,0,0,a^i,0])._cardinality_with_j_invariant_1728() for i in range(4)]
            [100160064, 100140050, 100120036, 100140050]

        Examples with `q=2^d`, d odd (3 isomorphism classes)::

            sage: F.<a> = GF(2**15,'a')
            sage: ais = [[0,0,1,0,0],[0,0,1,1,0],[0,0,1,1,1]]
            sage: curves=[EllipticCurve(F,ai) for ai in ais]
            sage: all([all([e1==e2 or not e1.is_isomorphic(e2) for e1 in curves]) for e2 in curves])
            True
            sage: [e._cardinality_with_j_invariant_1728() for e in curves]
            [32769, 33025, 32513]

        Examples with `q=2^d`, d even (7 isomorphism classes)::

            sage: F.<a> = GF(2**16,'a')
            sage: b = a^11 # trace 1
            sage: ais = [[0,0,1,0,0],[0,0,1,0,b],[0,0,1,b,0],[0,0,a,0,0],[0,0,a,0,a^2*b],[0,0,a^2,0,0],[0,0,a^2,0,a^4*b]]
            sage: curves=[EllipticCurve(F,ai) for ai in ais]
            sage: all([all([e1==e2 or not e1.is_isomorphic(e2) for e1 in curves]) for e2 in curves])
            True
            sage: [e._cardinality_with_j_invariant_1728() for e in curves]
            [65025, 66049, 65537, 65793, 65281, 65793, 65281]

        Examples with `q=3^d`, d odd (4 isomorphism classes)::

            sage: F.<a> = GF(3**15,'a')
            sage: b=a^7  # has trace 1
            sage: ais=[[0,0,0,1,0],[0,0,0,-1,0],[0,0,0,-1,b],[0,0,0,-1,-b]]
            sage: curves=[EllipticCurve(F,ai) for ai in ais]
            sage: all([all([e1==e2 or not e1.is_isomorphic(e2) for e1 in curves]) for e2 in curves])
            True
            sage: [e._cardinality_with_j_invariant_1728() for e in curves]
            [14348908, 14348908, 14342347, 14355469]

        Examples with `q=3^d`, d even (6 isomorphism classes)::

            sage: F.<g>=GF(3^18,'g')
            sage: i=F(-1).sqrt()
            sage: a=g^8  # has trace 1
            sage: ais= [[0,0,0,1,0],[0,0,0,1,i*a],[0,0,0,g,0],[0,0,0,g^3,0],[0,0,0,g^2,0], [0,0,0,g^2,i*a*g^3]]
            sage: curves=[EllipticCurve(F,ai) for ai in ais]
            sage: all([all([e1==e2 or not e1.is_isomorphic(e2) for e1 in curves]) for e2 in curves])
            True
            sage: [E._cardinality_with_j_invariant_1728() for E in curves]
            [387459856, 387400807, 387420490, 387420490, 387381124, 387440173]

        TESTS:

        Check that a bug noted at :trac:`15667` is fixed::

            sage: F.<a>=GF(3^6,'a')
            sage: EllipticCurve([a^5 + 2*a^3 + 2*a^2 + 2*a, a^4 + a^3 + 2*a + 1]).cardinality()
            784
            sage: EllipticCurve([a^5 + 2*a^3 + 2*a^2 + 2*a, a^4 + a^3 + 2*a + 1]).cardinality_exhaustive()
            784

        """
        try:
            return self._order
        except AttributeError:
            pass

        k = self.base_ring()
        assert self.j_invariant()==k(1728)
        q = k.cardinality()
        p = k.characteristic()
        d = k.degree()
        x=polygen(ZZ)

# p=2, j=0=1728
#
# Number of isomorphism classes is 3 (in odd degree) or 7 (in even degree)
#
        if p==2:
            if d%2==1:
                # The 3 classes are represented, independently of d,
                # by [0,0,1,0,0], [0,0,1,1,0], [0,0,1,1,1]
                E=EllipticCurve(k,[0,0,1,0,0])
                if self.is_isomorphic(E):
                    N = q+1
                else:
                    n = (d+1)//2
                    t = 2**n
                    n = n%4
                    if n==0 or n==1: t=-t
                    E=EllipticCurve(k,[0,0,1,1,1])
                    if self.is_isomorphic(E): t=-t
                    N = q+1-t
            else:
                # The 7 classes are represented by E1=[0,0,1,0,0],
                # E2=[0,0,1,0,b], E3=[0,0,1,b,0], E4=[0,0,a,0,0],
                # E4=[0,0,a,0,a^2*b], E6=[0,0,a^2,0,0],
                # E7=[0,0,a^2,0,a^4*b], where a is a non-cube and b
                # has trace 1.  E1's Frobenius is pi=(-2)**(d//2); the
                # Frobeniuses are then pi, -pi, 0; w*pi, -w*pi;
                # w^2*pi, -w^2*pi where w is either cube root of
                # unity, so the traces are 2*pi, -2*pi, 0, -pi, +pi;
                # -pi, +pi.
                delta = self.discriminant()
                discube = (delta**((q-1)//3) == k(1))
                pi = (-2)**(d//2)
                if discube:
                    a = k.gen()
                    b = a
                    while b.trace()==0: b*=a
                    if self.is_isomorphic(EllipticCurve(k,[0,0,1,b,0])):
                        t = 0
                    else:
                        t = 2*pi
                        if not self.is_isomorphic(EllipticCurve(k,[0,0,1,0,0])):
                            t = -t

                else:
                    t = pi
                    if self.is_isomorphic(EllipticCurve(k,[0,0,delta,0,0])):
                        t = -t
                N = q+1-t


# p=3, j=0=1728
#
# Number of isomorphism classes is 4 (odd degree) or 6 (even degree)
#
        elif p==3:
            if d%2==1:
                # The 4 classes are represented by [0,0,0,1,0],
                # [0,0,0,-1,0], [0,0,0,-1,a], [0,0,0,-1,-a] where a
                # has trace 1
                delta = self.discriminant()
                if (-delta).is_square():
                    t = 0
                else:
                    u = delta.sqrt()
                    if not u.is_square(): u=-u
                    tr = ((self.a3()**2+self.a6())/u).trace()
                    if tr==0:
                        t = 0
                    else:
                        d2 = (d+1)//2
                        t = 3**d2
                        if d2%2==1: t = -t
                        if tr==-1:  t = -t
                N = q+1-t
            else:
                # The 6 classes are represented by: [0,0,0,1,0],
                # [0,0,0,1,i*a]; [0,0,0,g,0], [0,0,0,g^3,0];
                # [0,0,0,g^2,0], [0,0,0,g^2,i*a*g^3]; where g
                # generates the multiplicative group modulo 4th
                # powers, and a has nonzero trace.

                # The curve is isomorphic to [0,0,0,A4,A6]

                A4 = self.a4() - self.a1()*self.a3() # = -b4 = 2*b4
                if A4.is_square():
                    u = A4.sqrt()
                    t = (-3)**(d//2)
                    i = k(-1).sqrt()
                    A6 = self.a3()**2 + self.a6()   # = b6
                    if (A6/(i*u*A4)).trace()==0:
                        t *= 2
                    else:
                        t *= -1
                    if not u.is_square():
                        t *= -1
                else:
                    t = 0
                N = q+1-t

# p>3, j=1728
#
# Number of isomorphism classes is 4 if q=1 (mod 4), else 2
#
        elif p%4==3:
            if d%2==1:
                t = 0
            else:
                t  = (-p)**(d//2)
                w = (self.c4()/k(48))**((q-1)//4)
                if   w==1:    t =  2*t
                elif w==-1:   t = -2*t
                else: t = 0

            N = q+1-t

# p=1 (mod 4).  First find Frobenius pi=a+b*i for [0,0,0,-1,0] over GF(p):
# N(pi)=p and N(pi-1)=0 (mod 8).
#
        else: # p%4==1
            R = ZZ.extension(x**2+1,'i')
            i = R.gen(1)
            pi = R.fraction_field().factor(p)[0][0].gens_reduced()[0]
            a,b = pi.list()
            if a%2==0:
                a,b = -b,a
            if (a+b+1)%4==0:
                a,b = -a,-b
            pi = a+b*i        # Now pi=a+b*i with (a,b)=(1,0),(3,2) mod 4

        # Lift to Frobenius for [0,0,0,-1,0] over GF(p^d):

            if d>1:
                pi = pi**d
                a,b = pi.list()

        # Compute appropriate quartic twist:

            w = (self.c4()/k(48))**((q-1)//4)
            if w==1:
                t = 2*a
            elif w==-1:
                t = -2*a
            elif k(b)==w*k(a):
                t = 2*b
            else:
                t = -2*b
            N = q+1-t

        self._order = Integer(N)
        return self._order

    def _cardinality_with_j_invariant_0(self):
        r"""
        Special function to compute cardinality when j=0.

        EXAMPLES: An example with q=p=1 (mod 6)

        ::

            sage: F=GF(1009)
            sage: [EllipticCurve(F,[0,0,0,0,11^i])._cardinality_with_j_invariant_0() for i in range(6)]
            [948, 967, 1029, 1072, 1053, 991]

        An example with q=p=5 (mod 6)

        ::

            sage: F=GF(1013)
            sage: [EllipticCurve(F,[0,0,0,0,3^i])._cardinality_with_j_invariant_0() for i in range(6)]
            [1014, 1014, 1014, 1014, 1014, 1014]

        An example with `q=p^2`, p=5 (mod 6)

        ::

            sage: F.<a>=GF(1013^2,'a')
            sage: [EllipticCurve(F,[0,0,0,0,a^i])._cardinality_with_j_invariant_0() for i in range(6)]
            [1028196, 1027183, 1025157, 1024144, 1025157, 1027183]

        For examples in characteristic 2 and 3, see the function
        _cardinality_with_j_invariant_1728()
        """

        try:
            return self._order
        except AttributeError:
            pass

        k = self.base_ring()
        assert self.j_invariant()==k(0)
        p = k.characteristic()
        if p==2 or p==3:  # then 0==1728
            return self._cardinality_with_j_invariant_1728()

        q = k.cardinality()
        d = k.degree()
        x=polygen(ZZ)

# p>3, j=0
#
# Number of isomorphism classes is 4 if q=1 (mod 4), else 2
#
        if p%6==5:
            if d%2==1:
                t = 0
            else:
                t  = (-p)**(d//2)
                w = (self.c6()/k(-864))**((q-1)//6)
                if   w==1:    t =  2*t
                elif w==-1:   t = -2*t
                elif w**3==1: t = -t

            N = q+1-t

# p=1 (mod 6).  First find Frobenius pi=a+b*w for [0,0,0,0,1] over GF(p):
# N(pi)=p and N(pi-1)=0 (mod 12).
#
        else: # p%6==1
            R = ZZ.extension(x**2-x+1,'zeta6')
            zeta6 = R.gen(1)
            pi = R.fraction_field().factor(p)[0][0].gens_reduced()[0]
            while (pi-1).norm()%12 !=0:  pi*=zeta6
            a,b = pi.list()
            z = k(-b)/k(a)  # a *specific* 6th root of unity in k

# Now pi=a+b*zeta6 with N(pi-1)=0 (mod 12)

# Lift to Frobenius for [0,0,0,0,1] over GF(p^d):

            if d>1:
                pi = pi**d
                a,b = pi.list()

# Compute appropriate sextic twist:

            w = (self.c6()/k(-864))**((q-1)//6)

            if   w==1:       t = 2*a+b  # = Trace(pi)
            elif w==-1:      t = -2*a-b # = Trace(-pi)
            elif w==z:       t = a-b    # = Trace(pi*zeta6)
            elif w==z**2:    t = -a-2*b # = Trace(pi*zeta6**2)
            elif w==z**4:    t = b-a    # = Trace(pi*zeta6**4)
            elif w==z**5:    t = a+2*b  # = Trace(pi*zeta6**5)

            N = q+1-t

        self._order = Integer(N)
        return self._order

    def cardinality(self, algorithm='pari', extension_degree=1):
        r"""
        Return the number of points on this elliptic curve.

        INPUT:

        - ``algorithm`` -- string (default: ``'pari'``), used only for
          point counting over prime fields:

          - ``'pari'`` -- use the baby-step giant-step or
            Schoof-Elkies-Atkin methods as implemented in the PARI
            C-library function ``ellap``

          - ``'bsgs'`` -- use the baby-step giant-step method as
            implemented in Sage, with the Cremona-Sutherland version
            of Mestre's trick

          - ``'all'`` -- compute cardinality with both ``'pari'`` and
            ``'bsgs'``; return result if they agree or raise a
            ``RuntimeError`` if they do not

        - ``extension_degree`` -- an integer `d` (default: 1): if the
          base field is `\GF{q}`, return the cardinality of ``self``
          over the extension `\GF{q^d}` of degree `d`.

        OUTPUT:

        The order of the group of rational points of ``self`` over its
        base field, or over an extension field of degree `d` as above.
        The result is cached.

        Over prime fields, one of the above algorithms is used. Over
        non-prime fields, the serious point counting is done on a standard
        curve with the same `j`-invariant over the field `\GF{p}(j)`, then
        lifted to the base field, and finally account is taken of twists.

        For `j = 0` and `j = 1728` special formulas are used instead.

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
            sage: E = EllipticCurve(GF(11^2, 'a'), [3,3])
            sage: E.cardinality()
            128
            sage: EllipticCurve(GF(11^100, 'a'), [3,3]).cardinality()
            137806123398222701841183371720896367762643312000384671846835266941791510341065565176497846502742959856128

        TESTS::

            sage: EllipticCurve(GF(10009), [1,2,3,4,5]).cardinality(algorithm='foobar')
            Traceback (most recent call last):
            ...
            ValueError: Algorithm is not known

        If the cardinality has already been computed, then the ``algorithm``
        keyword is ignored::

            sage: E = EllipticCurve(GF(10007), [1,2,3,4,5])
            sage: E.cardinality(algorithm='pari')
            10076
            sage: E.cardinality(algorithm='foobar')
            10076
        """
        if extension_degree>1:
            # A recursive call to cardinality() with
            # extension_degree=1, which will cache the cardinality, is
            # made by the call to frobenius_order() here:
            R=self.frobenius_order()
            if R.degree()==1:
                return (self.frobenius()**extension_degree-1)**2
            else:
                return (self.frobenius()**extension_degree-1).norm()

        # Now extension_degree==1
        try:
            return self._order
        except AttributeError:
            pass

        k = self.base_ring()
        q = k.cardinality()

        if q < 50:
            return self.cardinality_exhaustive()

        # use special code for j=0, 1728 (for any field)
        j = self.j_invariant()
        if j==k(0):
            return self._cardinality_with_j_invariant_0()
        if j==k(1728):
            return self._cardinality_with_j_invariant_1728()

        N = 0
        p = k.characteristic()
        d = k.degree()

        # Over prime fields, we have a variety of algorithms to choose from:

        if d == 1:
            if algorithm in ('heuristic', 'sea'):  # for backwards compatibility
                algorithm = 'pari'
            if algorithm == 'pari':
                N = self.cardinality_pari()
            elif algorithm == 'bsgs':
                N = self.cardinality_bsgs()
            elif algorithm == 'all':
                N1 = self.cardinality_pari()
                N2 = self.cardinality_bsgs()
                if N1 == N2:
                    N = N1
                else:
                    raise RuntimeError("BUG! Cardinality with pari=%s but with bsgs=%s"%(N1, N2))
            else:
                raise ValueError("Algorithm is not known")
            self._order = Integer(N)
            return self._order

        # now k is not a prime field and j is not 0, 1728

        # we count points on a standard curve with the same
        # j-invariant, defined over the field it generates, then lift
        # to the curve's own field, and finally allow for twists

        # Since j is not 0, 1728 the only twists are quadratic

        j_pol=j.minimal_polynomial()
        j_deg=j_pol.degree()

        # if not possible to work over a smaller field:
        if d==j_deg:
            self._order = self.cardinality_bsgs()
            return self._order

        # Let jkj be the j-invariant as element of the smallest finite
        # field over which j is defined.
        if j_deg == 1:
            # j_pol is of the form X - j
            jkj = -j_pol[0]
        else:
            jkj = GF(p**j_deg, name='a', modulus=j_pol).gen()

        # recursive call which will do all the real work:
        Ej = EllipticCurve_from_j(jkj)
        N=Ej.cardinality(extension_degree=d//j_deg)

        # if curve ia a (quadratic) twist of the "standard" one:
        if not self.is_isomorphic(EllipticCurve_from_j(j)):
            N = 2*(q+1) - N

        self._order = N
        return self._order

    order = cardinality # alias

    def frobenius_polynomial(self):
        r"""
        Return the characteristic polynomial of Frobenius.

        The Frobenius endomorphism of the elliptic curve has quadratic
        characteristic polynomial. In most cases this is irreducible and
        defines an imaginary quadratic order; for some supersingular
        curves, Frobenius is an integer a and the polynomial is
        `(x-a)^2`.

        .. note::

           This computes the curve cardinality, which may be
           time-consuming.

        EXAMPLES::

            sage: E=EllipticCurve(GF(11),[3,3])
            sage: E.frobenius_polynomial()
            x^2 - 4*x + 11

        For some supersingular curves, Frobenius is in Z and the polynomial
        is a square::

            sage: E=EllipticCurve(GF(25,'a'),[0,0,0,0,1])
            sage: E.frobenius_polynomial().factor()
            (x + 5)^2
        """
        x=polygen(ZZ)
        return x**2-self.trace_of_frobenius()*x+self.base_field().cardinality()

    def frobenius_order(self):
        r"""
        Return the quadratic order Z[phi] where phi is the Frobenius
        endomorphism of the elliptic curve

        .. note::

           This computes the curve cardinality, which may be
           time-consuming.

        EXAMPLES::

            sage: E=EllipticCurve(GF(11),[3,3])
            sage: E.frobenius_order()
            Order in Number Field in phi with defining polynomial x^2 - 4*x + 11

        For some supersingular curves, Frobenius is in Z and the Frobenius
        order is Z::

            sage: E=EllipticCurve(GF(25,'a'),[0,0,0,0,1])
            sage: R=E.frobenius_order()
            sage: R
            Order in Number Field in phi with defining polynomial x + 5
            sage: R.degree()
            1
        """
        f = self.frobenius_polynomial().factor()[0][0]
        return ZZ.extension(f,names='phi')

    def frobenius(self):
        r"""
        Return the frobenius of self as an element of a quadratic order

        .. note::

           This computes the curve cardinality, which may be
           time-consuming.

        Frobenius is only determined up to conjugacy.

        EXAMPLES::

            sage: E=EllipticCurve(GF(11),[3,3])
            sage: E.frobenius()
            phi
            sage: E.frobenius().minpoly()
            x^2 - 4*x + 11

        For some supersingular curves, Frobenius is in Z::

            sage: E=EllipticCurve(GF(25,'a'),[0,0,0,0,1])
            sage: E.frobenius()
            -5
        """
        R = self.frobenius_order()
        if R.degree()==1:
            return self.frobenius_polynomial().roots(multiplicities=False)[0]
        else:
            return R.gen(1)

    def cardinality_exhaustive(self):
        r"""
        Return the cardinality of self over the base field. Simply adds up
        the number of points with each x-coordinate: only used for small
        field sizes!

        EXAMPLES::

            sage: p=next_prime(10^3)
            sage: E=EllipticCurve(GF(p),[3,4])
            sage: E.cardinality_exhaustive()
            1020
            sage: E=EllipticCurve(GF(3^4,'a'),[1,1])
            sage: E.cardinality_exhaustive()
            64
        """
        self._order = Integer(1+sum([len(self.lift_x(x,all=True)) for x in self.base_field()]))
        return self._order

    def cardinality_pari(self):
        r"""
        Return the cardinality of self over the (prime) base field using PARI.

        The result is not cached.

        EXAMPLES::

            sage: p=next_prime(10^3)
            sage: E=EllipticCurve(GF(p),[3,4])
            sage: E.cardinality_pari()
            1020
            sage: K=GF(next_prime(10^6))
            sage: E=EllipticCurve(K,[1,0,0,1,1])
            sage: E.cardinality_pari()
            999945

        TESTS::

            sage: K.<a>=GF(3^20)
            sage: E=EllipticCurve(K,[1,0,0,1,a])
            sage: E.cardinality_pari()
            Traceback (most recent call last):
            ...
            ValueError: cardinality_pari() only works over prime fields.
            sage: E.cardinality()
            3486794310

        """
        k = self.base_ring()
        p = k.characteristic()
        if k.degree()==1:
            return ZZ(p + 1 - int(self._pari_().ellap(p)))
        else:
            raise ValueError("cardinality_pari() only works over prime fields.")

    def cardinality_bsgs(self, verbose=False):
        r"""
        Return the cardinality of self over the base field. Will be called
        by user function cardinality only when necessary, i.e. when the
        j_invariant is not in the prime field.

        ALGORITHM: A variant of "Mestre's trick" extended to all finite
        fields by Cremona and Sutherland, 2008.

        .. note::

           1. The Mestre-Schoof-Cremona-Sutherland algorithm may fail for
              a small finite number of curves over `F_q` for `q` at most 49, so
              for `q<50` we use an exhaustive count.

           2. Quadratic twists are not implemented in characteristic 2
              when `j=0 (=1728)`; but this case is treated separately.

       EXAMPLES::

            sage: p=next_prime(10^3)
            sage: E=EllipticCurve(GF(p),[3,4])
            sage: E.cardinality_bsgs()
            1020
            sage: E=EllipticCurve(GF(3^4,'a'),[1,1])
            sage: E.cardinality_bsgs()
            64
            sage: F.<a>=GF(101^3,'a')
            sage: E=EllipticCurve([2*a^2 + 48*a + 27, 89*a^2 + 76*a + 24])
            sage: E.cardinality_bsgs()
            1031352
        """
        E1 = self
        k = self.base_field()
        q = k.order()
        if q<50:
            if verbose:
                print "q=",q,"< 50 so using exhaustive count"
            return self.cardinality_exhaustive()

        # Construct the quadratic twist:
        E2 = E1.quadratic_twist()
        if verbose:
            print "Quadratic twist is ",E2.ainvs()

        bounds = Hasse_bounds(q)
        lower, upper = bounds
        B = upper-q-1 # = floor(2*sqrt(q))
        a = ZZ(0)
        N1 = N2 = M = ZZ(1)
        kmin = -B
        kmax = B
        q1 = q+1
        # Throughout, we have #E=q+1-t where |t|<=B and t=a+k*M = a
        # (mod M) where kmin <= k <= kmax.

        # M is the lcm of the orders of all the points found on E1 and
        # E2, which will eventually exceed 2*B, at which point
        # kmin=kmax.

        if q > 2**10:
            N1 = ZZ(2)**sum([e for P,e in E1._p_primary_torsion_basis(2)])
            N2 = ZZ(2)**sum([e for P,e in E2._p_primary_torsion_basis(2)])
            if q > 2**20:
                N1 *= ZZ(3)**sum([e for P,e in E1._p_primary_torsion_basis(3)])
                N2 *= ZZ(3)**sum([e for P,e in E2._p_primary_torsion_basis(3)])
                if q > 2**40:
                    N1 *= ZZ(5)**sum([e for P,e in E1._p_primary_torsion_basis(5)])
                    N2 *= ZZ(5)**sum([e for P,e in E2._p_primary_torsion_basis(5)])
            # We now know that t=q+1 (mod N1) and t=-(q+1) (mod N2)
            a = q1
            M = N1
            g,u,v = M.xgcd(N2) # g==u*M+v*N2
            if N2>g:
                a = (a*v*N2-q1*u*M)//g
                M *= (N2//g) # = lcm(M,N2)
                a = a%M
                if verbose:
                    print "(a,M)=",(a,M)
                kmin = ((-B-a)/M).ceil()
                kmax = ((B-a)/M).floor()
                if kmin==kmax:
                    self._order = q1-a-kmin*M
                    if verbose: print "no random points were needed"
                    return self._order
            if verbose: print "(2,3,5)-torsion subgroup gives M=",M

        # N1, N2 are divisors of the orders of E1, E2 separately,
        # which are used to speed up the computation of the orders of
        # points on each curve.  For large q it is worth initializing
        # these with the full order of the (2,3,5)-torsion which are
        # often non-trivial.

        while kmax!=kmin:
            # Get a random point on E1 and find its order, using the
            # Hasse bounds and the fact that we know that the group
            # order is a multiple of N1:
            n = generic.order_from_bounds(E1.random_point(),bounds,N1,operation='+')
            if verbose: print "New point on E has order ",n
            # update N1 and M
            N1 = N1.lcm(n)
            g,u,v = M.xgcd(n) # g==u*M+v*n
            if n>g:
                # update congruence a (mod M) with q+1 (mod n)
                a = (a*v*n+q1*u*M)//g
                M *= (n//g) # = lcm(M,n)
                a = a%M
                if verbose: print "(a,M)=",(a,M)
                kmin = ((-B-a)/M).ceil()
                kmax = ((B-a)/M).floor()
                if kmin==kmax:
                    self._order = q1-a-kmin*M
                    return self._order
                if verbose: print "number of possibilities is now ",kmax-kmin+1

            # Get a random point on E2 and find its order, using the
            # Hasse bounds and the fact that we know that the group
            # order is a multiple of N2:
            n = generic.order_from_bounds(E2.random_point(),bounds,N2,operation='+')
            if verbose:  print "New point on E' has order ",n
            # update N2 and M
            N2 = N2.lcm(n)
            g,u,v = M.xgcd(n) # g==u*M+v*n
            if n>g:
                # update congruence a (mod M) with -(q+1) (mod n)
                a = (a*v*n-q1*u*M)//g
                M *= (n//g) # = lcm(M,n)
                a = a%M
                if verbose: print "(a,M)=",(a,M)
                kmin = ((-B-a)/M).ceil()
                kmax = ((B-a)/M).floor()
                if kmin==kmax:
                    self._order = q1-a-kmin*M
                    return self._order
                if verbose: print "number of possibilities is now ",kmax-kmin+1

    def gens(self):
        """
        Returns a tuple of length up to 2 of points which generate the
        abelian group of points on this elliptic curve. See
        abelian_group() for limitations.

        The algorithm uses random points on the curve, and hence the
        generators are likely to differ from one run to another; but they
        are cached so will be consistent in any one run of Sage.

        AUTHORS:

        - John Cremona

        EXAMPLES::

            sage: E=EllipticCurve(GF(11),[2,5])
            sage: E.gens()                           # random output
            ((0 : 7 : 1),)
            sage: EllipticCurve(GF(41),[2,5]).gens() # random output
            ((21 : 1 : 1), (8 : 0 : 1))
            sage: F.<a>=GF(3^6,'a')
            sage: E=EllipticCurve([a,a+1])
            sage: pts=E.gens()
            sage: len(pts)
            1
            sage: pts[0].order()==E.cardinality()
            True
        """
        try:
            G =  self.abelian_group()
            return [x.element() for x in G.gens()]
        except AttributeError:
            pass

    def __iter__(self):
        """
        Return an iterator through the points of this elliptic curve.

        EXAMPLES::

            sage: E = EllipticCurve(GF(11), [1,2])
            sage: for P in E:  print P, P.order()
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
        Return the n'th point in self's __points list. This enables users
        to iterate over the curve's point set.

        EXAMPLE::

            sage: E=EllipticCurve(GF(97),[2,3])
            sage: S=E.points()
            sage: E[10]
            (10 : 76 : 1)
            sage: E[15]
            (17 : 10 : 1)
            sage: for P in E: print P.order()
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

    def abelian_group(self, debug=False):
        r"""
        Returns the abelian group structure of the group of points on this
        elliptic curve.

        .. warning::

           The algorithm is definitely *not* intended for use with
           *large* finite fields! The factorization of the orders of
           elements must be feasible. Also, baby-step-giant-step
           methods are used which have space and time requirements
           which are `O(\sqrt{q})`.

        Also, the algorithm uses random points on the curve and hence the
        generators are likely to differ from one run to another; but the
        group is cached so the generators will not change in any one run of
        Sage.

        INPUT:


        -  ``debug`` - (default: False): if True, print
           debugging messages


        OUTPUT:

        - an abelian group

        - tuple of images of each of the generators of the abelian
          group as points on this curve

        AUTHORS:

        - John Cremona

        EXAMPLES::

            sage: E=EllipticCurve(GF(11),[2,5])
            sage: E.abelian_group()
            Additive abelian group isomorphic to Z/10 embedded in Abelian group of points on Elliptic Curve defined by y^2 = x^3 + 2*x + 5 over Finite Field of size 11

        ::

            sage: E=EllipticCurve(GF(41),[2,5])
            sage: E.abelian_group()
            Additive abelian group isomorphic to Z/22 + Z/2 ...

        ::

            sage: F.<a>=GF(3^6,'a')
            sage: E=EllipticCurve([a^4 + a^3 + 2*a^2 + 2*a, 2*a^5 + 2*a^3 + 2*a^2 + 1])
            sage: E.abelian_group()
            Additive abelian group isomorphic to Z/26 + Z/26 ...

        ::

            sage: F.<a>=GF(101^3,'a')
            sage: E=EllipticCurve([2*a^2 + 48*a + 27, 89*a^2 + 76*a + 24])
            sage: E.abelian_group()
            Additive abelian group isomorphic to Z/1031352 ...

        The group can be trivial::

            sage: E=EllipticCurve(GF(2),[0,0,1,1,1])
            sage: E.abelian_group()
            Trivial group embedded in Abelian group of points on Elliptic Curve defined by y^2 + y = x^3 + x + 1 over Finite Field of size 2

        Of course, there are plenty of points if we extend the field::

            sage: E.cardinality(extension_degree=100)
            1267650600228231653296516890625

        This tests the patch for trac #3111, using 10 primes randomly
        selected::

            sage: E = EllipticCurve('389a')
            sage: for p in [5927, 2297, 1571, 1709, 3851, 127, 3253, 5783, 3499, 4817]:
            ...       G = E.change_ring(GF(p)).abelian_group()
            sage: for p in prime_range(10000):  # long time (19s on sage.math, 2011)
            ...       if p != 389:
            ...           G = E.change_ring(GF(p)).abelian_group()

        This tests that the bug reported in trac #3926 has been fixed::

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
        if not debug:
            # if we're in debug mode, always recalculate
            try:
                return self.__abelian_group
            except AttributeError:
                pass

        k = self.base_field()
        q = k.order()
        p = k.characteristic()
        d = k.degree()
        j = self.j_invariant()
        if d>1:
            d = j.minimal_polynomial().degree()


        # Before computing the group structure we compute the
        # cardinality.  While this is not strictly necessary, it makes
        # the code simpler and also makes computation of orders of
        # points faster.

        # j=0,1728

        if j==k(0):
            N = self._cardinality_with_j_invariant_0()
        if j==k(1728):
            N = self._cardinality_with_j_invariant_1728()

        bounds = Hasse_bounds(q)
        lower, upper = bounds
        if debug:
            print "Lower and upper bounds on group order: [",lower,",",upper,"]"

        try:
            N=self._order
            if debug:
                print "Group order already known to be ",N
        except Exception:
            if (q<50):
                if debug:
                    print "Computing group order naively"
                N=self.cardinality_exhaustive()
            elif d==1:
                if debug:
                    print "Computing group order using PARI"
                N=self.cardinality()
            else:
                if debug:
                    print "Computing group order using bsgs"
                N=self.cardinality_bsgs()
            if debug:
                print "... group order = ",N

        self._order=N
        plist = N.prime_factors()
        P1=self(0)
        P2=self(0)
        n1= Integer(1)
        n2= Integer(1)
        P1._order=n1
        P2._order=n2
        npts = 0

        # At all stages the current subgroup is generated by P1, P2 with
        # orders n1,n2 which are disjoint.  We stop when n1*n2 == N

        if debug:
            "About to start generating random points"

        while n1*n2 != N:
            if debug:
                "Getting a new random point"
            Q = self.random_point()
            while Q.is_zero(): Q = self.random_point()
            npts += 1
            if debug:
                print "Q = ",Q,":",
                print " Order(Q) = ", Q.order()

            Q1=n1*Q;

            if Q1.is_zero() and npts>=10: # then P1,n1 will not change but we may increase n2
                if debug: print "Case 2: n2 may increase"
                n1a = 1; n1b = n1
                P1a = P1
                n1a = n1.prime_to_m_part(N//n1)
                n1b = n1//n1a
                Q = n1a*Q       # has order | n1b
                P1a = n1a*P1    # has order = n1b
                if debug: print "n1a=",n1a
                a = None
                for m in n1b.divisors():
                    try:
                        a = generic.bsgs(m*P1a,m*Q,(0,(n1b//m)-1),operation='+')
                        break
                    except ValueError:
                        pass
                assert a is not None
                a *= (m*n1a)
                if debug: print "linear relation gives m=",m,", a=",a
                if debug: assert m*Q==a*P1
                if m>1: # else Q is in <P1>
                    Q=Q-(a//m)*P1; # has order m and is disjoint from P1
                    if debug: assert Q.order()==m
                    Q._order=m
                    if n2==1: # this is our first nontrivial P2
                        P2=Q
                        n2=m
                        if debug:
                            print "Adding second generator ",P2," of order ",n2
                            print "Subgroup order now ",n1*n2,"=",n1,"*",n2
                    else:     # we must merge P2 and Q:
                        oldn2=n2 # holds old value
                        P2,n2=generic.merge_points((P2,n2),(Q,m),operation='+', check=debug)
                        if debug: assert P2.order()==n2
                        P2._order=n2
                        if debug:
                            if n2>oldn2:
                                print "Replacing second generator by ",P2,
                                print " of order ",n2, "  gaining index ",n2//oldn2
                                print "Subgroup order now ",n1*n2,"=",n1,"*",n2
            elif not Q1.is_zero(): # Q1 nonzero: n1 will increase
                if debug:  print "Case 1: n1 may increase"
                oldn1=n1
                if n2>1:
                    P3=(n1//n2)*P1  # so P2,P3 are a basis for n2-torsion
                    if debug: assert P3.order()==n2
                    P3._order=n2
                    if debug: print "storing generator ",P3," of ",n2,"-torsion"
                m = generic.order_from_multiple(Q,N,plist,operation='+', check=debug)
                P1,n1=generic.merge_points((P1,n1),(Q,m), check=debug)
                if debug: assert P1.order()==n1
                P1._order=n1
                if debug:
                    print "Replacing first  generator by ",P1," of order ",
                    print n1,", gaining index ",n1//oldn1
                    print "Subgroup order now ",n1*n2,"=",n1,"*",n2

                # Now replace P2 by a point of order n2 s.t. it and
                # (n1//n2)*P1 are still a basis for n2-torsion:
                if n2>1:
                    a,m = generic.linear_relation(P1,P3,operation='+')
                    if debug: print "linear relation gives m=",m,", a=",a
                    P3=P3-(a//m)*P1
                    if debug: assert P3.order()==m
                    P3._order=m
                    if debug: print "First  P2 component =",P3
                    if m==n2:
                        P2=P3
                    else:
                        a,m = generic.linear_relation(P1,P2,operation='+')
                        if debug: print "linear relation gives m=",m,", a=",a
                        P2=P2-(a//m)*P1;
                        if debug: assert P2.order()==m
                        P2._order=m
                        if debug: print "Second  P2 component =",P2
                        P2,n2=generic.merge_points((P2,n2),(P3,m), check=debug)
                        if debug: assert P2.order()==n2
                        P2._order=n2
                        if debug: print "Combined P2 component =",P2

            if debug:
                if P1.order()!=n1:
                    print "Generator P1 = ",P1," has order ",P1.order(),
                    print " and not ",n1
                    raise ValueError
                if P2.order()!=n2:
                    print "Generator P2 = ",P2," has order ",P2.order()
                    print " and not ",n2
                    raise ValueError
                if n2>1:
                    if generic.linear_relation(P1,P2,operation='+')[1]!=n2:
                        print "Generators not independent!"
                        raise ValueError
                print "Generators: P1 = ",P1," of order ",n1,
                print ", P2 = ",P2," of order ",n2
                print "Subgroup order is now ",n1*n2,"=",n1,"*",n2

        # Finished: record group order, structure and generators

        from sage.groups.additive_abelian.additive_abelian_wrapper import AdditiveAbelianGroupWrapper
        self._order = n1*n2
        if n1==1:
            self.__abelian_group = AdditiveAbelianGroupWrapper(self.point_homset(), [], [])
        else:
            if n2==1:
                self.__abelian_group = AdditiveAbelianGroupWrapper(self.point_homset(), [P1], [n1])
            else:
                self.__abelian_group = AdditiveAbelianGroupWrapper(self.point_homset(), [P1, P2], [n1, n2])
        return self.__abelian_group

    def is_isogenous(self, other, field=None, proof=True):
        """
        Returns whether or not self is isogenous to other

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
        from ell_generic import is_EllipticCurve
        if not is_EllipticCurve(other):
            raise ValueError("Second argument is not an Elliptic Curve.")
        if self.is_isomorphic(other):
            return True
        elif self.base_field().characteristic() != other.base_field().characteristic():
            raise ValueError("The base fields must have the same characteristic.")
        elif field is None:
            if self.base_field().degree() == other.base_field().degree():
                if self.cardinality() == other.cardinality():
                    return True
                else:
                    return False
            elif self.base_field().degree() == gcd(self.base_field().degree(),other.base_field().degree()):
                if self.cardinality(extension_degree=other.base_field().degree()//self.base_field().degree()) == other.cardinality():
                    return True
                else:
                    return False
            elif other.base_field().degree() == gcd(self.base_field().degree(),other.base_field().degree()):
                if other.cardinality(extension_degree=self.base_field().degree()//other.base_field().degree()) == self.cardinality():
                    return True
                else:
                    return False
            else:
                raise ValueError("Curves have different base fields: use the field parameter.")
        else:
            if not lcm(self.base_field().degree(), other.base_field().degree()).divides(field.degree()):
                raise ValueError("Field must be an extension of the base fields of both curves")
            else:
                if \
self.cardinality(extension_degree=field.degree()//self.base_field().degree())\
 == other.cardinality(extension_degree=field.degree()//other.base_field().degree()):
                      return True
                else:
                      return False

    def is_supersingular(self, proof=True):
        r"""
        Return True if this elliptic curve is supersingular, else False.

        INPUT:

        - ``proof`` (boolean, default True) -- If True, returns a
          proved result.  If False, then a return value of False is
          certain but a return value of True may be based on a
          probabilistic test.  See the documentaion of the function
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
          probabilistic test.  See the documentaion of the function
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
            sage: E.order() * G  # This takes practically no time.
            (0 : 1 : 0)
            sage: proof.arithmetic(prev_proof_state) # restore state

        It is an error to pass a value which is not an integer in the
        Hasse-Weil range::

            sage: E = EllipticCurve(GF(7), [0, 1]) # This curve has order 6
            sage: E.set_order("hi")
            Traceback (most recent call last):
            ...
            ValueError: Value hi illegal (not an integer in the Hasse range)
            sage: E.set_order(3.14159)
            Traceback (most recent call last):
            ...
            ValueError: Value 3.14159000000000 illegal (not an integer in the Hasse range)
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
            WARNING: No checking done in set_order
            sage: E.order()
            4

        The value of num_checks must be an integer. Negative values
        are interpreted as zero, which means don't do any checking::

            sage: E = EllipticCurve(GF(7), [0, 1]) # This curve has order 6
            sage: E.set_order(4, num_checks=-12)
            WARNING: No checking done in set_order
            sage: E.order()
            4

        NOTES:

        The implementation is based on the fact that orders of elliptic curves
        are cached in the (pseudo-private) _order slot.

        AUTHORS:

         - Mariah Lenox (2011-02-16)
        """
        # Is value in the Hasse range?
        q = self.base_field().order()
        a,b = Hasse_bounds(q,1)
        #a = q + 1 - 2*q.isqrt()
        #b = q + 1 + 2*q.isqrt()
        if not value in ZZ:
            raise ValueError('Value %s illegal (not an integer in the Hasse range)'%value)
        if not a <= value <= b:
            raise ValueError('Value %s illegal (not an integer in the Hasse range)'%value)
        # Is value*random == identity?
        for i in range(num_checks):
            G = self.random_point()
            if value * G != self(0):
                raise ValueError('Value %s illegal (multiple of random point not the identity)'%value)
        if(num_checks <= 0):
            print 'WARNING: No checking done in set_order'
        self._order = value

def supersingular_j_polynomial(p):
    """
    Return a polynomial whose roots are the supersingular `j`-invariants
    in characteristic `p`, other than 0, 1728.

    INPUT:

    - `p` (integer) -- a prime number.

    ALGORITHM:

    First compute H(X) whose roots are the Legendre
    `\lambda`-invariants of supersingular curves (Silverman V.4.1(b))
    in charactersitic `p`.  Then, using a resultant computation with
    the polynomial relating `\lambda` and `j` (Silverman III.1.7(b)),
    we recover the polynomial (in variable ``j``) whose roots are the
    `j`-invariants.  Factors of `j` and `j-1728` are removed if
    present.

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

        sage: supersingular_j_polynomial(6)
        Traceback (most recent call last):
        ...
        ValueError: p (=6) should be a prime number

    """
    try:
        p = ZZ(p)
    except TypeError:
        raise ValueError("p (=%s) should be a prime number"%p)
    if not p.is_prime():
        raise ValueError("p (=%s) should be a prime number"%p)

    J = polygen(GF(p),'j')
    if p<13:
        return J.parent().one()
    from sage.rings.arith import binomial
    from sage.misc.all import prod
    m=(p-1)//2
    X,T = PolynomialRing(GF(p),2,names=['X','T']).gens()
    H = sum([binomial(m,i)**2 * T**i for i in xrange(m+1)])
    F = T**2 * (T-1)**2 * X - 256*(T**2-T+1)**3
    R = F.resultant(H,T)
    R =  prod([fi for fi,e in R([J,0]).factor()])
    if R(0)==0:
        R = R//J
    if R(1728)==0:
        R = R//(J-1728)
    return R

# For p in [13..300] we have precomputed these polynomials and store
# them (as lists of their coefficients in ZZ) in a dict:

supersingular_j_polynomials = dict()

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
        raise ValueError("%s must be an element of a finite field"%j)

    F = j.parent()
    p = F.characteristic()
    d = F.degree()

    if j.is_zero():
        return p==3 or p%3==2

    if (j-1728).is_zero():
        return p==2 or p%4==3

    # From now on we know that j != 0, 1728

    if p in (2,3,5,7,11):
        return False # since j=0, 1728 are the only s.s. invariants

    # supersingular j-invariants have degree at most 2:

    jpol = j.minimal_polynomial()
    degj = jpol.degree()
    if degj > 2:
        return False

    # if p occurs in the precomputed list, use that:

    try:
        coeffs = supersingular_j_polynomials[p]
        return PolynomialRing(F,'x')(coeffs)(j).is_zero()
    except KeyError:
        pass

    # Over GF(p), supersingular elliptic curves have cardinality
    # exactly p+1, so we check some random points in order to detect
    # non-supersingularity.  Over GF(p^2) (for p at least 5) the
    # cardinality is either (p-1)^2 or (p+1)^2, and the group has
    # exponent p+1 or p-1, so we can do a similar random check: unless
    # (p+1)*P=0 for all the random points, or (p-1)*P=0 for all of
    # them, we can certainly return False.

    # First we replace j by an element of GF(p) or GF(p^2) (since F
    # might be a proper extension of these):

    if degj==1:
        j = -jpol(0) # = j, but in GF(p)
    elif d>2:
        F = GF(p^2,'a')
        j = jpol.roots(F,multiplicities=False)[0] # j, but in GF(p^2)

    E = EllipticCurve(j=j)
    if degj==1:
        for i in range(10):
            P = E.random_element()
            if not ((p+1)*P).is_zero():
                return False
    else:
        n = None # will hold either p+1 or p-1 later
        for i in range(10):
            P = E.random_element()
            # avoid 2-torsion;  we know that a1=a3=0 and #E>4!
            while P[2].is_zero() or P[1].is_zero():
                P = E.random_element()

            if n is None:  # not yet decided between p+1 and p-1
                pP = p*P
                if not pP[0]==P[0]: # i.e. pP is neither P nor -P
                    return False
                if pP[1]==P[1]: # then p*P == P != -P
                    n=p-1
                else:           # then p*P == -P != P
                    n=p+1
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


