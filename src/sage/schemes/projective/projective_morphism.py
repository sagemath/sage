r"""
Morphisms on projective varieties

A morphism of schemes determined by rational functions that define
what the morphism does on points in the ambient projective space.


AUTHORS:

- David Kohel, William Stein

- William Stein (2006-02-11): fixed bug where P(0,0,0) was allowed as
  a projective point.

- Volker Braun (2011-08-08): Renamed classes, more documentation, misc
  cleanups.

- Ben Hutz (2013-03) iteration functionality and new directory structure
  for affine/projective, height functionality

- Brian Stout, Ben Hutz (Nov 2013) - added minimal model functionality

"""

# Historical note: in trac #11599, V.B. renamed
# * _point_morphism_class -> _morphism
# * _homset_class -> _point_homset

#*****************************************************************************
#       Copyright (C) 2011 Volker Braun <vbraun.name@gmail.com>
#       Copyright (C) 2006 David Kohel <kohel@maths.usyd.edu.au>
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.categories.homset        import Hom
from sage.functions.all            import sqrt
from sage.matrix.constructor       import matrix, identity_matrix
from sage.misc.cachefunc           import cached_method
from sage.misc.misc                import subsets
from sage.misc.mrange              import xmrange
from sage.modules.free_module_element import vector
from sage.rings.all                import Integer, moebius
from sage.rings.arith              import gcd, lcm, next_prime, binomial, primes
from sage.rings.complex_field      import ComplexField
from sage.rings.finite_rings.constructor import GF, is_PrimeFiniteField
from sage.rings.finite_rings.integer_mod_ring import Zmod
from sage.rings.fraction_field     import FractionField
from sage.rings.integer_ring       import ZZ
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.quotient_ring      import QuotientRing_generic
from sage.rings.rational_field     import QQ
from sage.rings.real_mpfr          import RealField
from sage.schemes.generic.morphism import SchemeMorphism_polynomial
from sage.symbolic.constants       import e
from copy import copy


class SchemeMorphism_polynomial_projective_space(SchemeMorphism_polynomial):
    """
    A morphism of schemes determined by rational functions that define
    what the morphism does on points in the ambient projective space.

    EXAMPLES::

        sage: R.<x,y> = QQ[]
        sage: P1 = ProjectiveSpace(R)
        sage: H = P1.Hom(P1)
        sage: H([y,2*x])
        Scheme endomorphism of Projective Space of dimension 1 over Rational Field
          Defn: Defined on coordinates by sending (x : y) to
                (y : 2*x)

    An example of a morphism between projective plane curves (see :trac:`10297`)::

        sage: P2.<x,y,z> = ProjectiveSpace(QQ,2)
        sage: f = x^3+y^3+60*z^3
        sage: g = y^2*z-( x^3 - 6400*z^3/3)
        sage: C = Curve(f)
        sage: E = Curve(g)
        sage: xbar,ybar,zbar = C.coordinate_ring().gens()
        sage: H = C.Hom(E)
        sage: H([zbar,xbar-ybar,-(xbar+ybar)/80])
        Scheme morphism:
          From: Projective Curve over Rational Field defined by x^3 + y^3 + 60*z^3
          To:   Projective Curve over Rational Field defined by -x^3 + y^2*z + 6400/3*z^3
          Defn: Defined on coordinates by sending (x : y : z) to
                (z : x - y : -1/80*x - 1/80*y)

    A more complicated example::

        sage: P2.<x,y,z> = ProjectiveSpace(2,QQ)
        sage: P1 = P2.subscheme(x-y)
        sage: H12 = P1.Hom(P2)
        sage: H12([x^2,x*z, z^2])
        Scheme morphism:
          From: Closed subscheme of Projective Space of dimension 2 over Rational Field defined by:
          x - y
          To:   Projective Space of dimension 2 over Rational Field
          Defn: Defined on coordinates by sending (x : y : z) to
              (y^2 : y*z : z^2)

    We illustrate some error checking::

        sage: R.<x,y> = QQ[]
        sage: P1 = ProjectiveSpace(R)
        sage: H = P1.Hom(P1)
        sage: f = H([x-y, x*y])
        Traceback (most recent call last):
        ...
        ValueError: polys (=[x - y, x*y]) must be of the same degree

        sage: H([x-1, x*y+x])
        Traceback (most recent call last):
        ...
        ValueError: polys (=[x - 1, x*y + x]) must be homogeneous

        sage: H([exp(x),exp(y)])
        Traceback (most recent call last):
        ...
        TypeError: polys (=[e^x, e^y]) must be elements of
        Multivariate Polynomial Ring in x, y over Rational Field
    """

    def __init__(self, parent, polys, check=True):
        """
        The Python constructor.

        See :class:`SchemeMorphism_polynomial` for details.

        EXAMPLES::

            sage: P1.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = P1.Hom(P1)
            sage: H([y,2*x])
            Scheme endomorphism of Projective Space of dimension 1 over Rational Field
              Defn: Defined on coordinates by sending (x : y) to
                    (y : 2*x)
        """
        SchemeMorphism_polynomial.__init__(self, parent, polys, check)
        if check:
            # morphisms from projective space are always given by
            # homogeneous polynomials of the same degree
            try:
                d = polys[0].degree()
            except AttributeError:
                polys = [f.lift() for f in polys]
            if not all([f.is_homogeneous() for f in polys]):
                raise  ValueError("polys (=%s) must be homogeneous"%polys)
            degs = [f.degree() for f in polys]
            if not all([d == degs[0] for d in degs[1:]]):
                raise ValueError("polys (=%s) must be of the same degree"%polys)

    def __eq__(self, right):
        """
        Tests the equality of two projective spaces.

        INPUT:

        - ``right`` - a map on projective space

        OUTPUT:

        - Boolean - True if ``self`` and ``right`` define the same projective map. False otherwise.

        EXAMPLES::

            sage: P1.<x,y> = ProjectiveSpace(RR,1)
            sage: P2.<x,y> = ProjectiveSpace(QQ,1)
            sage: P1==P2
            False

            ::

            sage: R.<x,y> = QQ[]
            sage: P1 = ProjectiveSpace(R)
            sage: P2.<x,y> = ProjectiveSpace(QQ,1)
            sage: P1==P2
            True
        """
        if not isinstance(right, SchemeMorphism_polynomial):
            return False
        else:
            n = len(self._polys)
            for i in range(0,n):
                for j in range(i+1,n):
                    if self._polys[i]*right._polys[j] != self._polys[j]*right._polys[i]:
                        return False
        return True

    def __ne__(self, right):
        """
        Tests the inequality of two projective spaces.

        INPUT:

        - ``right`` -- a map on projective space

        OUTPUT:

        - Boolean -- True if ``self`` and ``right`` define different projective maps. False otherwise.

        EXAMPLES::

            sage: P1.<x,y> = ProjectiveSpace(RR,1)
            sage: P2.<x,y> = ProjectiveSpace(QQ,1)
            sage: P1!=P2
            True

            ::

            sage: R.<x,y> = QQ[]
            sage: P1 = ProjectiveSpace(R)
            sage: P2.<x,y> = ProjectiveSpace(QQ,1)
            sage: P1!=P2
            False
        """
        if not isinstance(right, SchemeMorphism_polynomial):
            return True
        else:
            n=len(self._polys)
            for i in range(0,n):
                for j in range(i+1,n):
                    if self._polys[i]*right._polys[j] != self._polys[j]*right._polys[i]:
                        return True
        return False

    def scale_by(self, t):
        """
        Scales each coordinates by a factor of `t`.

        A ``TypeError`` occurs if the point is not in the coordinate_ring
        of the parent after scaling.

        INPUT:

        - ``t`` -- a ring element

        OUTPUT:

        - None.

        EXAMPLES::

            sage: A.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = Hom(A,A)
            sage: f = H([x^3-2*x*y^2,x^2*y])
            sage: f.scale_by(1/x)
            sage: f
            Scheme endomorphism of Projective Space of dimension 1 over Rational
            Field
              Defn: Defined on coordinates by sending (x : y) to
                    (x^2 - 2*y^2 : x*y)

        ::

            sage: R.<t> = PolynomialRing(QQ)
            sage: P.<x,y> = ProjectiveSpace(R,1)
            sage: H = Hom(P,P)
            sage: f = H([3/5*x^2,6*y^2])
            sage: f.scale_by(5/3*t); f
            Scheme endomorphism of Projective Space of dimension 1 over Univariate
            Polynomial Ring in t over Rational Field
              Defn: Defined on coordinates by sending (x : y) to
                    (t*x^2 : 10*t*y^2)

        ::

            sage: P.<x,y,z> = ProjectiveSpace(GF(7),2)
            sage: X = P.subscheme(x^2-y^2)
            sage: H = Hom(X,X)
            sage: f = H([x^2,y^2,z^2])
            sage: f.scale_by(x-y);f
            Scheme endomorphism of Closed subscheme of Projective Space of dimension
            2 over Finite Field of size 7 defined by:
              x^2 - y^2
              Defn: Defined on coordinates by sending (x : y : z) to
                    (x*y^2 - y^3 : x*y^2 - y^3 : x*z^2 - y*z^2)
        """
        if t == 0:
            raise ValueError("Cannot scale by 0")
        R=self.domain().coordinate_ring()
        if isinstance(R, QuotientRing_generic):
            phi=R.coerce_map_from(self.domain().ambient_space().coordinate_ring())
            for i in range(self.codomain().ambient_space().dimension_relative()+1):
                self._polys[i]=phi(self._polys[i]*t).lift()
        else:
            for i in range(self.codomain().ambient_space().dimension_relative()+1):
                self._polys[i]=R(self._polys[i]*t)

    def normalize_coordinates(self):
        """
        Scales by 1/gcd of the coordinate functions. Also, scales to clear any denominators from the coefficients.
        This is done in place.

        OUTPUT:

        - None.

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = Hom(P,P)
            sage: f = H([5/4*x^3,5*x*y^2])
            sage: f.normalize_coordinates(); f
            Scheme endomorphism of Projective Space of dimension 1 over Rational
            Field
              Defn: Defined on coordinates by sending (x : y) to
                    (x^2 : 4*y^2)

        ::

            sage: P.<x,y,z> = ProjectiveSpace(GF(7),2)
            sage: X = P.subscheme(x^2-y^2)
            sage: H = Hom(X,X)
            sage: f = H([x^3+x*y^2,x*y^2,x*z^2])
            sage: f.normalize_coordinates(); f
            Scheme endomorphism of Closed subscheme of Projective Space of dimension
            2 over Finite Field of size 7 defined by:
              x^2 - y^2
              Defn: Defined on coordinates by sending (x : y : z) to
                    (2*y^2 : y^2 : z^2)

        .. NOTE:: gcd raises an error if the base_ring does not support gcds.
        """
        GCD = gcd(self[0],self[1])
        index=2
        if self[0].lc()>0 or self[1].lc() >0:
            neg=0
        else:
            neg=1
        N=self.codomain().ambient_space().dimension_relative()+1
        while GCD!=1 and index < N:
            if self[index].lc()>0:
                neg=0
            GCD=gcd(GCD,self[index])
            index+=+1

        if GCD != 1:
            R=self.domain().base_ring()
            if neg == 1:
                self.scale_by(R(-1)/GCD)
            else:
                self.scale_by(R(1)/GCD)
        else:
            if neg == 1:
                self.scale_by(-1)

        #clears any denominators from the coefficients
        LCM = lcm([self[i].denominator() for i in range(N)])
        self.scale_by(LCM)

        #scales by 1/gcd of the coefficients.
        GCD = gcd([self[i].content() for i in range(N)])
        if GCD!=1:
            self.scale_by(1/GCD)

    def dynatomic_polynomial(self, period):
        r"""
        For a map `f:\mathbb{P}^1 \to \mathbb{P}^1` this function computes the dynatomic polynomial.

        The dynatomic polynomial is the analog of the cyclotomic
        polynomial and its roots are the points of formal period `n`.

        ALGORITHM:

        For a positive integer `n`, let `[F_n,G_n]` be the coordinates of the `nth` iterate of `f`.
        Then construct

        .. MATH::

            \Phi^{\ast}_n(f)(x,y) = \sum_{d \mid n} (yF_d(x,y) - xG_d(x,y))^{\mu(n/d)}

        where `\mu` is the Moebius function.

        For a pair `[m,n]`, let `f^m = [F_m,G_m]`. Compute

        .. MATH::

            \Phi^{\ast}_{m,n}(f)(x,y) = \Phi^{\ast}_n(f)(F_m,G_m)/\Phi^{\ast}_n(f)(F_{m-1},G_{m-1})

        REFERENCES:

        .. [Hutz] B. Hutz. Efficient determination of rational preperiodic
           points for endomorphisms of projective space.
           :arxiv:`1210.6246`, 2012.

        .. [MoPa] P. Morton and P. Patel. The Galois theory of periodic points
           of polynomial maps. Proc. London Math. Soc., 68 (1994), 225-263.

        INPUT:

        - ``period`` -- a positive integer or a list/tuple `[m,n]` where `m` is the preperiod and `n` is the period

        OUTPUT:

        - If possible, a two variable polynomial in the coordinate ring of ``self``.
          Otherwise a fraction field element of the coordinate ring of ``self``

        .. TODO::

            Do the division when the base ring is p-adic or a function field
            so that the output is a polynomial.

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = Hom(P,P)
            sage: f = H([x^2+y^2,y^2])
            sage: f.dynatomic_polynomial(2)
            x^2 + x*y + 2*y^2

        ::

            sage: P.<x,y> = ProjectiveSpace(ZZ,1)
            sage: H = Hom(P,P)
            sage: f = H([x^2+y^2,x*y])
            sage: f.dynatomic_polynomial(4)
            2*x^12 + 18*x^10*y^2 + 57*x^8*y^4 + 79*x^6*y^6 + 48*x^4*y^8 + 12*x^2*y^10 + y^12

        ::

            sage: P.<x,y> = ProjectiveSpace(CC,1)
            sage: H = Hom(P,P)
            sage: f = H([x^2+y^2,3*x*y])
            sage: f.dynatomic_polynomial(3)
            13.0000000000000*x^6 + 117.000000000000*x^4*y^2 +
            78.0000000000000*x^2*y^4 + y^6

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = Hom(P,P)
            sage: f = H([x^2-10/9*y^2,y^2])
            sage: f.dynatomic_polynomial([2,1])
            x^4*y^2 - 11/9*x^2*y^4 - 80/81*y^6

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = Hom(P,P)
            sage: f = H([x^2-29/16*y^2,y^2])
            sage: f.dynatomic_polynomial([2,3])
            x^12 - 95/8*x^10*y^2 + 13799/256*x^8*y^4 - 119953/1024*x^6*y^6 +
            8198847/65536*x^4*y^8 - 31492431/524288*x^2*y^10 +
            172692729/16777216*y^12

        ::

            sage: P.<x,y> = ProjectiveSpace(ZZ,1)
            sage: H = Hom(P,P)
            sage: f = H([x^2-y^2,y^2])
            sage: f.dynatomic_polynomial([1,2])
            x^2 - x*y

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = Hom(P,P)
            sage: f = H([x^3-y^3,3*x*y^2])
            sage: f.dynatomic_polynomial([0,4])==f.dynatomic_polynomial(4)
            True

        ::

            sage: P.<x,y,z> = ProjectiveSpace(QQ,2)
            sage: H = Hom(P,P)
            sage: f = H([x^2+y^2,x*y,z^2])
            sage: f.dynatomic_polynomial(2)
            Traceback (most recent call last):
            ...
            TypeError: Does not make sense in dimension >1

        ::

            sage: P.<x,y> = ProjectiveSpace(Qp(5),1)
            sage: H = Hom(P,P)
            sage: f = H([x^2+y^2,y^2])
            sage: f.dynatomic_polynomial(2)
            (x^4*y + (2 + O(5^20))*x^2*y^3 - x*y^4 + (2 + O(5^20))*y^5)/(x^2*y -
            x*y^2 + y^3)

        .. TODO:: It would be nice to get this to actually be a polynomial.

        ::

            sage: L.<t> = PolynomialRing(QQ)
            sage: P.<x,y> = ProjectiveSpace(L,1)
            sage: H = Hom(P,P)
            sage: f = H([x^2+t*y^2,y^2])
            sage: f.dynatomic_polynomial(2)
            x^2 + x*y + (t + 1)*y^2

        ::

            sage: K.<c> = PolynomialRing(ZZ)
            sage: P.<x,y> = ProjectiveSpace(K,1)
            sage: H = Hom(P,P)
            sage: f = H([x^2+ c*y^2,y^2])
            sage: f.dynatomic_polynomial([1,2])
            x^2 - x*y + (c + 1)*y^2
       """
        if self.domain() != self.codomain():
            raise TypeError("Must have same domain and codomain to iterate")
        from sage.schemes.projective.projective_space import is_ProjectiveSpace
        if is_ProjectiveSpace(self.domain()) is False:
            raise NotImplementedError("Not implemented for subschemes")
        if self.domain().dimension_relative()>1:
            raise TypeError("Does not make sense in dimension >1")

        if (isinstance(period,(list,tuple)) is False):
            period=[0,period]
        try:
            period[0]=ZZ(period[0])
            period[1]=ZZ(period[1])
        except TypeError:
            raise TypeError("Period and preperiod must be integers")
        if period[1]<=0:
                raise AttributeError("Period must be at least 1")

        if period[0]!=0:
            m=period[0]
            fm=self.nth_iterate_map(m)
            fm1=self.nth_iterate_map(m-1)
            n=period[1]
            PHI=1;
            x=self.domain().gen(0)
            y=self.domain().gen(1)
            F=self._polys
            f=F
            for d in range(1,n+1):
                if n%d == 0:
                    PHI=PHI*((y*F[0]-x*F[1])**moebius(n/d))
                if d !=n: #avoid extra iteration
                    F=[f[0](F[0],F[1]),f[1](F[0],F[1])]
            if m!=0:
                PHI=PHI(fm._polys)/PHI(fm1._polys)
        else:
            PHI=1;
            x=self.domain().gen(0)
            y=self.domain().gen(1)
            F=self._polys
            f=F
            for d in range(1,period[1]+1):
                if period[1]%d == 0:
                    PHI=PHI*((y*F[0]-x*F[1])**moebius(period[1]/d))
                if d !=period[1]: #avoid extra iteration
                    F=[f[0](F[0],F[1]),f[1](F[0],F[1])]
        from sage.rings.polynomial.polynomial_ring import PolynomialRing_general
        from sage.rings.polynomial.multi_polynomial_ring_generic import MPolynomialRing_generic
        if (self.domain().base_ring() == RealField()
            or self.domain().base_ring() == ComplexField()):
            PHI=PHI.numerator()._maxima_().divide(PHI.denominator())[0].sage()
        elif isinstance(self.domain().base_ring(),(PolynomialRing_general,MPolynomialRing_generic)):
            from sage.rings.padics.generic_nodes import is_pAdicField, is_pAdicRing
            from sage.rings.function_field.function_field import is_FunctionField
            BR=self.domain().base_ring().base_ring()
            if is_pAdicField(BR) or is_pAdicRing(BR) or is_FunctionField(BR):
                raise NotImplementedError("Not implemented")
            PHI=PHI.numerator()._maxima_().divide(PHI.denominator())[0].sage()
            #do it again to divide out by denominators of coefficients
            PHI=PHI.numerator()._maxima_().divide(PHI.denominator())[0].sage()
        if PHI.denominator() == 1:
            PHI=self.coordinate_ring()(PHI)
        return(PHI)

    def nth_iterate_map(self, n):
        r"""
        For a map ``self`` this function returns the nth iterate of ``self`` as a
        function on ``self.domain()``

        ALGORITHM:

        Uses a form of successive squaring to reducing computations.


        .. TODO:: This could be improved.

        INPUT:

        - ``n`` -- a positive integer.

        OUTPUT:

        - A map between projective spaces

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = Hom(P,P)
            sage: f = H([x^2+y^2,y^2])
            sage: f.nth_iterate_map(2)
            Scheme endomorphism of Projective Space of dimension 1 over Rational
            Field
              Defn: Defined on coordinates by sending (x : y) to
                    (x^4 + 2*x^2*y^2 + 2*y^4 : y^4)

        ::

            sage: P.<x,y> = ProjectiveSpace(CC,1)
            sage: H = Hom(P,P)
            sage: f = H([x^2-y^2,x*y])
            sage: f.nth_iterate_map(3)
            Scheme endomorphism of Projective Space of dimension 1 over Complex
            Field with 53 bits of precision
              Defn: Defined on coordinates by sending (x : y) to
                    (x^8 + (-7.00000000000000)*x^6*y^2 + 13.0000000000000*x^4*y^4 +
            (-7.00000000000000)*x^2*y^6 + y^8 : x^7*y + (-4.00000000000000)*x^5*y^3
            + 4.00000000000000*x^3*y^5 - x*y^7)

        ::

            sage: P.<x,y,z> = ProjectiveSpace(ZZ,2)
            sage: H = Hom(P,P)
            sage: f = H([x^2-y^2,x*y,z^2+x^2])
            sage: f.nth_iterate_map(2)
            Scheme endomorphism of Projective Space of dimension 2 over Integer Ring
              Defn: Defined on coordinates by sending (x : y : z) to
                    (x^4 - 3*x^2*y^2 + y^4 : x^3*y - x*y^3 : 2*x^4 - 2*x^2*y^2 + y^4
            + 2*x^2*z^2 + z^4)

        ::

            sage: P.<x,y,z> = ProjectiveSpace(QQ,2)
            sage: X = P.subscheme(x*z-y^2)
            sage: H = Hom(X,X)
            sage: f = H([x^2,x*z,z^2])
            sage: f.nth_iterate_map(2)
            Scheme endomorphism of Closed subscheme of Projective Space of dimension
            2 over Rational Field defined by:
              -y^2 + x*z
              Defn: Defined on coordinates by sending (x : y : z) to
                    (x^4 : x^2*z^2 : z^4)
        """
        if self.domain() != self.codomain():
            raise TypeError("Domain and Codomain of function not equal")
        if n <0:
            raise TypeError("Iterate number must be a positive integer")
        N=self.codomain().ambient_space().dimension_relative()+1
        F=list(self._polys)
        D=Integer(n).digits(2)  #need base 2
        Coord_ring=self.codomain().coordinate_ring()
        if isinstance(Coord_ring, QuotientRing_generic):
            PHI=[Coord_ring.gen(i).lift() for i in range(N)]
        else:
            PHI=[Coord_ring.gen(i) for i in range(N)]
        for i in range(len(D)):
            T=tuple([F[j] for j in range(N)])
            for k in range(D[i]):
                PHI=[PHI[j](T) for j in range(N)]
            if i!=len(D)-1: #avoid extra iterate
                F=[F[j](T) for j in range(N)] #'square'
        H=Hom(self.domain(),self.codomain())
        return(H(PHI))

    def nth_iterate(self, P, n, normalize=False):
        r"""
        For a map ``self`` and a point `P` in ``self.domain()``
        this function returns the nth iterate of `P` by ``self``.

        If ``normalize`` is ``True``, then the coordinates are
        automatically normalized.

        .. TODO:: Is there a more efficient way to do this?

        INPUT:

        - ``P`` -- a point in ``self.domain()``

        - ``n`` -- a positive integer.

        - ``normalize`` - Boolean (optional Default: ``False``)

        OUTPUT:

        - A point in ``self.codomain()``

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(ZZ,1)
            sage: H = Hom(P,P)
            sage: f = H([x^2+y^2,2*y^2])
            sage: Q = P(1,1)
            sage: f.nth_iterate(Q,4)
            (32768 : 32768)

        ::

            sage: P.<x,y> = ProjectiveSpace(ZZ,1)
            sage: H = Hom(P,P)
            sage: f = H([x^2+y^2,2*y^2])
            sage: Q = P(1,1)
            sage: f.nth_iterate(Q,4,1)
            (1 : 1)

        Is this the right behavior? ::

            sage: P.<x,y,z> = ProjectiveSpace(QQ,2)
            sage: H = Hom(P,P)
            sage: f = H([x^2,2*y^2,z^2-x^2])
            sage: Q = P(2,7,1)
            sage: f.nth_iterate(Q,2)
            (-16/7 : -2744 : 1)

        ::

            sage: R.<t> = PolynomialRing(QQ)
            sage: P.<x,y,z> = ProjectiveSpace(R,2)
            sage: H = Hom(P,P)
            sage: f = H([x^2+t*y^2,(2-t)*y^2,z^2])
            sage: Q = P(2+t,7,t)
            sage: f.nth_iterate(Q,2)
            (t^4 + 2507*t^3 - 6787*t^2 + 10028*t + 16 : -2401*t^3 + 14406*t^2 -
            28812*t + 19208 : t^4)

        ::

            sage: P.<x,y,z> = ProjectiveSpace(ZZ,2)
            sage: X = P.subscheme(x^2-y^2)
            sage: H = Hom(X,X)
            sage: f = H([x^2,y^2,z^2])
            sage: f.nth_iterate(X(2,2,3),3)
            (256 : 256 : 6561)

        ::

            sage: K.<c> = FunctionField(QQ)
            sage: P.<x,y> = ProjectiveSpace(K,1)
            sage: H = Hom(P,P)
            sage: f = H([x^3-2*x*y^2 - c*y^3,x*y^2])
            sage: f.nth_iterate(P(c,1),2)
            ((c^6 - 9*c^4 + 25*c^2 - c - 21)/(c^2 - 3) : 1)
        """
        return(P.nth_iterate(self,n,normalize))

    def degree(self):
        r"""
        This function returns the degree of ``self``.

        The degree is defined as the degree of the homogeneous
        polynomials that are the coordinates of ``self``.

        OUTPUT:

        - A positive integer

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = Hom(P,P)
            sage: f = H([x^2+y^2,y^2])
            sage: f.degree()
            2

        ::

            sage: P.<x,y,z> = ProjectiveSpace(CC,2)
            sage: H = Hom(P,P)
            sage: f = H([x^3+y^3,y^2*z,z*x*y])
            sage: f.degree()
            3

        ::

            sage: R.<t> = PolynomialRing(QQ)
            sage: P.<x,y,z> = ProjectiveSpace(R,2)
            sage: H = Hom(P,P)
            sage: f = H([x^2+t*y^2,(2-t)*y^2,z^2])
            sage: f.degree()
            2

        ::

            sage: P.<x,y,z> = ProjectiveSpace(ZZ,2)
            sage: X = P.subscheme(x^2-y^2)
            sage: H = Hom(X,X)
            sage: f = H([x^2,y^2,z^2])
            sage: f.degree()
            2
        """
        return(self._polys[0].degree())

    def dehomogenize(self, n):
        r"""
        Returns the standard dehomogenization at the nth coordinate `(\frac{self[0]}{self[n]},\frac{self[1]}{self[n]},...)`.

        Note that the new function is defined over the fraction field
        of the base ring of ``self``.

        INPUT:

        - ``n`` -- a nonnegative integer

        OUTPUT:

        - :class:`SchemeMorphism_polynomial_affine_space` (on nth affine patch)

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(ZZ,1)
            sage: H = Hom(P,P)
            sage: f = H([x^2+y^2,y^2])
            sage: f.dehomogenize(0)
            Scheme endomorphism of Affine Space of dimension 1 over Integer Ring
              Defn: Defined on coordinates by sending (x) to
                    (x^2/(x^2 + 1))

        ::

            sage: P.<x,y,z> = ProjectiveSpace(QQ,2)
            sage: H = Hom(P,P)
            sage: f = H([x^2+y^2,y^2-z^2,2*z^2])
            sage: f.dehomogenize(2)
            Scheme endomorphism of Affine Space of dimension 2 over Rational Field
              Defn: Defined on coordinates by sending (x0, x1) to
                    (1/2*x0^2 + 1/2*x1^2, 1/2*x1^2 - 1/2)

        ::

            sage: R.<t> = PolynomialRing(QQ)
            sage: P.<x,y,z> = ProjectiveSpace(FractionField(R),2)
            sage: H = Hom(P,P)
            sage: f = H([x^2+t*y^2,t*y^2-z^2,t*z^2])
            sage: f.dehomogenize(2)
            Scheme endomorphism of Affine Space of dimension 2 over Fraction Field
            of Univariate Polynomial Ring in t over Rational Field
              Defn: Defined on coordinates by sending (x0, x1) to
                    (1/t*x0^2 + x1^2, x1^2 - 1/t)

        ::

            sage: P.<x,y,z> = ProjectiveSpace(ZZ,2)
            sage: X = P.subscheme(x^2-y^2)
            sage: H = Hom(X,X)
            sage: f = H([x^2,y^2,x*z])
            sage: f.dehomogenize(2)
            Scheme endomorphism of Closed subscheme of Affine Space of dimension 2
            over Integer Ring defined by:
              x0^2 - x1^2
              Defn: Defined on coordinates by sending (x0, x1) to
                    (x1^2/x0, x1^2/x0)
        """
        PS=self.domain()
        A=PS.ambient_space()
        if self._polys[n].substitute({A.gen(n):1}) == 0:
            raise ValueError("Can't dehomogenize at 0 coordinate.")
        else:
            Aff=PS.affine_patch(n)
            S=Aff.ambient_space().coordinate_ring()
            R=A.coordinate_ring()
            phi=R.hom([S.gen(j) for j in range(0,n)] + [1] + [S.gen(j) for j in range(n,A.dimension_relative())],S)
            F=[]
            for i in range(0,A.dimension_relative()+1):
                if i !=n:
                    F.append(phi(self._polys[i])/phi(self._polys[n]))
            H=Hom(Aff,Aff)
            return(H(F))

    def orbit(self, P, N, **kwds):
        r"""
        Returns the orbit of `P` by ``self``. If `n` is an integer it returns `[P,self(P),\ldots,self^n(P)]`.
        If `n` is a list or tuple `n=[m,k]` it returns `[self^m(P),\ldots,self^k(P)]`.
        Automatically normalize the points if ``normalize=True``. Perform the checks on point initialize if
        ``check=True``

        INPUT:

        - ``P`` -- a point in ``self.domain()``

        - ``n`` -- a non-negative integer or list or tuple of two non-negative integers

        kwds:

        - ``check`` -- boolean (optional - default: ``True``)

        - ``normalize`` -- boolean (optional - default: ``False``)


        OUTPUT:

        - a list of points in ``self.codomain()``

        EXAMPLES::

            sage: P.<x,y,z> = ProjectiveSpace(ZZ,2)
            sage: H = Hom(P,P)
            sage: f = H([x^2+y^2,y^2-z^2,2*z^2])
            sage: f.orbit(P(1,2,1),3)
            [(1 : 2 : 1), (5 : 3 : 2), (34 : 5 : 8), (1181 : -39 : 128)]

        ::

            sage: P.<x,y,z> = ProjectiveSpace(ZZ,2)
            sage: H = Hom(P,P)
            sage: f = H([x^2+y^2,y^2-z^2,2*z^2])
            sage: f.orbit(P(1,2,1),[2,4])
            [(34 : 5 : 8), (1181 : -39 : 128), (1396282 : -14863 : 32768)]

        ::

            sage: P.<x,y,z> = ProjectiveSpace(ZZ,2)
            sage: X = P.subscheme(x^2-y^2)
            sage: H = Hom(X,X)
            sage: f = H([x^2,y^2,x*z])
            sage: f.orbit(X(2,2,3),3,normalize=True)
            [(2 : 2 : 3), (2 : 2 : 3), (2 : 2 : 3), (2 : 2 : 3)]

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = Hom(P,P)
            sage: f = H([x^2+y^2,y^2])
            sage: f.orbit(P.point([1,2],False),4,check=False)
            [(1 : 2), (5 : 4), (41 : 16), (1937 : 256), (3817505 : 65536)]

        ::

            sage: K.<c> = FunctionField(QQ)
            sage: P.<x,y> = ProjectiveSpace(K,1)
            sage: H = Hom(P,P)
            sage: f = H([x^2+c*y^2,y^2])
            sage: f.orbit(P(0,1),3)
            [(0 : 1), (c : 1), (c^2 + c : 1), (c^4 + 2*c^3 + c^2 + c : 1)]
        """
        return(P.orbit(self,N,**kwds))

    @cached_method
    def is_morphism(self):
        r"""
        returns ``True`` if self is a morphism (no common zero of defining polynomials).

        The map is a morphism if and only if the ideal generated by
        the defining polynomials is the unit ideal.

        OUTPUT:

        - Boolean

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = Hom(P,P)
            sage: f = H([x^2+y^2,y^2])
            sage: f.is_morphism()
            True

        ::

            sage: P.<x,y,z> = ProjectiveSpace(RR,2)
            sage: H = Hom(P,P)
            sage: f = H([x*z-y*z,x^2-y^2,z^2])
            sage: f.is_morphism()
            False

        ::

            sage: R.<t> = PolynomialRing(GF(5))
            sage: P.<x,y,z> = ProjectiveSpace(R,2)
            sage: H = Hom(P,P)
            sage: f = H([x*z-t*y^2,x^2-y^2,t*z^2])
            sage: f.is_morphism()
            True
        """
        from sage.schemes.projective.projective_space import is_ProjectiveSpace
        if is_ProjectiveSpace(self.domain()) is False or is_ProjectiveSpace(self.codomain()) is False:
            raise NotImplementedError
        R=self.coordinate_ring()
        F=self._polys
        if R.base_ring().is_field():
            J=R.ideal(F)
        else:
            S=PolynomialRing(R.base_ring().fraction_field(),R.gens(),R.ngens())
            J=S.ideal([S.coerce(F[i]) for i in range(R.ngens())])
        if J.dimension()>0:
            return False
        else:
            return True

    def resultant(self, normalize=False):
        r"""
        Computes the resultant of the defining polynomials of ``self`` if ``self`` is a map on the projective line.

        If ``normalize`` is ``True``, then first normalize the coordinate
        functions with :meth:`normalize_coordinates`.

        INPUT:

        - ``normalize`` -- Boolean (optional - default: ``False``)

        OUTPUT:

        - an element of ``self.codomain().base_ring()``

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = Hom(P,P)
            sage: f = H([x^2+y^2,6*y^2])
            sage: f.resultant()
            36

        ::

            sage: R.<t> = PolynomialRing(GF(17))
            sage: P.<x,y> = ProjectiveSpace(R,1)
            sage: H = Hom(P,P)
            sage: f = H([t*x^2+t*y^2,6*y^2])
            sage: f.resultant()
            2*t^2
        """
        if self.domain().dimension_relative() > 1:
            raise TypeError("Only for dimension 1, use self.primes_of_bad_reduction() to get bad primes")
        if normalize is True:
            F=copy(self)
            F.normalize_coordinates()
        else:
            F=self
        x=self.domain().gen(0)
        y=self.domain().gen(1)
        d=self.degree()
        f=F[0].substitute({y:1})
        g=F[1].substitute({y:1})
        res=(f.lc()**(d-g.degree())*g.lc()**(d-f.degree())*f._pari_().polresultant(g, x))
        return(self.codomain().base_ring()(res))

    @cached_method
    def primes_of_bad_reduction(self, check=True):
        r"""
        Determines the primes of bad reduction for a map `self: \mathbb{P}^N \to \mathbb{P}^N`
        defined over `\ZZ` or `\QQ`.

        If ``check`` is ``True``, each prime is verified to be of bad reduction.

        ALGORITHM:

        `p` is a prime of bad reduction if and only if the defining
        polynomials of self have a common zero. Or stated another way,
        `p` is a prime of bad reducion if and only if the radical of
        the ideal defined by the defining polynomials of self is not
        `(x_0,x_1,\ldots,x_N)`.  This happens if and only if some
        power of each `x_i` is not in the ideal defined by the
        defining polynomials of self. This last condition is what is
        checked. The lcm of the coefficients of the monomials `x_i` in
        a groebner basis is computed. This may return extra primes.

        INPUT:

        - ``check`` -- Boolean (optional - default: ``True``)

        OUTPUT:

        - a list of integer primes.

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = Hom(P,P)
            sage: f = H([1/3*x^2+1/2*y^2,y^2])
            sage: print f.primes_of_bad_reduction()
            [2, 3]

        ::

            sage: P.<x,y,z,w> = ProjectiveSpace(QQ,3)
            sage: H = Hom(P,P)
            sage: f = H([12*x*z-7*y^2,31*x^2-y^2,26*z^2,3*w^2-z*w])
            sage: f.primes_of_bad_reduction()
            [2, 3, 7, 13, 31]

        This is an example where check=False returns extra primes::

            sage: P.<x,y,z> = ProjectiveSpace(ZZ,2)
            sage: H = Hom(P,P)
            sage: f = H([3*x*y^2 + 7*y^3 - 4*y^2*z + 5*z^3, -5*x^3 + x^2*y + y^3 + 2*x^2*z, -2*x^2*y + x*y^2 + y^3 - 4*y^2*z + x*z^2])
            sage: f.primes_of_bad_reduction(False)
            [2, 5, 37, 2239, 304432717]
            sage: f.primes_of_bad_reduction()
            [5, 37, 2239, 304432717]
        """
        if self.base_ring() != ZZ and self.base_ring() != QQ:
            raise TypeError("Must be ZZ or QQ")
        from sage.schemes.projective.projective_space import is_ProjectiveSpace
        if is_ProjectiveSpace(self.domain()) is False or is_ProjectiveSpace(self.codomain()) is False:
            raise NotImplementedError
        R=self.coordinate_ring()
        F=self._polys
        if R.base_ring().is_field():
            J=R.ideal(F)
        else:
            S=PolynomialRing(R.base_ring().fraction_field(),R.gens(),R.ngens())
            J=S.ideal([S.coerce(F[i]) for i in range(R.ngens())])
        if J.dimension()>0:
            raise TypeError("Not a morphism.")
        #normalize to coefficients in the ring not the fraction field.
        F=[F[i]*lcm([F[j].denominator() for j in range(len(F))]) for i in range(len(F))]

        #move the ideal to the ring of integers
        if R.base_ring().is_field():
            S=PolynomialRing(R.base_ring().ring_of_integers(),R.gens(),R.ngens())
            F=[F[i].change_ring(R.base_ring().ring_of_integers()) for i in range(len(F))]
            J=S.ideal(F)
        else:
            J=R.ideal(F)
        GB=J.groebner_basis()
        badprimes=[]

        #get the primes dividing the coefficients of the monomials x_i^k_i
        for i in range(len(GB)):
            LT=GB[i].lt().degrees()
            power=0
            for j in range(R.ngens()):
                if LT[j]!=0:
                    power+=1
            if power == 1:
                badprimes=badprimes+GB[i].lt().coefficients()[0].support()
        badprimes=list(set(badprimes))
        badprimes.sort()

        #check to return only the truly bad primes
        if check == True:
            index=0
            while index < len(badprimes):  #figure out which primes are really bad primes...
                S=PolynomialRing(GF(badprimes[index]),R.gens(),R.ngens())
                J=S.ideal([S.coerce(F[j]) for j in range(R.ngens())])
                if J.dimension() == 0:
                    badprimes.pop(index)
                else:
                    index+=1
        return(badprimes)

    def conjugate(self, M):
        r"""
        Conjugates ``self`` by ``M``, i.e. `M^{-1} \circ f \circ M`.

        If possible the map will be defined over the same space as
        ``self``. Otherwise, will try to coerce to the base_ring of
        ``M``.

        INPUT:

        - ``M`` -- a square invertible matrix

        OUTPUT:

        - a map from ``self.domain()`` to ``self.codomain()``.

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(ZZ,1)
            sage: H = Hom(P,P)
            sage: f = H([x^2+y^2,y^2])
            sage: f.conjugate(matrix([[1,2],[0,1]]))
            Scheme endomorphism of Projective Space of dimension 1 over Integer Ring
              Defn: Defined on coordinates by sending (x : y) to
                    (x^2 + 4*x*y + 3*y^2 : y^2)

        ::

            sage: R.<x> = PolynomialRing(QQ)
            sage: K.<i> = NumberField(x^2+1)
            sage: P.<x,y> = ProjectiveSpace(ZZ,1)
            sage: H = Hom(P,P)
            sage: f = H([x^3+y^3,y^3])
            sage: f.conjugate(matrix([[i,0],[0,-i]]))
            Scheme endomorphism of Projective Space of dimension 1 over Integer Ring
              Defn: Defined on coordinates by sending (x : y) to
                    (-x^3 + y^3 : -y^3)

        ::

            sage: P.<x,y,z> = ProjectiveSpace(ZZ,2)
            sage: H = Hom(P,P)
            sage: f = H([x^2+y^2,y^2,y*z])
            sage: f.conjugate(matrix([[1,2,3],[0,1,2],[0,0,1]]))
            Scheme endomorphism of Projective Space of dimension 2 over Integer Ring
              Defn: Defined on coordinates by sending (x : y : z) to
                    (x^2 + 4*x*y + 3*y^2 + 6*x*z + 9*y*z + 7*z^2 : y^2 + 2*y*z : y*z + 2*z^2)

        ::

            sage: P.<x,y> = ProjectiveSpace(ZZ,1)
            sage: H = Hom(P,P)
            sage: f = H([x^2+y^2,y^2])
            sage: f.conjugate(matrix([[2,0],[0,1/2]]))
            Scheme endomorphism of Projective Space of dimension 1 over Multivariate
            Polynomial Ring in x, y over Rational Field
              Defn: Defined on coordinates by sending (x : y) to
                    (2*x^2 + 1/8*y^2 : 1/2*y^2)

        ::

            sage: R.<x> = PolynomialRing(QQ)
            sage: K.<i> = NumberField(x^2+1)
            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = Hom(P,P)
            sage: f = H([1/3*x^2+1/2*y^2,y^2])
            sage: f.conjugate(matrix([[i,0],[0,-i]]))
            Scheme endomorphism of Projective Space of dimension 1 over Multivariate
            Polynomial Ring in x, y over Number Field in i with defining polynomial
            x^2 + 1
              Defn: Defined on coordinates by sending (x : y) to
                    ((1/3*i)*x^2 + (1/2*i)*y^2 : (-i)*y^2)
        """

        if M.is_square() == 1 and M.determinant() != 0 and  M.ncols() == self.domain().ambient_space().dimension_relative()+1:
            X=M*vector(self[0].parent().gens())
            F=vector(self._polys)
            F=F(list(X))
            N=M.inverse()
            F=N*F
            R=self.codomain().ambient_space().coordinate_ring()
            try:
                F = [R(f) for f in F]
                PS=self.codomain()
            except TypeError: #no longer defined over same ring
                R=R.change_ring(M.base_ring())
                F = [R(f) for f in F]
                PS=self.codomain().change_ring(R)
            H=Hom(PS,PS)
            return(H(F))
        else:
            raise TypeError("matrix must be invertible and size dimension +1 ")

    def green_function(self, P, v, **kwds):
        r"""
        Evaluates the local Green's function at the place ``v`` for ``P`` with ``N`` terms of the
        series or, in dimension 1, to within a given error bound.

        Use ``v=0`` for the archimedean place. Must be over `\ZZ` or `\QQ`.

        ALGORITHM:

        See Exercise 5.29 and Figure 5.6 of ``The Arithmetic of Dynamics Systems``, Joseph H. Silverman, Springer, GTM 241, 2007.

        INPUT:

        - ``P`` - a projective point

        - ``v`` - non-negative integer. a place, use v=0 for the archimedean place

        kwds:

        - ``N`` - positive integer. number of terms of the series to use

        - ``prec`` - positive integer, float point or p-adic precision, default: 100

        - ``error_bound`` - a positive real number

        OUTPUT:

        - a real number

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = Hom(P,P)
            sage: f = H([x^2+y^2,x*y])
            sage: f.green_function(P.point([5,2],False),0,N=30)
            1.7315451844777407992085512000
            sage: f.green_function(P.point([2,1],False),0,N=30)
            0.86577259223181088325226209926
            sage: f.green_function(P.point([1,1],False),0,N=30)
            0.43288629610862338612700146098
        """
        if self.base_ring() != ZZ and self.base_ring() != QQ:
            raise TypeError("Must be ZZ or QQ")
        return(P.green_function(self,v,**kwds))

    def canonical_height(self, P, **kwds):
        r"""
        Evaluates the canonical height of ``P`` with respect to ``self``. Must be over `\ZZ` or `\QQ`.

        Specify either the number of terms of the series to evaluate
        or, in dimension 1, the error bound required.

        ALGORITHM:

        The sum of the Green's function at the archimedean place and the places of bad reduction.

        INPUT:

        - ``P`` -- a projective point

        kwds:

        - ``badprimes`` - a list of primes of bad reduction

        - ``N`` - positive integer. number of terms of the series to use in the local green functions

        - ``prec`` - positive integer, float point or p-adic precision, default: 100

        - ``error_bound`` - a positive real number

        OUTPUT:

        - a real number

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(ZZ,1)
            sage: H = Hom(P,P)
            sage: f = H([x^2+y^2,2*x*y]);
            sage: f.canonical_height(P.point([5,4]),error_bound=0.001)
            2.1970553519503404898926835324
            sage: f.canonical_height(P.point([2,1]),error_bound=0.001)
            1.0984430632822307984974382955

        Notice that preperiodic points may not be exactly 0::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = Hom(P,P)
            sage: f = H([x^2-29/16*y^2,y^2]);
            sage: f.canonical_height(P.point([1,4]),N=60)
            1.2024186864216154694752186858e-18

        ::

            sage: P.<x,y,z> = ProjectiveSpace(QQ,2)
            sage: X = P.subscheme(x^2-y^2);
            sage: H = Hom(X,X)
            sage: f = H([x^2,y^2,4*z^2]);
            sage: Q = X([4,4,1])
            sage: f.canonical_height(Q,badprimes=[2])
            0.0013538030870311431824555314882
        """
        if self.base_ring() != ZZ and self.base_ring() != QQ:
            raise TypeError("Must be ZZ or QQ")
        return(P.canonical_height(self,**kwds))

    def global_height(self, prec=None):
        r"""
        Returns the maximum of the heights of the coefficients in any of the coordinate functions of ``self``.

        INPUT:

        - ``prec`` -- desired floating point precision (default:
          default RealField precision).

        OUTPUT:

        - a real number

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = Hom(P,P)
            sage: f = H([1/1331*x^2+1/4000*y^2,210*x*y]);
            sage: f.global_height()
            8.29404964010203

        This function does not automatically normalize::

            sage: P.<x,y,z> = ProjectiveSpace(ZZ,2)
            sage: H = Hom(P,P)
            sage: f = H([4*x^2+100*y^2,210*x*y,10000*z^2]);
            sage: f.global_height()
            9.21034037197618
            sage: f.normalize_coordinates()
            sage: f.global_height()
            8.51719319141624

        ::

            sage: R.<z> = PolynomialRing(QQ)
            sage: K.<w> = NumberField(z^2-2)
            sage: O = K.maximal_order()
            sage: P.<x,y> = ProjectiveSpace(O,1)
            sage: H = Hom(P,P)
            sage: f = H([2*x^2 + 3*O(w)*y^2,O(w)*y^2])
            sage: f.global_height()
            1.44518587894808

        .. TODO:: add heights to integer.pyx and remove special case
        """
        if self.domain().base_ring() == ZZ:
            if prec is None:
                R = RealField()
            else:
                R = RealField(prec)
            H=R(0)
            for i in range(self.domain().ambient_space().dimension_relative()+1):
                C=self[i].coefficients()
                h=max([c.abs() for c in C])
                H=max(H,R(h).log())
            return(H)
        H=0
        for i in range(self.domain().ambient_space().dimension_relative()+1):
            C=self[i].coefficients()
            h=max([c.global_height(prec) for c in C])
            H=max(H,h)
        return(H)

    def height_difference_bound(self, prec=None):
        r"""
        Returns an upper bound on the different bewtween the canonical height of a point with
        respect to ``self`` and the height of the point. ``self`` must be a morphism.

        ALGORITHM:

            Uses a Nullstellensatz argument to compute the constant.
            For details: B. Hutz, Efficient determination of rational preperiodic points for endomorphisms of projective
            space, arxiv:1210.6246 (2012).

        INPUT:

        - ``prec`` - positive integer, float point, default: RealField default

        OUTPUT:

        - a real number

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = End(P)
            sage: f = H([x^2+y^2,x*y]);
            sage: f.height_difference_bound()
            1.38629436111989

        This function does not automatically normalize. ::

            sage: P.<x,y,z> = ProjectiveSpace(ZZ,2)
            sage: H = End(P)
            sage: f = H([4*x^2+100*y^2,210*x*y,10000*z^2]);
            sage: f.height_difference_bound()
            11.0020998412042
            sage: f.normalize_coordinates()
            sage: f.height_difference_bound()
            10.7632079329219
        """
        if self.domain().base_ring() != ZZ and self.domain().base_ring() != QQ:
            raise NotImplementedError("Must be over ZZ or QQ")
        if not self.is_endomorphism():
            raise NotImplementedError("Must be an endomorphism of projective space")
        if prec is None:
            R=RealField()
        else:
            R=RealField(prec)
        N=self.domain().dimension_relative()
        d=self.degree()
        D=(N+1)*(d-1) +1
        #compute upper bound
        U = self.global_height(prec) + R(binomial(N+d,d)).log()
        #compute lower bound - from explicit polynomials of Nullstellensatz
        CR=self.domain().coordinate_ring()
        CR=CR.change_ring(QQ) #.lift() only works over fields
        I=CR.ideal(self.defining_polynomials())
        MCP=[]
        for k in range(N+1):
            CoeffPolys=(CR.gen(k)**D).lift(I)
            Res=1
            h=1
            for j in range(len(CoeffPolys)):
                if CoeffPolys[j]!=0:
                    for i in range(len(CoeffPolys[j].coefficients())):
                        Res=lcm(Res,abs(CoeffPolys[j].coefficients()[i].denominator()))
                        h=max(h,abs(CoeffPolys[j].coefficients()[i].numerator()))
            MCP.append([Res,Res*h]) #since we need to clear denominators
        maxh=1
        gcdRes=0
        for k in range(len(MCP)):
            gcdRes=gcd(gcdRes,MCP[k][0])
            maxh=max(maxh,MCP[k][1])
        L=abs(R(gcdRes/((N+1)*binomial(N + D-d,D-d)*maxh)).log())

        C=max(U,L) #height difference dh(P) - L <= h(f(P)) <= dh(P) +U
        return(C/(d-1))

    def multiplier(self,P,n, check=True):
        r"""
        Returns the multiplier of ``self`` at the `QQ`-rational point ``P`` of period ``n``.
        ``self`` must be an endomorphism of projective space

        INPUT:

        - ``P`` - a point on domain of ``self``

        - ``n`` - a positive integer, the period of ``P``

        - ``check`` -- verify that ``P`` has period ``n``, Default:True

        OUTPUT:

        - a square matrix of size ``self.codomain().dimension_relative()`` in the ``base_ring`` of ``self``

        EXAMPLES::

            sage: P.<x,y,z> = ProjectiveSpace(QQ,2)
            sage: H = End(P)
            sage: f = H([x^2,y^2,4*z^2]);
            sage: Q = P.point([4,4,1],False);
            sage: f.multiplier(Q,1)
            [2 0]
            [0 2]

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = End(P)
            sage: f = H([7*x^2 - 28*y^2,24*x*y])
            sage: f.multiplier(P(2,5),4)
            [231361/20736]

        ::

            sage: P.<x,y> = ProjectiveSpace(CC,1)
            sage: H = End(P)
            sage: f = H([x^3 - 25*x*y^2 + 12*y^3,12*y^3])
            sage: f.multiplier(P(1,1),5)
            [0.389017489711935]

        ::

            sage: P.<x,y> = ProjectiveSpace(RR,1)
            sage: H = End(P)
            sage: f = H([x^2-2*y^2,y^2])
            sage: f.multiplier(P(2,1),1)
            [4.00000000000000]

        ::

            sage: P.<x,y> = ProjectiveSpace(Qp(13),1)
            sage: H = End(P)
            sage: f = H([x^2-29/16*y^2,y^2])
            sage: f.multiplier(P(5,4),3)
            [6 + 8*13 + 13^2 + 8*13^3 + 13^4 + 8*13^5 + 13^6 + 8*13^7 + 13^8 +
            8*13^9 + 13^10 + 8*13^11 + 13^12 + 8*13^13 + 13^14 + 8*13^15 + 13^16 +
            8*13^17 + 13^18 + 8*13^19 + O(13^20)]

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = End(P)
            sage: f = H([x^2-y^2,y^2])
            sage: f.multiplier(P(0,1),1)
            Traceback (most recent call last):
            ...
            ValueError: (0 : 1) is not periodic of period 1

        ..  TODO:: would be better to keep the dehomogenizations for reuse
        """
        if not self.is_endomorphism():
            raise NotImplementedError("Must be an endomorphism of projective space")
        if check:
            if self.nth_iterate(P,n)!=P:
                raise ValueError("%s is not periodic of period %s"% (P,n))
        N=self.domain().dimension_relative()
        l=identity_matrix(FractionField(self.codomain().base_ring()),N,N)
        Q=P
        Q.normalize_coordinates()
        index=N
        indexlist=[]
        while Q[index]==0:
            index-=1
        indexlist.append(index)
        for i in range(0,n):
            F=[]
            R=self(Q)
            R.normalize_coordinates()
            index=N
            while R[index]==0:
                index-=1
            indexlist.append(index)
            S=PolynomialRing(FractionField(self.codomain().base_ring()),N,'x')
            CR=self.coordinate_ring()
            map_vars=list(S.gens())
            map_vars.insert(indexlist[i],1)
            phi=CR.hom(map_vars,S)
            #make map between correct affine patches
            for j in range(N+1):
                if j != indexlist[i+1]:
                    F.append(phi(self._polys[j])/phi(self._polys[indexlist[i+1]]))
                J = matrix(FractionField(S),N,N)
            for j1 in range(0,N):
                for j2 in range(0,N):
                    J[j1,j2]=F[j1].derivative(S.gen(j2))
            l=J(tuple(Q.dehomogenize(indexlist[i])))*l #get the correct order for chain rule matrix multiplication
            Q=R
        return l

    def _multipliermod(self,P,n,p,k):
        r"""
        Returns the multiplier of ``self`` at the point ``P`` of period ``n`` modulo `p^k`.
        self must be an endomorphism of projective space defined over `\QQ` or '\ZZ'.
        This function should not be used at the top level as it does not perform input checks.
        It is used primarily for the rational preperiodic and periodic point algorithms.

        INPUT:

        - ``P`` - a point on domain of ``self``

        - ``n`` - a positive integer, the period of ``P``

        - ``p`` - a positive integer

        - ``k`` - a positive integer

        OUTPUT:

        - a square matrix of size ``self.codomain().dimension_relative()`` in `Zmod(p^k)`.

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = End(P)
            sage: f = H([x^2-29/16*y^2,y^2])
            sage: f._multipliermod(P(5,4),3,11,1)
            [3]

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = End(P)
            sage: f = H([x^2-29/16*y^2,y^2])
            sage: f._multipliermod(P(5,4),3,11,2)
            [80]

        .. TODO:: would be better to keep the dehomogenizations for reuse
        """
        N=self.domain().dimension_relative()
        BR=FractionField(self.codomain().base_ring())
        l=identity_matrix(BR,N,N)
        Q=copy(P)
        g=gcd(Q._coords) #we can't use normalize_coordinates since it can cause denominators
        Q.scale_by(1/g)
        index=N
        indexlist=[]
        while Q[index]%p==0:
            index-=1
        indexlist.append(index)
        for i in range(0,n):
            F=[]
            R=self(Q,False)
            g=gcd(R._coords)
            R.scale_by(1/g)
            for index in range(N+1):
                R._coords[index]=R._coords[index]%(p**k)
            index=N
            while R[index]%p==0:
                index-=1
            indexlist.append(index)
            S=PolynomialRing(BR,N,'x')
            CR=self.coordinate_ring()
            map_vars=list(S.gens())
            map_vars.insert(indexlist[i],1)
            phi=CR.hom(map_vars,S)
            for j in range(N+1):
                if j != indexlist[i+1]:
                    F.append(phi(self._polys[j])/phi(self._polys[indexlist[i+1]]))
            J = matrix(FractionField(S),N,N)
            for j1 in range(0,N):
                for j2 in range(0,N):
                    J[j1,j2]=F[j1].derivative(S.gen(j2))
            l=(J(tuple(Q.dehomogenize(indexlist[i])))*l)%(p**k)
            Q=R
        return(l)

    def possible_periods(self,**kwds):
        r"""
        Returns the set of possible periods for rational periodic points of self.
        Must be defined over `\ZZ` or `\QQ`.

        ALGORITHM:
            Calls ``self.possible_periods()`` modulo all primes of good reduction in range ``prime_bound``.
            Returns the intersection of those lists.

        INPUT:

        kwds:

        - ``prime_bound`` - a list or tuple of two positive integers. Or an integer for the upper bound. (optional)
            default: [1,20].

        - ``bad_primes`` - a list or tuple of integer primes, the primes of bad reduction.  (optional)

        OUTPUT:

        - a list of positive integers.

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = End(P)
            sage: f = H([x^2-29/16*y^2,y^2])
            sage: f.possible_periods()
            [1, 3]

        ::

            sage: PS.<x,y> = ProjectiveSpace(1,QQ)
            sage: H = End(PS)
            sage: f = H([5*x^3 - 53*x*y^2 + 24*y^3, 24*y^3])
            sage: f.possible_periods(prime_bound=[1,5])
            Traceback (most recent call last):
            ...
            ValueError: No primes of good reduction in that range
            sage: f.possible_periods(prime_bound=[1,10])
            [1, 4, 12]
            sage: f.possible_periods(prime_bound=[1,20])
            [1, 4]

        ::

            sage: P.<x,y,z> = ProjectiveSpace(ZZ,2)
            sage: H = End(P)
            sage: f = H([2*x^3 - 50*x*z^2 + 24*z^3,5*y^3 - 53*y*z^2 + 24*z^3,24*z^3])
            sage: f.possible_periods(prime_bound=10)
            [1, 2, 6, 20, 42, 60, 140, 420]
            sage: f.possible_periods(prime_bound=20) # long time
            [1, 20]
        """
        if not self.is_endomorphism():
            raise NotImplementedError("Must be an endomorphism of projective space")
        if self.domain().base_ring()!=ZZ and self.domain().base_ring()!=QQ:
            raise NotImplementedError("Must be ZZ or QQ")


        primebound = kwds.pop("prime_bound",[1,20])
        badprimes = kwds.pop("bad_primes",None)

        if (isinstance(primebound,(list,tuple))==False):
            try:
                primebound=[1,ZZ(primebound)]
            except TypeError:
                raise TypeError("Bound on primes must be an integer")
        else:
            try:
                primebound[0]=ZZ(primebound[0])
                primebound[1]=ZZ(primebound[1])
            except TypeError:
                raise TypeError("Prime bounds must be integers")

        if badprimes==None:
            badprimes=self.primes_of_bad_reduction()
        firstgood=0
        for q in primes(primebound[0],primebound[1]+1):
            if q in badprimes: #bad reduction
                t=0
                #print "Bad Reduction at ", q
            elif firstgood==0:
                F=self.change_ring(GF(q))
                periods=set(F.possible_periods())
                firstgood=1
            else:
                F=self.change_ring(GF(q))
                periodsq=set(F.possible_periods())
                periods=periods.intersection(periodsq)
        if firstgood==0:
            raise ValueError("No primes of good reduction in that range")
        else:
            return(sorted(periods))

    def _preperiodic_points_to_cyclegraph(self,preper):
        r"""
        Given the complete set of periodic or preperiodic points returns the
        digraph representing the orbit. If it is not the complete set, this function
        will not fill in the gaps.


        INPUT:

        - ``preper`` - a list or tuple of projective points. The complete set of rational periodic or preperiodic points.

        OUTPUT:

        - a digraph representing the orbit the rational preperiodic points ``preper`` in projective space.

        Examples::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = End(P)
            sage: f = H([x^2-2*y^2,y^2])
            sage: preper = [P(-2, 1), P(1, 0), P(0, 1), P(1, 1), P(2, 1), P(-1, 1)]
            sage: f._preperiodic_points_to_cyclegraph(preper)
            Looped digraph on 6 vertices
        """
        V=[]
        E=[]
        for i in range(0,len(preper)):
            V.append(str(preper[i]))
            Q=self(preper[i])
            Q.normalize_coordinates()
            E.append([str(Q)])
        from sage.graphs.digraph import DiGraph
        g=DiGraph(dict(zip(V,E)), loops=True)
        return(g)

    def is_PGL_minimal(self, prime_list=None):
        r"""
        Checks if ``self`` is a minimal model in its conjugacy class.  See [Bruin-Molnar]
        and [Molnar] for a description of the algorithm.

        INPUT:

        - ``prime_list`` -- list of primes to check minimality, if None, check all places

        OUTPUT:

        - Boolean - True if ``self`` is minimal, False otherwise.

        EXAMPLES::

            sage: PS.<X,Y> = ProjectiveSpace(QQ,1)
            sage: H = End(PS)
            sage: f = H([X^2+3*Y^2,X*Y])
            sage: f.is_PGL_minimal()
            True

        ::

            sage: PS.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = End(PS)
            sage: f = H([6*x^2+12*x*y+7*y^2,12*x*y])
            sage: f.is_PGL_minimal()
            False

        ::

            sage: PS.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = End(PS)
            sage: f = H([6*x^2+12*x*y+7*y^2,y^2])
            sage: f.is_PGL_minimal()
            Traceback (most recent call last):
            ...
            TypeError: Affine minimality is only considered for maps not of the form
            f or 1/f for a polynomial f.
        """
        if self.base_ring()!=QQ and self.base_ring()!=ZZ:
            raise NotImplementedError("Minimal models only implemented over ZZ or QQ")
        if not self.is_morphism():
            raise TypeError("The function is not a morphism")
        if self.degree()==1:
            raise NotImplementedError("Minimality is only for degree 2 or higher")

        from endPN_minimal_model import affine_minimal
        return(affine_minimal(self, False ,prime_list ,True))

    def minimal_model(self, return_transformation = False,prime_list=None):
        r"""
        Given ``self`` a scheme morphism on the projective line over the rationals,
        determine if ``self`` is minimal. In particular, determine
        if ``self`` is affine minimal, which is enough to decide if it is minimal
        or not. See Proposition 2.10 in [Bruin-Molnar].

        REFERENCES:

        .. [Bruin-Molnar] N. Bruin and A. Molnar, Minimal models for rational
           functions in a dynamical setting
           LMS Journal of Computation and Mathematics, Volume 15 (2012), pp 400-417.

        .. [Molnar] A. Molnar, Fractional Linear Minimal Models of Rational Functions,
           M.Sc. Thesis.

        INPUT:

        - ``self`` -- scheme morphism on the projective line defined over `QQ`.

        - ``return_transformation`` -- a boolean value, default value True. This
                                    signals a return of the ``PGL_2`` transformation
                                    to conjugate ``self`` to the calculated minimal
                                    model. default: False

        - ``prime_list`` -- a list of primes, in case one only wants to determine minimality
                   at those specific primes.

        OUTPUT:

        - a scheme morphism on the projective line which is a minimal model of ``self``.

        - a `PGL(2,QQ)` element which conjugates ``self`` to a minimal model

        EXAMPLES::

            sage: PS.<X,Y> = ProjectiveSpace(QQ,1)
            sage: H = End(PS)
            sage: f = H([X^2+3*Y^2,X*Y])
            sage: f.minimal_model(return_transformation=True)
            (
            Scheme endomorphism of Projective Space of dimension 1 over Rational
            Field
              Defn: Defined on coordinates by sending (X : Y) to
                    (X^2 + 3*Y^2 : X*Y)
            ,
            [1 0]
            [0 1]
            )

        ::

            sage: PS.<X,Y> = ProjectiveSpace(QQ,1)
            sage: H = End(PS)
            sage: f = H([7365/2*X^4 + 6282*X^3*Y + 4023*X^2*Y^2 + 1146*X*Y^3 + 245/2*Y^4, -12329/2*X^4 - 10506*X^3*Y - 6723*X^2*Y^2 - 1914*X*Y^3 - 409/2*Y^4])
            sage: f.minimal_model(return_transformation=True)
            (
            Scheme endomorphism of Projective Space of dimension 1 over Rational
            Field
              Defn: Defined on coordinates by sending (X : Y) to
                    (22176*X^4 + 151956*X^3*Y + 390474*X^2*Y^2 + 445956*X*Y^3 +
            190999*Y^4 : -12329*X^4 - 84480*X^3*Y - 217080*X^2*Y^2 - 247920*X*Y^3 -
            106180*Y^4),
            [2 3]
            [0 1]
            )

        ::

            sage: PS.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = End(PS)
            sage: f = H([6*x^2+12*x*y+7*y^2,12*x*y])
            sage: f.minimal_model()
            Scheme endomorphism of Projective Space of dimension 1 over Rational
            Field
              Defn: Defined on coordinates by sending (x : y) to
                    (x^2 + 12*x*y + 42*y^2 : 2*x*y)

        ::

            sage: PS.<x,y> = ProjectiveSpace(ZZ,1)
            sage: H = End(PS)
            sage: f = H([6*x^2+12*x*y+7*y^2,12*x*y + 42*y^2])
            sage: g,M=f.minimal_model(return_transformation=True)
            sage: f.conjugate(M)==g
            True

        ::

            sage: PS.<X,Y> = ProjectiveSpace(QQ,1)
            sage: H = End(PS)
            sage: f = H([X+Y,X-3*Y])
            sage: f.minimal_model()
            Traceback (most recent call last):
            ...
            NotImplementedError: Minimality is only for degree 2 or higher

        ::

            sage: PS.<X,Y> = ProjectiveSpace(QQ,1)
            sage: H = End(PS)
            sage: f = H([X^2-Y^2,X^2+X*Y])
            sage: f.minimal_model()
            Traceback (most recent call last):
            ...
            TypeError: The function is not a morphism

        """
        if self.base_ring()!=QQ and self.base_ring()!=ZZ:
            raise NotImplementedError("Minimal models only implemented over ZZ or QQ")
        if not self.is_morphism():
            raise TypeError("The function is not a morphism")
        if self.degree()==1:
            raise NotImplementedError("Minimality is only for degree 2 or higher")

        from endPN_minimal_model import affine_minimal
        return(affine_minimal(self, return_transformation,prime_list,False))

class SchemeMorphism_polynomial_projective_space_field(SchemeMorphism_polynomial_projective_space):

    def lift_to_rational_periodic(self,points_modp,B=None):
        r"""
        Given a list of points in projective space over `GF(p)`, determine if they lift to
        `\QQ`-rational periodic points. ``self`` must be an endomorphism of projective
        space defined over `\QQ`

        ALGORITHM:
            Use Hensel lifting to find a `p`-adic approximation for that rational point. The accuracy needed
            is determined by the height bound `B`. Then apply the the LLL algorithm to determine if the lift
            corresponds to a rational point.

            If the point is a point of high multiplicity (multiplier 1) then procedure can be very slow.


        INPUT:

        - ``points_modp`` - a list or tuple of pairs containing a point in projective space
          over `GF(p)` and the possible period.

        - ``B`` - a positive integer - the height bound for a rational preperiodic point. (optional)

        OUTPUT:

        - a list of projective points.

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = End(P)
            sage: f = H([x^2 - y^2,y^2])
            sage: f.lift_to_rational_periodic([[P(0,1).change_ring(GF(7)),4]])
            [[(0 : 1), 2]]

        ::

            There may be multiple points in the lift.
            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = End(P)
            sage: f = H([-5*x^2 + 4*y^2,4*x*y])
            sage: f.lift_to_rational_periodic([[P(1,0).change_ring(GF(3)),1]]) # long time
            [[(1 : 0), 1], [(2/3 : 1), 1], [(-2/3 : 1), 1]]

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = End(P)
            sage: f = H([16*x^2 - 29*y^2,16*y^2])
            sage: f.lift_to_rational_periodic([[P(3,1).change_ring(GF(13)), 3]])
            [[(-1/4 : 1), 3]]

        ::

            sage: P.<x,y,z> = ProjectiveSpace(QQ,2)
            sage: H = End(P)
            sage: f = H([76*x^2 - 180*x*y + 45*y^2 + 14*x*z + 45*y*z - 90*z^2, 67*x^2 - 180*x*y - 157*x*z + 90*y*z, -90*z^2])
            sage: f.lift_to_rational_periodic([[P(14,19,1).change_ring(GF(23)), 9]]) # long time
            [[(-9 : -4 : 1), 9]]
        """
        if points_modp==[]:
            return([])
        else:
            if B==None:
                B=e**self.height_difference_bound()

            p=points_modp[0][0].codomain().base_ring().characteristic()
            if p==0:
                raise TypeError("Must be positive characteristic")
            PS=self.domain()
            N=PS.dimension_relative()
            R=RealField()
            #compute the maximum p-adic precision needed to conclusively determine
            #if the rational point exists
            L=R((R(2**(N/2+1)*sqrt(N+1)*B**2).log())/R(p).log()+1).trunc()

            points=[]
            for i in range(len(points_modp)):
                #[point mod p, period, current p-adic precision]
                points.append([points_modp[i][0].change_ring(QQ,False),points_modp[i][1],1])
            good_points=[]
            #shifts is used in non-Hensel lifting
            shifts=None
            #While there are still points to consider try to lift to next precision
            while points!=[]:
                q=points.pop()
                qindex=N
                #Find the last non-zero coordinate to use for normalizations
                while q[0][qindex]%p==0:
                    qindex-=1
                T=q[0]
                n=q[1]
                k=q[2]
                T.scale_by(1/T[qindex]) #normalize
                bad=0
                #stop where we reach the needed precision or the point is bad
                while k< L and bad==0:
                    l=self._multipliermod(T,n,p,2*k)
                    l-=l.parent().one() #f^n(x) - x
                    lp=l.change_ring(Zmod(p**k))
                    ldet=lp.determinant()
                    # if the matrix is invertible then we can Hensel lift
                    if ldet%p!=0:
                        RQ=ZZ.quo(p**(2*k))
                        T.clear_denominators()
                        newT = T.change_ring(RQ,False)
                        fp=self.change_ring(RQ,False)
                        S=newT.nth_iterate(fp,n,False).change_ring(QQ,False)
                        T.scale_by(1/T[qindex])
                        S.scale_by(1/S[qindex])
                        for i in range(N+1):
                            S._coords[i]=S._coords[i]-T._coords[i]
                            if S[i]%(p**k) !=0 and i!=N:
                                bad=1
                                break
                        if bad==1:
                            break
                        S.scale_by(-1/p**k)
                        vecs=[Zmod(p**k)(S._coords[iS]) for iS in range(N+1)]
                        vecs.pop(qindex)
                        newvecs=list((lp.inverse())*vector(vecs)) #l.inverse should be mod p^k!!
                        newS=[]
                        [newS.append(QQ(newvecs[i])) for i in range(qindex)]
                        newS.append(0)
                        [newS.append(QQ(newvecs[i])) for i in range(qindex,N)]
                        S=PS.point(newS, False) #don't check for [0,...,0]
                        for i in range(N+1):
                            S._coords[i]=S._coords[i]%(p**k)
                        for i in range(N+1):
                            T._coords[i]+=S._coords[i]*(p**k)
                        T.normalize_coordinates()
                        #Hensel gives us 2k for the newprecision
                        k=min(2*k,L)
                    else:
                        #we are unable to Hensel Lift so must try all possible lifts
                        #to the next precision (k+1)
                        first=0
                        newq=[]
                        RQ=Zmod(p**(k+1))
                        fp=self.change_ring(RQ,False)
                        if shifts is None:
                            shifts=xmrange([p for i in range(N)])
                        for shift in shifts:
                            newT=T.change_ring(RQ,False)
                            shiftindex=0
                            for i in range(N+1):
                                if i != qindex:
                                    newT._coords[i]=newT[i]+shift[shiftindex]*p**k
                                    shiftindex+=1
                            TT=fp.nth_iterate(newT,n,False)
                            if TT== newT:
                                if first==0:
                                    newq.append(newT.change_ring(QQ,False))
                                    newq.append(n)
                                    newq.append(k+1)
                                    first=1
                                else:
                                    points.append([newT.change_ring(QQ,False),n,k+1])
                        if newq==[]:
                            bad=1
                            break
                        else:
                            T=newq[0]
                            k+=1
                #given a p-adic lift of appropriate precision
                #perform LLL to find the "smallest" rational approximation
                #If this height is small enough, then it is a valid rational point
                if bad==0:
                    M=matrix(N+2,N+1)
                    T.clear_denominators()
                    for i in range(N+1):
                        M[0,i]=T[i]
                        M[i+1,i]=p**L
                    M[N+1,N]=p**L
                    M=M.LLL()
                    Q=[]
                    [Q.append(M[1,i]) for i in range(N+1)]
                    g=gcd(Q)
                    #remove gcds since this is a projective point
                    newB=B*g
                    for i in range(N+1):
                        if abs(Q[i]) >newB:
                            #height too big, so not a valid point
                            bad=1
                            break
                    if bad==0:
                        P=PS.point(Q,False)
                        #check that it is actually periodic
                        newP=copy(P)
                        k=1
                        done=False
                        while done==False and k <= n:
                              newP=self(newP)
                              if newP==P:
                                  if ([P,k] in good_points)==False:
                                      good_points.append([newP,k])
                                  done=True
                              k+=1

            return(good_points)

    def rational_periodic_points(self,**kwds):
        r"""
        Determine the set of rational periodic points for self an endomorphism of projective space.
        Must be defined over `\QQ`.

        The default parameter values are typically good choices for `\mathbb{P}^1`. If you are having
        trouble getting a partiuclar map to finish, try first computing the possible periods, then
        try various different ``lifting_prime``.

        ALGORITHM:
            Modulo each prime of good reduction `p` determine the set of periodic points modulo `p`.
            For each cycle modulo `p` compute the set of possible periods (`mrp^e`). Take the intersection
            of the list of possible periods modulo several primes of good reduction to get a possible list
            of minimal periods of rational periodic points. Take each point modulo `p` associated to each
            of these possible periods and try to lift it to a rational point with a combination of
            `p`-adic approximation and the LLL basis reducion algorithm.

            See B. Hutz, Determination of all rational preperiodic points for morphisms of Pn, submitted, 2012.

        INPUT:

        kwds:

        - ``prime_bound`` - a pair (list or tuple) of positive integers that represent the
            limits of primes to use in the reduction step. Or an integer that represents the upper bound. (optional)
            default: [1,20]

        -  ``lifting_prime`` - a prime integer. (optional) argument that specifies modulo which prime to try and perform the
            lifting. default: 23

        - ``periods`` - a list of positive integers which is the list of possible periods. (optional)

        - ``bad_primes`` - a list or tuple of integer primes, the primes of bad reduction.  (optional)

        OUTPUT:

        - a list of rational points in projective space.

        Examples::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = End(P)
            sage: f = H([x^2-3/4*y^2,y^2])
            sage: f.rational_periodic_points(prime_bound=20,lifting_prime=7) # long time
            [(-1/2 : 1), (1 : 0), (3/2 : 1)]

        ::

            sage: P.<x,y,z> = ProjectiveSpace(QQ,2)
            sage: H = End(P)
            sage: f = H([2*x^3 - 50*x*z^2 + 24*z^3,5*y^3 - 53*y*z^2 + 24*z^3,24*z^3])
            sage: f.rational_periodic_points(prime_bound=[1,20]) # long time
            [(-3 : 1 : 1), (3 : 1 : 1), (5 : 1 : 1), (-1 : 0 : 1), (3 : 3 : 1), (-3
            : 3 : 1), (-1 : 3 : 1), (1 : 3 : 1), (-3 : -1 : 1), (5 : 3 : 1), (-1 :
            -1 : 1), (1 : 1 : 1), (3 : 0 : 1), (-3 : 0 : 1), (5 : 0 : 1), (3 : -1 :
            1), (1 : 0 : 0), (5 : -1 : 1), (-1 : 1 : 1), (1 : -1 : 1), (0 : 1 : 0),
            (1 : 0 : 1)]

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = End(P)
            sage: f = H([-5*x^2 + 4*y^2,4*x*y])
            sage: f.rational_periodic_points() # long time
            [(2/3 : 1), (-2 : 1), (1 : 0), (2 : 1), (-2/3 : 1)]

        .. TODO::

            - move some of this to Cython so that it is faster especially the possible periods mod `p`.

            - have the last prime of good redution used also return the list of points instead of getting the
              information again for all_points.
        """
        if not self.is_endomorphism():
            raise NotImplementedError("Must be an endomorphism of projective space")
        if self.domain().base_ring()!=QQ:
            raise NotImplementedError("Must be QQ") #for p-adic lifting

        primebound = kwds.pop("prime_bound",[1,20])
        p = kwds.pop("lifting_prime",23)
        periods = kwds.pop("periods",None)
        badprimes = kwds.pop("bad_primes",None)

        if (isinstance(primebound,(list,tuple))==False):
            try:
                primebound=[1,ZZ(primebound)]
            except TypeError:
                raise TypeError("Bound on primes must be an integer")
        else:
            try:
                primebound[0]=ZZ(primebound[0])
                primebound[1]=ZZ(primebound[1])
            except TypeError:
                raise TypeError("Prime bounds must be integers")

        if badprimes==None:
            badprimes=self.primes_of_bad_reduction()
        if periods==None:
            periods=self.possible_periods(prime_bound=primebound,bad_primes=badprimes)
        PS=self.domain()
        R=PS.base_ring()
        periodic=set()
        while p in badprimes:
            p = next_prime(p+1)
        B=e**self.height_difference_bound()

        f=self.change_ring(GF(p))
        all_points=f.possible_periods(True) #return the list of points and their periods.
        pos_points=[]
        for i in range(len(all_points)):
            if all_points[i][1] in periods and  (all_points[i] in pos_points)==False:  #check period, remove duplicates
                pos_points.append(all_points[i])
        periodic_points=self.lift_to_rational_periodic(pos_points,B)
        for P,n in periodic_points:
            for k in range(n):
                  P.normalize_coordinates()
                  periodic.add(P)
                  P=self(P)
        return(list(periodic))

    def rational_preimages(self,Q):
        r"""
        Given a rational point `Q` in the domain of ``self``, return all the rational points `P`
        in the domain of ``self`` with `self(P)==Q`. In other words, the set of first pre-images of `Q`.
        ``self`` must be defined over `\QQ` and be an endomorphism of projective space.

        ALGORITHM:
            Use elimination via groebner bases to find the rational pre-images

        INPUT:

        - ``Q`` - a rational point in the domain of ``self``.

        OUTPUT:

        - a list of rational points in the domain of ``self``.

        Examples::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = End(P)
            sage: f = H([16*x^2 - 29*y^2,16*y^2])
            sage: f.rational_preimages(P(-1,4))
            [(5/4 : 1), (-5/4 : 1)]

        ::

            sage: P.<x,y,z> = ProjectiveSpace(QQ,2)
            sage: H = End(P)
            sage: f = H([76*x^2 - 180*x*y + 45*y^2 + 14*x*z + 45*y*z - 90*z^2, 67*x^2 - 180*x*y - 157*x*z + 90*y*z, -90*z^2])
            sage: f.rational_preimages(P(-9,-4,1))
            [(0 : 4 : 1)]

        ::

            A non-periodic example.
            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = End(P)
            sage: f = H([x^2 + y^2,2*x*y])
            sage: f.rational_preimages(P(17,15))
            [(5/3 : 1), (3/5 : 1)]
        """
        if not self.is_endomorphism():
            raise NotImplementedError("Must be an endomorphism of projective space")
        if self.domain().base_ring()!=QQ:
            raise NotImplementedError("Must be QQ")
        PS=self.domain()
        R=PS.coordinate_ring()
        N=PS.dimension_relative()
        #need a lexicographic ordering for elimination
        R=PolynomialRing(R.base_ring(),N+1,R.gens(),order='lex')
        I=[]
        preimages=set()
        for i in range(N+1):
            for j in range(i+1,N+1):
                I.append(Q[i]*self[j]-Q[j]*self[i])
        I = I*R
        #Determine the points through elimination
        #This is much faster than using the I.variety() function on each affine chart.
        for k in range(N+1):
            #create the elimination ideal for the kth affine patch
            G=I.substitute({R.gen(k):1}).groebner_basis()
            if G!=[1]:
                P={}
                #keep track that we know the kth coordinate is 1
                P.update({R.gen(k):1})
                points=[P]
                #work backwards from solving each equation for the possible
                #values of the next coordiante
                for i in range(len(G)-1,-1,-1):
                    new_points=[]
                    good=0
                    for P in points:
                        #subsitute in our dictionary entry that has the values
                        #of coordinates known so far. This results in a single
                        #variable polynomial (by elimination)
                        L=G[i].substitute(P)
                        if L!=0:
                            L=L.factor()
                            #the linear factors give the possible rational values of
                            #this coordinate
                            for pol,pow in L:
                                if pol.degree()==1 and len(pol.variables()) ==1:
                                    good=1
                                    r=pol.variable()
                                    varindex=R.gens().index(r)
                                    #add this coordinates information to th
                                    #each dictionary entry
                                    P.update({R.gen(varindex):-pol.coefficient({r:0})/pol.coefficient({r:1})})
                                    new_points.append(copy(P))
                    if good==1:
                        points=new_points
                #the dictionary entries now have values for all coordiantes
                #they are the rational solutions to the equations
                #make them into projective points
                for i in range(len(points)):
                    if len(points[i])==N+1:
                        S=PS([points[i][R.gen(j)] for j in range(N+1)])
                        S.normalize_coordinates()
                        preimages.add(S)
        return(list(preimages))


    def all_rational_preimages(self,points):
        r"""
        Given a set of rational points in the domain of ``self``, return all the rational
        pre-images of those points. In others words, all the rational points which have some
        iterate in the set points. This function repeatedly calls ``rational_preimages``.
        If the degree is at least two, by Northocott, this is always a finite set.
        ``self`` must be defined over `\QQ` and be an endomorphism of projective space.

        INPUT:

        - ``points`` - a list of rational points in the domain of ``self``

        OUTPUT:

        - a list of rational points in the domain of ``self``.

        Examples::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = End(P)
            sage: f = H([16*x^2 - 29*y^2,16*y^2])
            sage: f.all_rational_preimages([P(-1,4)])
            [(-1/4 : 1), (5/4 : 1), (-3/4 : 1), (3/4 : 1), (-5/4 : 1), (1/4 : 1),
            (7/4 : 1), (-7/4 : 1)]

        ::

            sage: P.<x,y,z> = ProjectiveSpace(QQ,2)
            sage: H = End(P)
            sage: f = H([76*x^2 - 180*x*y + 45*y^2 + 14*x*z + 45*y*z - 90*z^2, 67*x^2 - 180*x*y - 157*x*z + 90*y*z, -90*z^2])
            sage: f.all_rational_preimages([P(-9,-4,1)])
            [(-9 : -4 : 1), (1 : 3 : 1), (0 : 4 : 1), (0 : -1 : 1), (0 : 0 : 1), (1
            : 1 : 1), (0 : 1 : 1), (1 : 2 : 1), (1 : 0 : 1)]

        ::

            A non-periodic example.
            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = End(P)
            sage: f = H([x^2 + y^2,2*x*y])
            sage: f.all_rational_preimages([P(17,15)])
            [(5/3 : 1), (1/3 : 1), (3/5 : 1), (3 : 1)]
        """
        if not self.is_endomorphism():
            raise NotImplementedError("Must be an endomorphism of projective space")
        if self.domain().base_ring()!=QQ:
            raise NotImplementedError("Must be QQ")

        PS=self.domain()
        RPS=PS.base_ring()
        all_preimages=set()
        while points!=[]:
            P=points.pop()
            preimages=self.rational_preimages(P)
            for i in range(len(preimages)):
                if not preimages[i] in all_preimages:
                    points.append(preimages[i])
                    all_preimages.add(preimages[i])
        return(list(all_preimages))

    def rational_preperiodic_points(self,**kwds):
        r"""
        Determined the set of rational preperiodic points for ``self``.
        ``self`` must be defined over `\QQ` and be an endomorphism of projective space.

        The default parameter values are typically good choices for `\mathbb{P}^1`. If you are having
        trouble getting a partiuclar map to finish, try first computing the possible periods, then
        try various different ``lifting_prime``.

        ALGORITHM:

        - Determines the list of possible periods.

        - Determines the rational periodic points from the possible periods.

        - Determines the rational preperiodic points from the rational periodic points
          by determining rational preimages.

        INPUT:

        kwds:

        - ``prime_bound`` - a pair (list or tuple) of positive integers that represent the
          limits of primes to use in the reduction step. Or an integer that represents the upper bound. (optional)
          default: [1,20]

        - ``lifting_prime`` - a prime integer. (optional) argument that specifies modulo which prime to try and perform the
          lifting. default: 23

        - ``periods`` - a list of positive integers which is the list of possible periods. (optional)

        - ``bad_primes`` - a list or tuple of integer primes, the primes of bad reduction.  (optional)

        OUTPUT:

        - a list of rational points in projective space.

        Examples::

            sage: PS.<x,y> = ProjectiveSpace(1,QQ)
            sage: H = End(PS)
            sage: f = H([x^2 -y^2,3*x*y])
            sage: f.rational_preperiodic_points()
            [(-2 : 1), (0 : 1), (-1/2 : 1), (1 : 0), (1 : 1), (2 : 1), (-1 : 1),
            (1/2 : 1)]

        ::

            sage: PS.<x,y> = ProjectiveSpace(1,QQ)
            sage: H = End(PS)
            sage: f = H([5*x^3 - 53*x*y^2 + 24*y^3, 24*y^3])
            sage: f.rational_preperiodic_points(prime_bound=10)
            [(1 : 1), (1 : 0), (-1 : 1), (3 : 1), (0 : 1)]

        ::

            sage: PS.<x,y,z> = ProjectiveSpace(2,QQ)
            sage: H = End(PS)
            sage: f = H([x^2 - 21/16*z^2,y^2-2*z^2,z^2])
            sage: f.rational_preperiodic_points(prime_bound=[1,8],lifting_prime=7,periods=[2]) # long time
            [(5/4 : 1 : 1), (1/4 : 1 : 1), (1/4 : -2 : 1), (1/4 : 2 : 1), (5/4 : -1 : 1), (5/4 : 2 : 1),
            (-5/4 : 2 : 1), (1/4 : -1 : 1), (-1/4 : 2 : 1), (-1/4 : -2 : 1), (-5/4 : 1 : 1), (5/4 : 0 : 1),
            (-5/4 : 0 : 1), (-1/4 : 1 : 1), (-1/4 : 0 : 1), (-5/4 : -1 : 1), (5/4 : -2 : 1), (-1/4 : -1 : 1),
            (-5/4 : -2 : 1), (1/4 : 0 : 1)]
        """
        #input error checking done in possible_periods and rational_peridic_points
        badprimes = kwds.pop("bad_primes",None)
        periods = kwds.pop("periods",None)
        primebound = kwds.pop("prime_bound",[1,20])
        if badprimes==None:
            badprimes=self.primes_of_bad_reduction()
        if periods==None:
            periods=self.possible_periods(prime_bound=primebound,bad_primes=badprimes) #determine the set of possible periods
        if periods==[]:
            return([]) #no rational preperiodic points
        else:
            p = kwds.pop("lifting_prime",23)
            T=self.rational_periodic_points(prime_bound=primebound,lifting_prime=p,periods=periods,bad_primes=badprimes) #find the rationla preperiodic points
            preper=self.all_rational_preimages(T) #find the preperiodic points
            preper=list(preper)
            return(preper)

    def rational_preperiodic_graph(self,**kwds):
        r"""
        Determines the set of rational preperiodic points for ``self``.
        self must be defined over `\QQ` and be an endomorphism of projective space.

        ALGORITHM:
        - Determines the list of possible periods.

        - Determines the rational periodic points from the possible periods.

        - Determines the rational preperiodic points from the rational periodic points
          by determining rational preimages.


        INPUT:

        kwds:

        - ``prime_bound`` - a pair (list or tuple) of positive integers that represent the
            limits of primes to use in the reduction step. Or an integer that represents the upper bound. (optional)
            default: [1,20]

        -  ``lifting_prime`` - a prime integer. (optional) argument that specifies modulo which prime to try and perform the
            lifting. default: 23

        - ``periods`` - a list of positive integers which is the list of possible periods. (optional)

        - ``bad_primes`` - a list or tuple of integer primes, the primes of bad reduction.  (optional)

        OUTPUT:

        - a digraph representing the orbits of the rational preperiodic points in projective space.

        Examples::

            sage: PS.<x,y> = ProjectiveSpace(1,QQ)
            sage: H = End(PS)
            sage: f = H([7*x^2 - 28*y^2,24*x*y])
            sage: f.rational_preperiodic_graph()
            Looped digraph on 12 vertices

        ::

            sage: PS.<x,y> = ProjectiveSpace(1,QQ)
            sage: H = End(PS)
            sage: f = H([-3/2*x^3 +19/6*x*y^2,y^3])
            sage: f.rational_preperiodic_graph(prime_bound=[1,8])
            Looped digraph on 12 vertices

        ::

            sage: PS.<x,y,z> = ProjectiveSpace(2,QQ)
            sage: H = End(PS)
            sage: f = H([2*x^3 - 50*x*z^2 + 24*z^3,5*y^3 - 53*y*z^2 + 24*z^3,24*z^3])
            sage: f.rational_preperiodic_graph(prime_bound=[1,11],lifting_prime=13) # long time
            Looped digraph on 30 vertices
        """
        #input checking done in .rational_preperiodic_points()
        preper=self.rational_preperiodic_points(**kwds)
        g=self._preperiodic_points_to_cyclegraph(preper)
        return(g)


class SchemeMorphism_polynomial_projective_space_finite_field(SchemeMorphism_polynomial_projective_space_field):

    def orbit_structure(self, P):
        r"""
        Every point is preperiodic over a finite field. This funtion returns the pair `[m,n]` where `m` is the
        preperiod and `n` the period of the point ``P`` by ``self``.

        INPUT:

        - ``P`` -- a point in ``self.domain()``

        OUTPUT:

        - a list `[m,n]` of integers

        EXAMPLES::

            sage: P.<x,y,z> = ProjectiveSpace(GF(5),2)
            sage: H = Hom(P,P)
            sage: f = H([x^2+y^2,y^2,z^2 + y*z])
            sage: f.orbit_structure(P(2,1,2))
            [0, 6]

        ::

            sage: P.<x,y,z> = ProjectiveSpace(GF(7),2)
            sage: X = P.subscheme(x^2-y^2)
            sage: H = Hom(X,X)
            sage: f = H([x^2,y^2,z^2])
            sage: f.orbit_structure(X(1,1,2))
            [0, 2]

        ::

            sage: P.<x,y> = ProjectiveSpace(GF(13),1)
            sage: H = Hom(P,P)
            sage: f = H([x^2-y^2,y^2])
            sage: f.orbit_structure(P(3,4))
            [2, 3]
        """
        return(P.orbit_structure(self))

    def cyclegraph(self):
        r"""
        returns Digraph of all orbits of ``self`` mod `p`.

        For subschemes, only points on the subscheme whose image are
        also on the subscheme are in the digraph.

        OUTPUT:

        - a digraph

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(GF(13),1)
            sage: H = Hom(P,P)
            sage: f = H([x^2-y^2,y^2])
            sage: f.cyclegraph()
            Looped digraph on 14 vertices

        ::

            sage: P.<x,y,z> = ProjectiveSpace(GF(5^2,'t'),2)
            sage: H = Hom(P,P)
            sage: f = H([x^2+y^2,y^2,z^2+y*z])
            sage: f.cyclegraph()
            Looped digraph on 651 vertices

        ::

            sage: P.<x,y,z> = ProjectiveSpace(GF(7),2)
            sage: X = P.subscheme(x^2-y^2)
            sage: H = Hom(X,X)
            sage: f = H([x^2,y^2,z^2])
            sage: f.cyclegraph()
            Looped digraph on 15 vertices
        """
        if self.domain() != self.codomain():
            raise NotImplementedError("Domain and Codomain must be equal")
        V=[]
        E=[]
        from sage.schemes.projective.projective_space import is_ProjectiveSpace
        if is_ProjectiveSpace(self.domain()) is True:
            for P in self.domain():
                V.append(str(P))
                Q=self(P)
                Q.normalize_coordinates()
                E.append([str(Q)])
        else:
            X=self.domain()
            for P in X.ambient_space():
                try:
                    XP=X.point(P)
                    V.append(str(XP))
                    Q=self(XP)
                    Q.normalize_coordinates()
                    E.append([str(Q)])
                except TypeError:  # not a point on the scheme
                    pass
        from sage.graphs.digraph import DiGraph
        g=DiGraph(dict(zip(V,E)), loops=True)
        return g

    def possible_periods(self,return_points=False):
        r"""
        Returns the list of possible minimal periods of a periodic point
        over `\QQ` and (optionally) a point in each cycle.

        ALGORITHM:

        The list comes from: Hutz, Good reduction of periodic points, Illinois Journal of
        Mathematics 53 (Winter 2009), no. 4, 1109-1126.

        INPUT:

        - ``return_points`` - Boolean (optional) - a value of True returns the points as well as the possible periods.

        OUTPUT:

        - a list of positive integers, or a list of pairs of projective points and periods if ``flag`` is 1.

        Examples::

            sage: P.<x,y> = ProjectiveSpace(GF(23),1)
            sage: H = End(P)
            sage: f = H([x^2-2*y^2,y^2])
            sage: f.possible_periods()
            [1, 5, 11, 22, 110]

        ::

            sage: P.<x,y> = ProjectiveSpace(GF(13),1)
            sage: H = End(P)
            sage: f = H([x^2-y^2,y^2])
            sage: f.possible_periods(True)
            [[(1 : 0), 1], [(0 : 1), 2], [(3 : 1), 3], [(3 : 1), 36]]

        ::

            sage: PS.<x,y,z> = ProjectiveSpace(2,GF(7))
            sage: H = End(PS)
            sage: f = H([-360*x^3 + 760*x*z^2, y^3 - 604*y*z^2 + 240*z^3, 240*z^3])
            sage: f.possible_periods()
            [1, 2, 4, 6, 12, 14, 28, 42, 84]

        .. TODO::

            - do not reutrn duplicate points

            - check == False to speed up?

            - move to Cython

        """
        if not is_PrimeFiniteField(self.domain().base_ring()):
            raise TypeError("Must be prime field")
        if not self.is_endomorphism():
            raise NotImplementedError("Must be an endomorphism of projective space")

        PS=self.domain()
        p=PS.base_ring().order()
        N=PS.dimension_relative()
        pointsdict=PS.rational_points_dictionary() #assume p is prime
        pointslist=list(pointsdict)
        hashlist=pointsdict.values()
        pointtable=[[0,0] for i in range(len(pointsdict))]
        index=1
        periods=set()
        points_periods=[]
        for j in range(len(pointsdict)):
            hashP=hashlist[j]
            if pointtable[hashP][1]==0:
                startindex=index
                P=pointslist[j]
                while pointtable[hashP][1]==0:
                    pointtable[hashP][1]=index
                    Q=self(P)
                    Q.normalize_coordinates()
                    hashQ=pointsdict[Q]
                    pointtable[hashP][0]=hashQ
                    P=Q
                    hashP=hashQ
                    index+=1
                if pointtable[hashP][1]>= startindex:
                    period=index-pointtable[hashP][1]
                    periods.add(period)
                    points_periods.append([P,period])
                    l=P.multiplier(self,period,False)
                    lorders=set()
                    for poly,_ in l.charpoly().factor():
                        if poly.degree() == 1:
                            eig = -poly.constant_coefficient()
                            if not eig:
                                continue # exclude 0
                        else:
                            eig = GF(p ** poly.degree(), 't', modulus=poly).gen()
                        if eig:
                            lorders.add(eig.multiplicative_order())
                    S = subsets(lorders)
                    S.next()   # get rid of the empty set
                    rvalues=set()
                    for s in S:
                        rvalues.add(lcm(s))
                    rvalues=list(rvalues)
                    if N==1:
                        for k in range(len(rvalues)):
                            r=rvalues[k]
                            periods.add(period*r)
                            points_periods.append([P,period*r])
                            if p==2 or p==3: #need e=1 for N=1, QQ
                                periods.add(period*r*p)
                                points_periods.append([P,period*r*p])
                    else:
                        for k in range(len(rvalues)):
                            r=rvalues[k]
                            periods.add(period*r)
                            periods.add(period*r*p)
                            points_periods.append([P,period*r])
                            points_periods.append([P,period*r*p])
                            if p==2:  #need e=3 for N>1, QQ
                                periods.add(period*r*4)
                                points_periods.append([P,period*r*4])
                                periods.add(period*r*8)
                                points_periods.append([P,period*r*8])

        if return_points==False:
            return(sorted(periods))
        else:
            return(points_periods)

