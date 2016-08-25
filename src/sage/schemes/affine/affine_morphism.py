r"""
Morphisms on affine varieties

A morphism of schemes determined by rational functions that define \
what the morphism does on points in the ambient affine space.


AUTHORS:

- David Kohel, William Stein

- Volker Braun (2011-08-08): Renamed classes, more documentation, misc
  cleanups.

- Ben Hutz (2013-03) iteration functionality and new directory structure
  for affine/projective
"""

#*****************************************************************************
#       Copyright (C) 2011 Volker Braun <vbraun.name@gmail.com>
#       Copyright (C) 2006 David Kohel <kohel@maths.usyd.edu.au>
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


from sage.calculus.functions import jacobian
from sage.categories.fields import Fields
from sage.categories.homset import Hom
from sage.matrix.constructor import matrix, identity_matrix
from sage.misc.cachefunc import cached_method
from sage.misc.all import prod
from sage.rings.all import Integer
from sage.arith.all import lcm, gcd
from sage.rings.complex_field import ComplexField
from sage.rings.finite_rings.finite_field_constructor import GF, is_PrimeFiniteField
from sage.rings.fraction_field import FractionField
from sage.rings.fraction_field_element import FractionFieldElement
from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.real_mpfr import RealField
from sage.schemes.generic.morphism import SchemeMorphism_polynomial
from sage.misc.lazy_attribute import lazy_attribute
from sage.ext.fast_callable import fast_callable
import sys

class SchemeMorphism_polynomial_affine_space(SchemeMorphism_polynomial):
    """
    A morphism of schemes determined by rational functions.

    EXAMPLES::

        sage: RA.<x,y> = QQ[]
        sage: A2 = AffineSpace(RA)
        sage: RP.<u,v,w> = QQ[]
        sage: P2 = ProjectiveSpace(RP)
        sage: H = A2.Hom(P2)
        sage: f = H([x, y, 1])
        sage: f
        Scheme morphism:
          From: Affine Space of dimension 2 over Rational Field
          To:   Projective Space of dimension 2 over Rational Field
          Defn: Defined on coordinates by sending (x, y) to
                (x : y : 1)
    """
    def __init__(self, parent, polys, check=True):
        r"""
        The Python constructor.

        See :class:`SchemeMorphism_polynomial` for details.

        INPUT:

        - ``parent`` -- Hom.

        - ``polys`` -- list or tuple of polynomial or rational functions.

        - ``check`` -- Boolean.

        OUTPUT:

        - :class:`SchemeMorphism_polynomial_affine_space`.

        EXAMPLES::

            sage: A.<x,y> = AffineSpace(ZZ, 2)
            sage: H = Hom(A, A)
            sage: H([3/5*x^2, y^2/(2*x^2)])
            Traceback (most recent call last):
            ...
            TypeError: polys (=[3/5*x^2, y^2/(2*x^2)]) must be rational functions in
            Multivariate Polynomial Ring in x, y over Integer Ring

        ::

            sage: A.<x,y> = AffineSpace(ZZ, 2)
            sage: H = Hom(A, A)
            sage: H([3*x^2/(5*y), y^2/(2*x^2)])
            Scheme endomorphism of Affine Space of dimension 2 over Integer Ring
              Defn: Defined on coordinates by sending (x, y) to
                    (3*x^2/(5*y), y^2/(2*x^2))

        ::

            sage: A.<x,y> = AffineSpace(QQ, 2)
            sage: H = Hom(A, A)
            sage: H([3/2*x^2, y^2])
            Scheme endomorphism of Affine Space of dimension 2 over Rational Field
              Defn: Defined on coordinates by sending (x, y) to
                    (3/2*x^2, y^2)

        ::

            sage: A.<x,y> = AffineSpace(QQ, 2)
            sage: X = A.subscheme([x-y^2])
            sage: H = Hom(X, X)
            sage: H([9/4*x^2, 3/2*y])
            Scheme endomorphism of Closed subscheme of Affine Space of dimension 2
            over Rational Field defined by:
              -y^2 + x
              Defn: Defined on coordinates by sending (x, y) to
                    (9/4*x^2, 3/2*y)

            sage: P.<x,y,z> = ProjectiveSpace(ZZ, 2)
            sage: H = Hom(P, P)
            sage: f = H([5*x^3 + 3*x*y^2-y^3, 3*z^3 + y*x^2, x^3-z^3])
            sage: f.dehomogenize(2)
            Scheme endomorphism of Affine Space of dimension 2 over Integer Ring
              Defn: Defined on coordinates by sending (x0, x1) to
                    ((5*x0^3 + 3*x0*x1^2 - x1^3)/(x0^3 - 1), (x0^2*x1 + 3)/(x0^3 - 1))

            If you pass in quotient ring elements, they are reduced::

            sage: A.<x,y,z> = AffineSpace(QQ, 3)
            sage: X = A.subscheme([x-y])
            sage: H = Hom(X,X)
            sage: u,v,w = X.coordinate_ring().gens()
            sage: H([u, v, u+v])
            Scheme endomorphism of Closed subscheme of Affine Space of dimension 3
            over Rational Field defined by:
              x - y
              Defn: Defined on coordinates by sending (x, y, z) to
                    (y, y, 2*y)

            You must use the ambient space variables to create rational functions::

            sage: A.<x,y,z> = AffineSpace(QQ, 3)
            sage: X = A.subscheme([x^2-y^2])
            sage: H = Hom(X,X)
            sage: u,v,w = X.coordinate_ring().gens()
            sage: H([u, v, (u+1)/v])
            Traceback (most recent call last):
            ...
            ArithmeticError: Division failed. The numerator is not a multiple of the denominator.
            sage: H([x, y, (x+1)/y])
            Scheme endomorphism of Closed subscheme of Affine Space of dimension 3
            over Rational Field defined by:
              x^2 - y^2
              Defn: Defined on coordinates by sending (x, y, z) to
                    (x, y, (x + 1)/y)

            ::

            sage: R.<t> = PolynomialRing(QQ)
            sage: A.<x,y,z> = AffineSpace(R, 3)
            sage: X = A.subscheme(x^2-y^2)
            sage: H = End(X)
            sage: H([x^2/(t*y), t*y^2, x*z])
            Scheme endomorphism of Closed subscheme of Affine Space of dimension 3
            over Univariate Polynomial Ring in t over Rational Field defined by:
              x^2 - y^2
              Defn: Defined on coordinates by sending (x, y, z) to
                    (x^2/(t*y), t*y^2, x*z)
        """
        if check:
            if not isinstance(polys, (list, tuple)):
                raise TypeError("polys (=%s) must be a list or tuple"%polys)
            source_ring = parent.domain().ambient_space().coordinate_ring()
            target = parent.codomain().ambient_space()
            if len(polys) != target.ngens():
                raise ValueError("there must be %s polynomials"%target.ngens())
            try:
                polys = [source_ring(poly) for poly in polys]
            except TypeError: #maybe given quotient ring elements
                try:
                   polys = [source_ring(poly.lift()) for poly in polys]
                except (TypeError, AttributeError):
                    #must be a rational function since we cannot have
                    #rational functions for quotient rings
                    try:
                        if not all(p.base_ring()==source_ring.base_ring() for p in polys):
                            raise TypeError("polys (=%s) must be rational functions in %s"%(polys, source_ring))
                        polys = [source_ring(poly.numerator())/source_ring(poly.denominator()) for poly in polys]
                    except TypeError: #can't seem to coerce
                        raise TypeError("polys (=%s) must be rational functions in %s"%(polys, source_ring))
        self._is_prime_finite_field = is_PrimeFiniteField(polys[0].base_ring()) # Needed for _fast_eval and _fastpolys
        SchemeMorphism_polynomial.__init__(self, parent, polys, False)

    def __call__(self, x, check=True):
        """
        Evaluate affine morphism at point described by ``x``.

        EXAMPLES::

            sage: P.<x,y,z> = AffineSpace(QQ, 3)
            sage: H = Hom(P, P)
            sage: f = H([x^2+y^2, y^2, z^2 + y*z])
            sage: f(P([1, 1, 1]))
            (2, 1, 2)
        """
        from sage.schemes.affine.affine_point import SchemeMorphism_point_affine
        if check:
            if not isinstance(x, SchemeMorphism_point_affine):
                try:
                    x = self.domain()(x)
                except (TypeError, NotImplementedError):
                    raise TypeError("%s fails to convert into the map's domain %s, but a `pushforward` method is not properly implemented"%(x, self.domain()))
            elif self.domain() != x.codomain():
                raise TypeError("%s fails to convert into the map's domain %s,but a `pushforward` method is not properly implemented"%(x, self.domain()))

        # Passes the array of args to _fast_eval
        P = self._fast_eval(x._coords)
        return self.codomain().point(P, check)

    def __eq__(self, right):
        """
        Tests the equality of two affine maps.

        INPUT:

        - ``right`` - a map on affine space.

        OUTPUT:

        - Boolean - True if the two affine maps define the same map.

        EXAMPLES::

            sage: A.<x,y> = AffineSpace(QQ, 2)
            sage: A2.<u,v> = AffineSpace(QQ, 2)
            sage: H = End(A)
            sage: H2 = End(A2)
            sage: f = H([x^2 - 2*x*y, y/(x+1)])
            sage: g = H2([u^3 - v, v^2])
            sage: f == g
            False

            ::

            sage: A.<x,y,z> = AffineSpace(CC, 3)
            sage: H = End(A)
            sage: f = H([x^2 - CC.0*x*y + z*x, 1/z^2 - y^2, 5*x])
            sage: f == f
            True
        """
        if not isinstance(right, SchemeMorphism_polynomial):
            return False
        if self.parent() != right.parent():
            return False
        return all([self[i] == right[i] for i in range(len(self._polys))])

    def __ne__(self, right):
        """
        Tests the inequality of two affine maps.

        INPUT:

        - ``right`` -  a map on affine space.

        OUTPUT:

        - Boolean - True if the two affine maps define the same map.

        EXAMPLES::

            sage: A.<x,y> = AffineSpace(RR, 2)
            sage: H = End(A)
            sage: f = H([x^2 - y, y^2])
            sage: g = H([x^3-x*y, x*y^2])
            sage: f != g
            True
            sage: f != f
            False
        """
        if not isinstance(right, SchemeMorphism_polynomial):
            return True
        if self.parent() != right.parent():
            return True
        if all([self[i] == right[i] for i in range(len(self._polys))]):
            return False
        return True

    @lazy_attribute
    def _fastpolys(self):
        """
        Lazy attribute for fast_callable polynomials for affine morphsims.

        EXAMPLES::

            sage: P.<x,y> = AffineSpace(QQ, 2)
            sage: H = Hom(P, P)
            sage: f = H([x^2+y^2, y^2/(1+x)])
            sage: [t.op_list() for g in f._fastpolys for t in g]
            [[('load_const', 0), ('load_const', 1), ('load_arg', 1), ('ipow', 2),
            'mul', 'add', ('load_const', 1), ('load_arg', 0), ('ipow', 2), 'mul',
            'add', 'return'], [('load_const', 0), ('load_const', 1), ('load_arg',
            1), ('ipow', 2), 'mul', 'add', 'return'], [('load_const', 0),
            ('load_const', 1), 'add', 'return'], [('load_const', 0), ('load_const',
            1), ('load_arg', 0), ('ipow', 1), 'mul', 'add', ('load_const', 1),
            'add', 'return']]
        """
        polys = self._polys

        R = self.domain().ambient_space().coordinate_ring()
        # fastpolys[0] corresponds to the numerator polys, fastpolys[1] corresponds to denominator polys
        fastpolys = [[], []]
        for poly in polys:
            # Determine if truly polynomials. Store the numerator and denominator as separate polynomials
            # and repeat the normal process for both.
            try:
                poly_numerator = R(poly)
                poly_denominator = R.one()
            except TypeError:
                poly_numerator = R(poly.numerator())
                poly_denominator = R(poly.denominator())

            # These tests are in place because the float and integer domain evaluate
            # faster than using the base_ring
            if self._is_prime_finite_field:
                prime = polys[0].base_ring().characteristic()
                degree = max(poly_numerator.degree(), poly_denominator.degree())
                height = max([abs(c.lift()) for c in poly_numerator.coefficients()]\
                              + [abs(c.lift()) for c in poly_denominator.coefficients()])
                num_terms = max(len(poly_numerator.coefficients()), len(poly_denominator.coefficients()))
                largest_value = num_terms * height * (prime - 1) ** degree
                # If the calculations will not overflow the float data type use domain float
                # Else use domain integer
                if largest_value < (2 ** sys.float_info.mant_dig):
                    fastpolys[0].append(fast_callable(poly_numerator, domain=float))
                    fastpolys[1].append(fast_callable(poly_denominator, domain=float))
                else:
                    fastpolys[0].append(fast_callable(poly_numerator, domain=ZZ))
                    fastpolys[1].append(fast_callable(poly_denominator, domain=ZZ))
            else:
                fastpolys[0].append(fast_callable(poly_numerator, domain=poly.base_ring()))
                fastpolys[1].append(fast_callable(poly_denominator, domain=poly.base_ring()))
        return fastpolys

    def _fast_eval(self, x):
        """
        Evaluate affine morphism at point described by ``x``.

        EXAMPLES::

            sage: P.<x,y,z> = AffineSpace(QQ, 3)
            sage: H = Hom(P, P)
            sage: f = H([x^2+y^2, y^2, z^2 + y*z])
            sage: f._fast_eval([1, 1, 1])
            [2, 1, 2]

        ::

            sage: P.<x,y,z> = AffineSpace(QQ, 3)
            sage: H = Hom(P, P)
            sage: f = H([x^2/y, y/x, (y^2+z)/(x*y)])
            sage: f._fast_eval([2, 3, 1])
            [4/3, 3/2, 5/3]
        """
        R = self.domain().ambient_space().coordinate_ring()

        P = []
        for i in range(len(self._fastpolys[0])):
            # Check if denominator is the identity;
            #if not, then must append the fraction evaluated at the point
            if self._fastpolys[1][i] is R.one():
                P.append(self._fastpolys[0][i](*x))
            else:
                P.append(self._fastpolys[0][i](*x)/self._fastpolys[1][i](*x))
        return P

    def homogenize(self, n):
        r"""
        Return the homogenization of this map.

        If it's domain is a subscheme, the domain of
        the homogenized map is the projective embedding of the domain. The domain and codomain
        can be homogenized at different coordinates: ``n[0]`` for the domain and ``n[1]`` for the codomain.

        INPUT:

        - ``n`` -- a tuple of nonnegative integers. If ``n`` is an integer,
          then the two values of the tuple are assumed to be the same.

        OUTPUT:

        - :class:`SchemMorphism_polynomial_projective_space`.

        EXAMPLES::

            sage: A.<x,y> = AffineSpace(ZZ, 2)
            sage: H = Hom(A, A)
            sage: f = H([(x^2-2)/x^5, y^2])
            sage: f.homogenize(2)
            Scheme endomorphism of Projective Space of dimension 2 over Integer Ring
              Defn: Defined on coordinates by sending (x0 : x1 : x2) to
                    (x0^2*x2^5 - 2*x2^7 : x0^5*x1^2 : x0^5*x2^2)

        ::

            sage: A.<x,y> = AffineSpace(CC, 2)
            sage: H = Hom(A, A)
            sage: f = H([(x^2-2)/(x*y), y^2-x])
            sage: f.homogenize((2, 0))
            Scheme morphism:
              From: Projective Space of dimension 2 over Complex Field with 53 bits of precision
              To:   Projective Space of dimension 2 over Complex Field with 53 bits of precision
              Defn: Defined on coordinates by sending (x0 : x1 : x2) to
                    (x0*x1*x2^2 : x0^2*x2^2 + (-2.00000000000000)*x2^4 : x0*x1^3 - x0^2*x1*x2)

        ::

            sage: A.<x,y> = AffineSpace(ZZ, 2)
            sage: X = A.subscheme([x-y^2])
            sage: H = Hom(X, X)
            sage: f = H([9*y^2, 3*y])
            sage: f.homogenize(2)
            Scheme endomorphism of Closed subscheme of Projective Space
            of dimension 2 over Integer Ring defined by:
                x1^2 - x0*x2
                Defn: Defined on coordinates by sending (x0 : x1 : x2) to
                      (9*x1^2 : 3*x1*x2 : x2^2)

        ::

            sage: R.<t> = PolynomialRing(ZZ)
            sage: A.<x,y> = AffineSpace(R, 2)
            sage: H = Hom(A, A)
            sage: f = H([(x^2-2)/y, y^2-x])
            sage: f.homogenize((2, 0))
            Scheme morphism:
              From: Projective Space of dimension 2 over Univariate Polynomial Ring in t over Integer Ring
              To:   Projective Space of dimension 2 over Univariate Polynomial Ring in t over Integer Ring
              Defn: Defined on coordinates by sending (x0 : x1 : x2) to
                    (x1*x2^2 : x0^2*x2 + (-2)*x2^3 : x1^3 - x0*x1*x2)

        ::

            sage: A.<x> = AffineSpace(QQ, 1)
            sage: H = End(A)
            sage: f = H([x^2-1])
            sage: f.homogenize((1, 0))
            Scheme morphism:
              From: Projective Space of dimension 1 over Rational Field
              To:   Projective Space of dimension 1 over Rational Field
              Defn: Defined on coordinates by sending (x0 : x1) to
                    (x1^2 : x0^2 - x1^2)

        ::

            R.<a> = PolynomialRing(QQbar)
            A.<x,y> = AffineSpace(R, 2)
            H = End(A)
            f = H([QQbar(sqrt(2))*x*y, a*x^2])
            f.homogenize(2)
            Scheme endomorphism of Projective Space of dimension 2 over Univariate
            Polynomial Ring in a over Algebraic Field
              Defn: Defined on coordinates by sending (x0 : x1 : x2) to
                    (1.414213562373095?*x0*x1 : a*x0^2 : x2^2)

        ::

            sage: P.<x,y,z> = AffineSpace(QQ, 3)
            sage: H = End(P)
            sage: f = H([x^2 - 2*x*y + z*x, z^2 -y^2 , 5*z*y])
            sage: f.homogenize(2).dehomogenize(2) == f
            True

        ::

            sage: K.<c> = FunctionField(QQ)
            sage: A.<x> = AffineSpace(K, 1)
            sage: f = Hom(A, A)([x^2 + c])
            sage: f.homogenize(1)
            Scheme endomorphism of Projective Space of
            dimension 1 over Rational function field in c over Rational Field
              Defn: Defined on coordinates by sending (x0 : x1) to
                    (x0^2 + c*x1^2 : x1^2)

        ::

            sage: A.<z> = AffineSpace(QQbar, 1)
            sage: H = End(A)
            sage: f = H([2*z / (z^2+2*z+3)])
            sage: f.homogenize(1)
            Scheme endomorphism of Projective Space of dimension 1 over Algebraic
            Field
              Defn: Defined on coordinates by sending (x0 : x1) to
                    (x0*x1 : 1/2*x0^2 + x0*x1 + 3/2*x1^2)
        """
        #it is possible to homogenize the domain and codomain at different coordinates
        if isinstance(n, (tuple, list)):
            ind = tuple(n)
        else:
            ind = (n, n)

        #homogenize the domain and codomain
        A = self.domain().projective_embedding(ind[0]).codomain()
        B = self.codomain().projective_embedding(ind[1]).codomain()
        H = Hom(A, B)
        newvar = A.ambient_space().coordinate_ring().gen(ind[0])

        N = A.ambient_space().dimension_relative()  
        M = B.ambient_space().dimension_relative()  

        #create dictionary for mapping of coordinate rings
        R = self.domain().ambient_space().coordinate_ring()
        S = A.ambient_space().coordinate_ring()
        Rvars = R.gens()
        vars = list(S.gens())
        vars.remove(S.gen(ind[0]))
        D = dict([[Rvars[i],vars[i]] for i in range(N)])

        #clear the denominators if a rational function
        L = [self[i].denominator() for i in range(M)]
        l = [prod(L[:j] + L[j+1:M]) for j in range(M)]
        F = [S(R(self[i].numerator()*l[i]).subs(D)) for i in range(M)]

        #homogenize
        F.insert(ind[1], S(R(prod(L)).subs(D))) #coerce in case l is a constant
        try:
            #remove possible gcd of the polynomials
            g = gcd(F)
            F = [S(f/g) for f in F]
            #remove possible gcd of coefficients
            gc = gcd([f.content() for f in F])
            F = [S(f/gc) for f in F]
        except (AttributeError, ValueError): #no gcd
            pass
        d = max([F[i].degree() for i in range(M+1)])
        F = [F[i].homogenize(str(newvar))*newvar**(d-F[i].degree()) for i in range(M+1)]
        return(H(F))

    def dynatomic_polynomial(self, period):
        r"""
        For a map `f:\mathbb{A}^1 \to \mathbb{A}^1` this function computes 
        the (affine) dynatomic polynomial.

        The dynatomic polynomial is the analog of the cyclotomic polynomial and its roots are the points
        of formal period `n`.

        ALGORITHM:

        Homogenize to a map `f:\mathbb{P}^1 \to \mathbb{P}^1` and compute the dynatomic polynomial there.
        Then, dehomogenize.

        INPUT:

        - ``period`` -- a positive integer or a list/tuple `[m,n]`,
          where `m` is the preperiod and `n` is the period.

        OUTPUT:

        - If possible, a single variable polynomial in the coordinate ring of the polynomial. \
          Otherwise a fraction field element of the coordinate ring of the polynomial.

        EXAMPLES::

            sage: A.<x,y> = AffineSpace(QQ, 2)
            sage: H = Hom(A, A)
            sage: f = H([x^2+y^2, y^2])
            sage: f.dynatomic_polynomial(2)
            Traceback (most recent call last):
            ...
            TypeError: does not make sense in dimension >1

        ::

            sage: A.<x> = AffineSpace(ZZ, 1)
            sage: H = Hom(A, A)
            sage: f = H([(x^2+1)/x])
            sage: f.dynatomic_polynomial(4)
            2*x^12 + 18*x^10 + 57*x^8 + 79*x^6 + 48*x^4 + 12*x^2 + 1

        ::

            sage: A.<x> = AffineSpace(CC, 1)
            sage: H = Hom(A, A)
            sage: f = H([(x^2+1)/(3*x)])
            sage: f.dynatomic_polynomial(3)
            13.0000000000000*x^6 + 117.000000000000*x^4 + 78.0000000000000*x^2 +
            1.00000000000000

        ::

            sage: A.<x> = AffineSpace(QQ, 1)
            sage: H = Hom(A, A)
            sage: f = H([x^2-10/9])
            sage: f.dynatomic_polynomial([2, 1])
            531441*x^4 - 649539*x^2 - 524880

        ::

            sage: A.<x> = AffineSpace(CC, 1)
            sage: H = Hom(A, A)
            sage: f = H([x^2+CC.0])
            sage: f.dynatomic_polynomial(2)
            x^2 + x + 1.00000000000000 + 1.00000000000000*I

        ::

            sage: K.<c> = FunctionField(QQ)
            sage: A.<x> = AffineSpace(K, 1)
            sage: f = Hom(A, A)([x^2 + c])
            sage: f.dynatomic_polynomial(4)
            x^12 + 6*c*x^10 + x^9 + (15*c^2 + 3*c)*x^8 + 4*c*x^7 + (20*c^3 + 12*c^2 + 1)*x^6
            + (6*c^2 + 2*c)*x^5 + (15*c^4 + 18*c^3 + 3*c^2 + 4*c)*x^4 + (4*c^3 + 4*c^2 + 1)*x^3
            + (6*c^5 + 12*c^4 + 6*c^3 + 5*c^2 + c)*x^2 + (c^4 + 2*c^3 + c^2 + 2*c)*x
            + c^6 + 3*c^5 + 3*c^4 + 3*c^3 + 2*c^2 + 1

        ::

            sage: A.<z> = AffineSpace(QQ, 1)
            sage: H = End(A)
            sage: f = H([z^2+3/z+1/7])
            sage: f.dynatomic_polynomial(1).parent()
            Multivariate Polynomial Ring in z over Rational Field
        """
        if self.domain() != self.codomain():
            raise TypeError("must have same domain and codomain to iterate")
        from sage.schemes.affine.affine_space import is_AffineSpace
        if is_AffineSpace(self.domain())==False:
            raise NotImplementedError("not implemented for subschemes")
        if self.domain().dimension_relative()>1:
            raise TypeError("does not make sense in dimension >1")
        F = self.homogenize(1).dynatomic_polynomial(period)
        S = self.domain().coordinate_ring()
        if S(F.denominator()).degree() == 0:
            R = F.parent()
            phi = R.hom([S.gen(0), 1], S)
            return(phi(F))
        else:
            R = F.numerator().parent()
            phi = R.hom([S.gen(0), 1], S)
            return(phi(F.numerator())/phi(F.denominator()))

    def nth_iterate_map(self, n):
        r"""
        This function returns the ``n``-th iterate of the map.

        ALGORITHM:

        Uses a form of successive squaring to reducing computations.

        .. TODO::

            This could be improved.

        INPUT:

        - ``n`` - a positive integer.

        OUTPUT:

        - A map between Affine spaces.

        EXAMPLES::

            sage: A.<x,y> = AffineSpace(ZZ, 2)
            sage: H = Hom(A, A)
            sage: f = H([(x^2-2)/(2*y), y^2-3*x])
            sage: f.nth_iterate_map(2)
            Scheme endomorphism of Affine Space of dimension 2 over Integer Ring
              Defn: Defined on coordinates by sending (x, y) to
                    ((x^4 - 4*x^2 - 8*y^2 + 4)/(8*y^4 - 24*x*y^2), (2*y^5 - 12*x*y^3
            + 18*x^2*y - 3*x^2 + 6)/(2*y))

        ::

            sage: A.<x> = AffineSpace(QQ, 1)
            sage: H = Hom(A, A)
            sage: f = H([(3*x^2-2)/(x)])
            sage: f.nth_iterate_map(3)
            Scheme endomorphism of Affine Space of dimension 1 over Rational Field
              Defn: Defined on coordinates by sending (x) to
                    ((2187*x^8 - 6174*x^6 + 6300*x^4 - 2744*x^2 + 432)/(81*x^7 -
            168*x^5 + 112*x^3 - 24*x))

        ::

            sage: A.<x,y> = AffineSpace(ZZ, 2)
            sage: X = A.subscheme([x-y^2])
            sage: H = Hom(X, X)
            sage: f = H([9*x^2, 3*y])
            sage: f.nth_iterate_map(2)
            Scheme endomorphism of Closed subscheme of Affine Space of dimension 2
            over Integer Ring defined by:
              -y^2 + x
              Defn: Defined on coordinates by sending (x, y) to
                    (729*x^4, 9*y)
        """
        if self.domain() != self.codomain():
            raise TypeError("domain and codomain of function not equal")
        N = self.codomain().ambient_space().dimension_relative()
        F = list(self._polys)
        R = F[0].parent()
        Coord_ring = self.codomain().ambient_space().coordinate_ring()
        D = Integer(n).digits(2)
        PHI = [Coord_ring.gen(i) for i in range(N)]
        for i in range(len(D)):
            T = [F[j] for j in range(N)]
            for k in range(D[i]):
                PHI = [PHI[j](T) for j in range(N)]
            if i != len(D)-1: #avoid extra iterate
                F = [R(F[j](T)) for j in range(N)] #'square'
        H = Hom(self.domain(), self.codomain())
        return(H(PHI))

    def nth_iterate(self, P, n):
        r"""
        Returns the ``n``-th iterate of the point ``P`` by this map.

        INPUT:

        - ``P`` -- a point in the map's domain.

        - ``n`` -- a positive integer.

        OUTPUT:

        - a point in the map's codomain.

        EXAMPLES::

            sage: A.<x,y> = AffineSpace(QQ, 2)
            sage: H = Hom(A, A)
            sage: f = H([(x-2*y^2)/x, 3*x*y])
            sage: f.nth_iterate(A(9, 3), 3)
            (-104975/13123, -9566667)

        ::

            sage: A.<x,y> = AffineSpace(ZZ, 2)
            sage: X = A.subscheme([x-y^2])
            sage: H = Hom(X, X)
            sage: f = H([9*y^2, 3*y])
            sage: f.nth_iterate(X(9, 3), 4)
            (59049, 243)

        ::

            sage: R.<t> = PolynomialRing(QQ)
            sage: A.<x,y> = AffineSpace(FractionField(R), 2)
            sage: H = Hom(A, A)
            sage: f = H([(x-t*y^2)/x, t*x*y])
            sage: f.nth_iterate(A(1, t), 3)
            ((-t^16 + 3*t^13 - 3*t^10 + t^7 + t^5 + t^3 - 1)/(t^5 + t^3 - 1), -t^9 - t^7 + t^4)

        ::

            sage: A.<x,y,z> = AffineSpace(QQ, 3)
            sage: X = A.subscheme([x^2-y^2])
            sage: H = End(X)
            sage: u,v,w = X.coordinate_ring().gens()
            sage: f = H([x^2, y^2, x+y])
            sage: f.nth_iterate_map(2)
            Scheme endomorphism of Closed subscheme of Affine Space of dimension 3
            over Rational Field defined by:
              x^2 - y^2
              Defn: Defined on coordinates by sending (x, y, z) to
                    (x^4, y^4, x^2 + y^2)
        """
        return(P.nth_iterate(self, n))

    def orbit(self, P, n):
        r"""
        Returns the orbit of ``P`` by the map.

        If `n` is an integer it returns `[P,self(P),\ldots,self^n(P)]`.
        If `n` is a list or tuple `n=[m,k]` it returns `[self^m(P),\ldots,self^k(P)]`.

        INPUT:

        - ``P`` -- a point in the map's domain.

        - ``n`` -- a non-negative integer or list or tuple of two non-negative integers.

        OUTPUT:

        - a list of points in the map's codomain.

        EXAMPLES::

            sage: A.<x,y> = AffineSpace(QQ, 2)
            sage: H = Hom(A, A)
            sage: f = H([(x-2*y^2)/x, 3*x*y])
            sage: f.orbit(A(9, 3), 3)
            [(9, 3), (-1, 81), (13123, -243), (-104975/13123, -9566667)]

        ::

            sage: A.<x> = AffineSpace(QQ, 1)
            sage: H = Hom(A, A)
            sage: f = H([(x-2)/x])
            sage: f.orbit(A(1/2), [1, 3])
            [(-3), (5/3), (-1/5)]

        ::

            sage: A.<x,y> = AffineSpace(ZZ, 2)
            sage: X = A.subscheme([x-y^2])
            sage: H = Hom(X, X)
            sage: f = H([9*y^2, 3*y])
            sage: f.orbit(X(9, 3), (0, 4))
            [(9, 3), (81, 9), (729, 27), (6561, 81), (59049, 243)]

        ::

            sage: R.<t> = PolynomialRing(QQ)
            sage: A.<x,y> = AffineSpace(FractionField(R), 2)
            sage: H = Hom(A, A)
            sage: f = H([(x-t*y^2)/x, t*x*y])
            sage: f.orbit(A(1, t), 3)
            [(1, t), (-t^3 + 1, t^2), ((-t^5 - t^3 + 1)/(-t^3 + 1), -t^6 + t^3),
            ((-t^16 + 3*t^13 - 3*t^10 + t^7 + t^5 + t^3 - 1)/(t^5 + t^3 - 1), -t^9 -
            t^7 + t^4)]
        """
        return(P.orbit(self, n))

    def global_height(self, prec=None):
        r"""
        Returns the maximum of the heights of the coefficients in any
        of the coordinate functions of the affine morphism.

        INPUT:

        - ``prec`` -- desired floating point precision (default:
          default RealField precision).

        OUTPUT: A real number.

        EXAMPLES::

            sage: A.<x> = AffineSpace(QQ, 1)
            sage: H = Hom(A, A)
            sage: f = H([1/1331*x^2 + 4000]);
            sage: f.global_height()
            8.29404964010203

        ::

            sage: R.<x> = PolynomialRing(QQ)
            sage: k.<w> = NumberField(x^2 + 5)
            sage: A.<x,y> = AffineSpace(k, 2)
            sage: H = Hom(A, A)
            sage: f = H([13*w*x^2 + 4*y, 1/w*y^2]);
            sage: f.global_height(prec=100)
            3.3696683136785869233538671082

        ::

            sage: A.<x> = AffineSpace(ZZ, 1)
            sage: H = Hom(A, A)
            sage: f = H([7*x^2 + 1513]);
            sage: f.global_height()
            7.32184971378836
        """
        H=0
        for i in range(self.domain().ambient_space().dimension_relative()):
            C = self[i].coefficients()
            if C == []: #to deal with the case self[i]=0
                h=0
            else:
                h = max([c.global_height(prec) for c in C])
            H = max(H,h)
        return(H)

    def jacobian (self):
        r"""
        Returns the Jacobian matrix of partial derivitive of this map.

        The `(i, j)` entry of the Jacobian matrix is the partial derivative
        `diff(functions[i], variables[j])`.

        OUTPUT:

        - matrix with coordinates in the coordinate ring of the map.

        EXAMPLES::

            sage: A.<z> = AffineSpace(QQ, 1)
            sage: H = End(A)
            sage: f = H([z^2 - 3/4])
            sage: f.jacobian()
            [2*z]

        ::

            sage: A.<x,y> = AffineSpace(QQ, 2)
            sage: H = End(A)
            sage: f = H([x^3 - 25*x + 12*y, 5*y^2*x - 53*y + 24])
            sage: f.jacobian()
            [ 3*x^2 - 25          12]
            [      5*y^2 10*x*y - 53]

        ::

            sage: A.<x,y> = AffineSpace(ZZ, 2)
            sage: H = End(A)
            sage: f = H([(x^2 - x*y)/(1+y), (5+y)/(2+x)])
            sage: f.jacobian()
            [         (2*x - y)/(y + 1) (-x^2 - x)/(y^2 + 2*y + 1)]
            [  (-y - 5)/(x^2 + 4*x + 4)                  1/(x + 2)]
        """
        try:
            return self.__jacobian
        except AttributeError:
            pass
        self.__jacobian = jacobian(list(self),self.domain().ambient_space().gens())
        return self.__jacobian

    def multiplier(self, P, n, check=True):
        r"""
        Returns the multiplier of the point ``P`` of period ``n`` by the map.

        The map must be an endomorphism.

        INPUT:

        - ``P`` - a point on domain of the map.

        - ``n`` - a positive integer, the period of ``P``.

        - ``check`` -- verify that ``P`` has period ``n``, Default:True.

        OUTPUT:

        - a square matrix of size ``self.codomain().dimension_relative()`` in
          the ``base_ring`` of the map.

        EXAMPLES::

            sage: P.<x,y> = AffineSpace(QQ, 2)
            sage: H = End(P)
            sage: f = H([x^2, y^2])
            sage: f.multiplier(P([1, 1]), 1)
            [2 0]
            [0 2]

        ::

            sage: P.<x,y,z> = AffineSpace(QQ, 3)
            sage: H = End(P)
            sage: f = H([x, y^2, z^2 - y])
            sage: f.multiplier(P([1/2, 1, 0]), 2)
            [1 0 0]
            [0 4 0]
            [0 0 0]

        ::

            sage: P.<x> = AffineSpace(CC, 1)
            sage: H = End(P)
            sage: f = H([x^2 + 1/2])
            sage: f.multiplier(P([0.5 + 0.5*I]), 1)
            [1.00000000000000 + 1.00000000000000*I]

        ::

            sage: R.<t> = PolynomialRing(CC, 1)
            sage: P.<x> = AffineSpace(R, 1)
            sage: H = End(P)
            sage: f = H([x^2 - t^2 + t])
            sage: f.multiplier(P([-t + 1]), 1)
            [(-2.00000000000000)*t + 2.00000000000000]

        ::

            sage: P.<x,y> = AffineSpace(QQ, 2)
            sage: X = P.subscheme([x^2-y^2])
            sage: H = End(X)
            sage: f = H([x^2, y^2])
            sage: f.multiplier(X([1, 1]), 1)
            [2 0]
            [0 2]
        """
        if not self.is_endomorphism():
            raise TypeError("must be an endomorphism")
        if check:
            if self.nth_iterate(P, n) != P:
                raise ValueError("%s is not periodic of period %s" % (P, n))
            if n < 1:
                raise ValueError("period must be a positive integer")
        N = self.domain().ambient_space().dimension_relative()
        l = identity_matrix(FractionField(self.codomain().base_ring()), N, N)
        Q = P
        J = self.jacobian()
        for i in range(0, n):
            R = self(Q)
            l = J(tuple(Q))*l #chain rule matrix multiplication
            Q = R
        return l

class SchemeMorphism_polynomial_affine_space_field(SchemeMorphism_polynomial_affine_space):

    @cached_method
    def weil_restriction(self):
        r"""
        Compute the Weil restriction of this morphism over some extension field.

        If the field is a finite field, then this computes
        the Weil restriction to the prime subfield.

        A Weil restriction of scalars - denoted `Res_{L/k}` - is a
        functor which, for any finite extension of fields `L/k` and
        any algebraic variety `X` over `L`, produces another
        corresponding variety `Res_{L/k}(X)`, defined over `k`. It is
        useful for reducing questions about varieties over large
        fields to questions about more complicated varieties over
        smaller fields. Since it is a functor it also applied to morphisms.
        In particular, the functor applied to a morphism gives the equivalent
        morphism from the Weil restriction of the domain to the Weil restriction
        of the codomain.

        OUTPUT: Scheme morphism on the Weil restrictions of the domain
                and codomain of the map.

        EXAMPLES::

            sage: K.<v> = QuadraticField(5)
            sage: A.<x,y> = AffineSpace(K, 2)
            sage: H = End(A)
            sage: f = H([x^2-y^2, y^2])
            sage: f.weil_restriction()
            Scheme endomorphism of Affine Space of dimension 4 over Rational Field
              Defn: Defined on coordinates by sending (z0, z1, z2, z3) to
                    (z0^2 + 5*z1^2 - z2^2 - 5*z3^2, 2*z0*z1 - 2*z2*z3, z2^2 + 5*z3^2, 2*z2*z3)

        ::

            sage: K.<v> = QuadraticField(5)
            sage: PS.<x,y> = AffineSpace(K, 2)
            sage: H = Hom(PS, PS)
            sage: f = H([x, y])
            sage: F = f.weil_restriction()
            sage: P = PS(2, 1)
            sage: Q = P.weil_restriction()
            sage: f(P).weil_restriction() == F(Q)
            True
        """
        if any([isinstance(f,FractionFieldElement) for f in self]):
            raise TypeError("coordinate functions must be polynomials")

        DS = self.domain()
        R = DS.coordinate_ring()
        #using the Weil restriction on ideal generators to not duplicate code
        result = R.ideal(self._polys).weil_restriction().gens()
        H = Hom(DS.weil_restriction(), self.codomain().weil_restriction())

        return(H(result))

class SchemeMorphism_polynomial_affine_space_finite_field(SchemeMorphism_polynomial_affine_space_field):

    def orbit_structure(self, P):
        r"""
        Every point is preperiodic over a finite field.

        This function returns the pair `[m,n]` where `m` is the
        preperiod and `n` is the period of the point ``P`` by this map.

        INPUT:

        - ``P`` -- a point in the map's domain.

        OUTPUT:

        - a list `[m, n]` of integers.

        EXAMPLES::

            sage: A.<x,y> = AffineSpace(GF(13), 2)
            sage: H = Hom(A, A)
            sage: f = H([x^2 - 1, y^2])
            sage: f.orbit_structure(A(2, 3))
            [1, 6]

        ::

            sage: A.<x,y,z> = AffineSpace(GF(49, 't'), 3)
            sage: H = Hom(A, A)
            sage: f = H([x^2 - z, x - y + z, y^2 - x^2])
            sage: f.orbit_structure(A(1, 1, 2))
            [7, 6]
        """
        return(P.orbit_structure(self))

    def _fast_eval(self, x):
        """
        Evaluate affine morphism at point described by ``x``.

        EXAMPLES::

            sage: P.<x,y,z> = AffineSpace(GF(7), 3)
            sage: H = Hom(P, P)
            sage: f = H([x^2+y^2,y^2, z^2 + y*z])
            sage: f._fast_eval([1, 1, 1])
            [2, 1, 2]

        ::

            sage: P.<x,y,z> = AffineSpace(GF(19), 3)
            sage: H = Hom(P, P)
            sage: f = H([x/(y+1), y, (z^2 + y^2)/(x^2 + 1)])
            sage: f._fast_eval([2, 1, 3])
            [1, 1, 2]
        """

        R = self.domain().ambient_space().coordinate_ring()

        P=[]
        for i in range(len(self._fastpolys[0])):
            r = self._fastpolys[0][i](*x)
            if self._fastpolys[1][i] is R.one():
                if self._is_prime_finite_field:
                    p = self.base_ring().characteristic()
                    r = Integer(r) % p
                P.append(r)
            else:
                s = self._fastpolys[1][i](*x)
                if self._is_prime_finite_field:
                    p = self.base_ring().characteristic()
                    r = Integer(r) % p
                    s = Integer(s) % p
                P.append(r/s)
        return P

    def cyclegraph(self):
        r"""
        Returns digraph of all orbits of this morphism mod `p`.

        For subschemes, only points on the subscheme whose
        image are also on the subscheme are in the digraph.

        OUTPUT: A digraph.

        EXAMPLES::

            sage: P.<x,y> = AffineSpace(GF(5), 2)
            sage: H = Hom(P, P)
            sage: f = H([x^2-y, x*y+1])
            sage: f.cyclegraph()
            Looped digraph on 25 vertices

        ::

            sage: P.<x> = AffineSpace(GF(3^3, 't'), 1)
            sage: H = Hom(P, P)
            sage: f = H([x^2-1])
            sage: f.cyclegraph()
            Looped digraph on 27 vertices

        ::

            sage: P.<x,y> = AffineSpace(GF(7), 2)
            sage: X = P.subscheme(x-y)
            sage: H = Hom(X, X)
            sage: f = H([x^2, y^2])
            sage: f.cyclegraph()
            Looped digraph on 7 vertices
        """
        if self.domain() != self.codomain():
            raise NotImplementedError("domain and codomain must be equal")
        V = []
        E = []
        from sage.schemes.affine.affine_space import is_AffineSpace
        if is_AffineSpace(self.domain()) == True:
            for P in self.domain():
                V.append(str(P))
                Q = self(P)
                E.append([str(Q)])
        else:
            X = self.domain()
            for P in X.ambient_space():
                try:
                    XP = X.point(P)
                    V.append(str(XP))
                    Q = self(XP)
                    E.append([str(Q)])
                except TypeError:  # not on the scheme
                    pass
        from sage.graphs.digraph import DiGraph
        g = DiGraph(dict(zip(V, E)), loops=True)
        return g
