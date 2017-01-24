# -*- coding: utf-8 -*-
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

- Dillon Rose (2014-01):  Speed enhancements

- Ben Hutz (2015-11): iteration of subschemes


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
from __future__ import print_function
from __future__ import absolute_import

from sage.calculus.functions import jacobian
from sage.categories.number_fields import NumberFields
from sage.categories.homset import Hom, End
from sage.combinat.sf.sf import SymmetricFunctions
from sage.functions.all import sqrt
from sage.libs.pari.all import PariError
from sage.matrix.constructor import matrix, identity_matrix
from sage.misc.all import prod
from sage.misc.cachefunc import cached_method
from sage.misc.misc import subsets
from sage.misc.mrange import xmrange
from sage.modules.free_module_element import vector
from sage.rings.all import Integer, CIF
from sage.arith.all import gcd, lcm, next_prime, binomial, primes, moebius
from sage.rings.complex_field import ComplexField_class,ComplexField
from sage.rings.complex_interval_field import ComplexIntervalField_class
from sage.categories.finite_fields import FiniteFields
from sage.rings.finite_rings.finite_field_constructor import GF, is_PrimeFiniteField
from sage.rings.finite_rings.integer_mod_ring import Zmod
from sage.rings.fraction_field import FractionField
from sage.rings.fraction_field_element import is_FractionFieldElement, FractionFieldElement
from sage.rings.integer_ring import ZZ
from sage.rings.number_field.order import is_NumberFieldOrder
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.qqbar import QQbar, number_field_elements_from_algebraics
from sage.rings.quotient_ring import QuotientRing_generic
from sage.rings.qqbar import QQbar
from sage.rings.rational_field import QQ
from sage.rings.real_mpfr import RealField_class,RealField
from sage.rings.real_mpfi import RealIntervalField_class
from sage.schemes.generic.morphism import SchemeMorphism_polynomial
from sage.symbolic.constants import e
from copy import copy
from sage.parallel.ncpus import ncpus
from sage.parallel.use_fork import p_iter_fork
from sage.ext.fast_callable import fast_callable
from sage.misc.lazy_attribute import lazy_attribute
from sage.schemes.projective.projective_morphism_helper import _fast_possible_periods
import sys
from sage.sets.set import Set
from sage.combinat.permutation import Arrangements
from sage.combinat.subset import Subsets
from sage.categories.number_fields import NumberFields
_NumberFields = NumberFields()
from sage.rings.number_field.number_field import NumberField

class SchemeMorphism_polynomial_projective_space(SchemeMorphism_polynomial):
    r"""
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
          From: Projective Plane Curve over Rational Field defined by x^3 + y^3 + 60*z^3
          To:   Projective Plane Curve over Rational Field defined by -x^3 + y^2*z + 6400/3*z^3
          Defn: Defined on coordinates by sending (x : y : z) to
                (z : x - y : -1/80*x - 1/80*y)

    A more complicated example::

        sage: P2.<x,y,z> = ProjectiveSpace(2, QQ)
        sage: P1 = P2.subscheme(x-y)
        sage: H12 = P1.Hom(P2)
        sage: H12([x^2, x*z, z^2])
        Scheme morphism:
          From: Closed subscheme of Projective Space of dimension 2 over Rational Field defined by:
          x - y
          To:   Projective Space of dimension 2 over Rational Field
          Defn: Defined on coordinates by sending (x : y : z) to
              (x^2 : x*z : z^2)

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

    We can also compute the forward image of subschemes through
    elimination. In particular, let `X = V(h_1,\ldots, h_t)` and define the ideal
    `I = (h_1,\ldots,h_t,y_0-f_0(\bar{x}), \ldots, y_n-f_n(\bar{x}))`.
    Then the elimination ideal `I_{n+1} = I \cap K[y_0,\ldots,y_n]` is a homogeneous
    ideal and `f(X) = V(I_{n+1})`::

        sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
        sage: H = End(P)
        sage: f = H([(x-2*y)^2, (x-2*z)^2, x^2])
        sage: X = P.subscheme(y-z)
        sage: f(f(f(X)))
        Closed subscheme of Projective Space of dimension 2 over Rational Field
        defined by:
          y - z

    ::

        sage: P.<x,y,z,w> = ProjectiveSpace(QQ, 3)
        sage: H = End(P)
        sage: f = H([(x-2*y)^2, (x-2*z)^2, (x-2*w)^2, x^2])
        sage: f(P.subscheme([x,y,z]))
        Closed subscheme of Projective Space of dimension 3 over Rational Field
        defined by:
          w,
          y,
          x
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

        ::

            sage: R.<t> = PolynomialRing(QQ)
            sage: P.<x,y,z> = ProjectiveSpace(R, 2)
            sage: X = P.subscheme([x])
            sage: H = End(X)
            sage: H([x^2, t*y^2, x*z])
            Scheme endomorphism of Closed subscheme of Projective Space of dimension
            2 over Univariate Polynomial Ring in t over Rational Field defined by:
              x
              Defn: Defined on coordinates by sending (x : y : z) to
                    (x^2 : t*y^2 : x*z)

        When elements of the quotient ring is used, they are reduced::

            sage: P.<x,y,z> = ProjectiveSpace(CC, 2)
            sage: X = P.subscheme([x-y])
            sage: u,v,w = X.coordinate_ring().gens()
            sage: H = End(X)
            sage: H([u^2, v^2, w*u])
            Scheme endomorphism of Closed subscheme of Projective Space of dimension
            2 over Complex Field with 53 bits of precision defined by:
              x - y
              Defn: Defined on coordinates by sending (x : y : z) to
                    (y^2 : y^2 : y*z)
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
                raise  ValueError("polys (=%s) must be homogeneous" % polys)
            degs = [f.degree() for f in polys]
            if not all([d == degs[0] for d in degs[1:]]):
                raise ValueError("polys (=%s) must be of the same degree" % polys)
        self._is_prime_finite_field = is_PrimeFiniteField(polys[0].base_ring())

    def __call__(self, x, check=True):
        """
        Compute the forward image of the point or subscheme ``x`` by this map.

        For subschemes, the forward image is computed through elimination.
        In particular, let `X = V(h_1,\ldots, h_t)` and define the ideal
        `I = (h_1,\ldots,h_t,y_0-f_0(\bar{x}), \ldots, y_n-f_n(\bar{x}))`.
        Then the elimination ideal `I_{n+1} = I \cap K[y_0,\ldots,y_n]` is a homogeneous
        ideal and `self(X) = V(I_{n+1})`.

        The input boolean ``check`` can be set to false when fast iteration of
        points is desired. It bypasses all input checking and passes ``x`` straight
        to the fast evaluation of points function.

        INPUT:

        - ``x`` - a point or subscheme in domain of this map.

        - ``check`` - Boolean - if `False` assume that ``x`` is a point.

        EXAMPLES::

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: H = End(P)
            sage: f = H([x^2+y^2, y^2, z^2 + y*z])
            sage: f(P([1,1,1]))
            (1 : 1/2 : 1)

        ::

            sage: PS.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: P1.<u,v> = ProjectiveSpace(QQ,1)
            sage: H = End(P1)
            sage: f = H([u^2, v^2])
            sage: f(PS([0,1,1]))
            Traceback (most recent call last):
            ...
            TypeError: (0 : 1 : 1) fails to convert into the map's domain Projective Space of
            dimension 1 over Rational Field, but a `pushforward` method is not properly implemented

        ::

            sage: PS.<x,y> = ProjectiveSpace(QQ, 1)
            sage: P1.<u,v> = ProjectiveSpace(QQ, 1)
            sage: H = End(P1)
            sage: f = H([u^2, v^2])
            sage: f([0,1])
            (0 : 1)
            sage: f(PS([0,1]))
            (0 : 1)

        ::

            sage: PS.<x,y,z,w> = ProjectiveSpace(QQ, 3)
            sage: H = End(PS)
            sage: f = H([y^2, x^2, w^2, z^2])
            sage: X = PS.subscheme([z^2+y*w])
            sage: f(X)
            Closed subscheme of Projective Space of dimension 3 over Rational Field
            defined by:
              x*z - w^2

        ::

            sage: PS.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: P1.<u,v> = ProjectiveSpace(ZZ, 1)
            sage: H = End(PS)
            sage: f = H([x^2, y^2, z^2])
            sage: X = P1.subscheme([u-v])
            sage: f(X)
            Traceback (most recent call last):
            ...
            TypeError: subscheme must be in ambient space of domain of map

        ::

            sage: PS.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: P1.<u,v> = ProjectiveSpace(ZZ, 1)
            sage: H = End(P1)
            sage: f = H([u^2, v^2])
            sage: f([u-v])
            Closed subscheme of Projective Space of dimension 1 over Integer Ring defined by:
              u - v
            sage: X = PS.subscheme([x-z])
            sage: f([x-z])
            Traceback (most recent call last):
            ...
            TypeError: [x - z] fails to convert into the map's domain Projective Space of
            dimension 1 over Integer Ring, but a `pushforward` method is not properly implemented
        """
        from sage.schemes.projective.projective_point import SchemeMorphism_point_projective_ring
        if check:
            from sage.schemes.generic.algebraic_scheme import AlgebraicScheme_subscheme_projective
            if isinstance(x, SchemeMorphism_point_projective_ring):
                if self.domain() != x.codomain():
                    try:
                        x = self.domain()(x)
                    except (TypeError, NotImplementedError):
                        raise TypeError("%s fails to convert into the map's domain %s, but a `pushforward` method is not properly implemented"%(x, self.domain()))
                #else pass it onto the eval below
            elif isinstance(x, AlgebraicScheme_subscheme_projective):
                return x._forward_image(self) #call subscheme eval
            else: #not a projective point or subscheme
                try:
                    x = self.domain()(x)
                except (TypeError, NotImplementedError):
                    try:
                        x = self.domain().subscheme(x)
                        return x._forward_image(self) #call subscheme eval
                    except (TypeError, NotImplementedError):
                        raise TypeError("%s fails to convert into the map's domain %s, but a `pushforward` method is not properly implemented"%(x, self.domain()))

        # Passes the array of args to _fast_eval
        P = self._fast_eval(x._coords)
        return self.codomain().point(P, check)

    @lazy_attribute
    def _fastpolys(self):
        """
        Lazy attribute for fast_callable polynomials for this map.

        EXAMPLES::

            sage: P.<x,y>=ProjectiveSpace(QQ,1)
            sage: H = Hom(P,P)
            sage: f = H([x^2+y^2,y^2])
            sage: [g.op_list() for g in f._fastpolys]
            [[('load_const', 0), ('load_const', 1), ('load_arg', 1), ('ipow', 2), 'mul', 'add', ('load_const', 1), ('load_arg', 0), ('ipow', 2), 'mul', 'add', 'return'], [('load_const', 0), ('load_const', 1), ('load_arg', 1), ('ipow', 2), 'mul', 'add', 'return']]
        """
        polys = self._polys

        fastpolys = []
        for poly in polys:
            # These tests are in place because the float and integer domain evaluate
            # faster than using the base_ring
            if self._is_prime_finite_field:
                prime = polys[0].base_ring().characteristic()
                degree = polys[0].degree()
                coefficients = poly.coefficients()
                height = max(abs(c.lift()) for c in coefficients)
                num_terms = len(coefficients)
                largest_value = num_terms * height * (prime - 1) ** degree
                # If the calculations will not overflow the float data type use domain float
                # Else use domain integer
                if largest_value < (2 ** sys.float_info.mant_dig):
                    fastpolys.append(fast_callable(poly, domain=float))
                else:
                    fastpolys.append(fast_callable(poly, domain=ZZ))
            else:
                fastpolys.append(fast_callable(poly, domain=poly.base_ring()))
        return fastpolys

    def _fast_eval(self, x):
        """
        Evaluate projective morphism at point described by ``x``.

        EXAMPLES::

            sage: P.<x,y,z> = ProjectiveSpace(QQ,2)
            sage: H = Hom(P,P)
            sage: f = H([x^2+y^2, y^2, z^2 + y*z])
            sage: f._fast_eval([1,1,1])
            [2, 1, 2]

            ::

            sage: T.<z> = LaurentSeriesRing(ZZ)
            sage: P.<x,y> = ProjectiveSpace(T,1)
            sage: H = End(P)
            sage: f = H([x^2+x*y, y^2])
            sage: Q = P(z,1)
            sage: f._fast_eval(list(Q))
            [z + z^2, 1]

            ::

            sage: T.<z> = PolynomialRing(CC)
            sage: I = T.ideal(z^3)
            sage: P.<x,y> = ProjectiveSpace(T.quotient_ring(I),1)
            sage: H = End(P)
            sage: f = H([x^2+x*y, y^2])
            sage: Q = P(z^2, 1)
            sage: f._fast_eval(list(Q))
            [zbar^2, 1.00000000000000]

            ::

            sage: T.<z> = LaurentSeriesRing(CC)
            sage: R.<t> = PolynomialRing(T)
            sage: P.<x,y> = ProjectiveSpace(R,1)
            sage: H = End(P)
            sage: f = H([x^2+x*y, y^2])
            sage: Q = P(t^2, z)
            sage: f._fast_eval(list(Q))
            [t^4 + z*t^2, z^2]
        """
        P = [f(*x) for f in self._fastpolys]
        return P

    def __eq__(self, right):
        """
        Tests the equality of two projective morphisms.

        INPUT:

        - ``right`` - a map on projective space.

        OUTPUT:

        - Boolean - True if ``self`` and ``right`` define the same projective map. False otherwise.

        EXAMPLES::

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: H = Hom(P,P)
            sage: f = H([x^2 - 2*x*y + z*x, z^2 -y^2 , 5*z*y])
            sage: g = H([x^2, y^2, z^2])
            sage: f == g
            False

            ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: P2.<u,v> = ProjectiveSpace(CC, 1)
            sage: H = End(P)
            sage: H2 = End(P2)
            sage: f = H([x^2 - 2*x*y, y^2])
            sage: g = H2([u^2 - 2*u*v, v^2])
            sage: f == g
            False

            ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: H = End(P)
            sage: f = H([x^2 - 2*x*y, y^2])
            sage: g = H([x^2*y - 2*x*y^2, y^3])
            sage: f == g
            True
        """
        if not isinstance(right, SchemeMorphism_polynomial):
            return False
        if self.parent() != right.parent():
            return False
        n = len(self._polys)
        return all([self[i]*right[j] == self[j]*right[i] for i in range(0, n) for j in range(i+1, n)])

    def __ne__(self, right):
        """
        Tests the inequality of two projective morphisms.

        INPUT:

        - ``right`` -- a map on projective space.

        OUTPUT:

        - Boolean -- True if ``self`` and ``right`` define different projective maps. False otherwise.

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: H = Hom(P,P)
            sage: f = H([x^3 - 2*x^2*y , 5*x*y^2])
            sage: g = f.change_ring(GF(7))
            sage: f != g
            True

            ::

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: H = Hom(P, P)
            sage: f = H([x^2 - 2*x*y + z*x, z^2 -y^2 , 5*z*y])
            sage: f != f
            False
        """
        if not isinstance(right, SchemeMorphism_polynomial):
            return True
        if self.parent() != right.parent():
            return True
        n = len(self._polys)
        for i in range(0, n):
            for j in range(i + 1, n):
                if self._polys[i] * right._polys[j] != self._polys[j] * right._polys[i]:
                    return True
        return False

    def scale_by(self, t):
        """
        Scales each coordinate by a factor of ``t``.

        A ``TypeError`` occurs if the point is not in the coordinate_ring
        of the parent after scaling.

        INPUT:

        - ``t`` -- a ring element.

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
        R = self.domain().coordinate_ring()
        if isinstance(R, QuotientRing_generic):
            phi = R._internal_coerce_map_from(self.domain().ambient_space().coordinate_ring())
            for i in range(self.codomain().ambient_space().dimension_relative() + 1):
                self._polys[i] = phi(self._polys[i] * t).lift()
        else:
            for i in range(self.codomain().ambient_space().dimension_relative() + 1):
                self._polys[i] = R(self._polys[i] * t)

    def normalize_coordinates(self):
        """
        Scales by 1/gcd of the coordinate functions.

        Also, scales to clear any denominators from the coefficients. This is done in place.

        OUTPUT:

        - None.

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = Hom(P,P)
            sage: f = H([5/4*x^3, 5*x*y^2])
            sage: f.normalize_coordinates(); f
            Scheme endomorphism of Projective Space of dimension 1 over Rational
            Field
              Defn: Defined on coordinates by sending (x : y) to
                    (x^2 : 4*y^2)

        ::

            sage: P.<x,y,z> = ProjectiveSpace(GF(7),2)
            sage: X = P.subscheme(x^2-y^2)
            sage: H = Hom(X,X)
            sage: f = H([x^3+x*y^2, x*y^2, x*z^2])
            sage: f.normalize_coordinates(); f
            Scheme endomorphism of Closed subscheme of Projective Space of dimension
            2 over Finite Field of size 7 defined by:
              x^2 - y^2
              Defn: Defined on coordinates by sending (x : y : z) to
                    (2*y^2 : y^2 : z^2)

        ::

            sage: R.<a,b> = QQ[]
            sage: P.<x,y,z> = ProjectiveSpace(R, 2)
            sage: H = End(P)
            sage: f = H([a*(x*z+y^2)*x^2, a*b*(x*z+y^2)*y^2, a*(x*z+y^2)*z^2])
            sage: f.normalize_coordinates(); f
            Scheme endomorphism of Projective Space of dimension 2 over Multivariate
            Polynomial Ring in a, b over Rational Field
              Defn: Defined on coordinates by sending (x : y : z) to
                    (x^2 : b*y^2 : z^2)

        .. NOTE:: gcd raises an error if the base_ring does not support gcds.
        """
        GCD = gcd(self[0], self[1])
        index = 2
        if self[0].lc() > 0 or self[1].lc() > 0:
            neg = 0
        else:
            neg = 1
        N = self.codomain().ambient_space().dimension_relative() + 1
        while GCD != 1 and index < N:
            if self[index].lc() > 0:
                neg = 0
            GCD = gcd(GCD, self[index])
            index += +1

        if GCD != 1:
            R = self.domain().base_ring()
            if neg == 1:
                self.scale_by(R(-1) / GCD)
            else:
                self.scale_by(R(1) / GCD)
        else:
            if neg == 1:
                self.scale_by(-1)

        #clears any denominators from the coefficients
        LCM = lcm([self[i].denominator() for i in range(N)])
        self.scale_by(LCM)

        #scales by 1/gcd of the coefficients.
        GCD = gcd([self[i].content() for i in range(N)])
        if GCD != 1:
            self.scale_by(1 / GCD)

    def dynatomic_polynomial(self, period):
        r"""
        For a map `f:\mathbb{P}^1 \to \mathbb{P}^1` this function computes the dynatomic polynomial.

        The dynatomic polynomial is the analog of the cyclotomic
        polynomial and its roots are the points of formal period `period`. If possible the division is
        done in the coordinate ring of this map and a polynomial is returned. In rings where that is not possible,
        a FractionField element will be returned. In certain cases, when the conversion back to a polynomial
        fails, a SymbolRing element will be returned.

        ALGORITHM:

        For a positive integer `n`, let `[F_n,G_n]` be the coordinates of the `nth` iterate of `f`.
        Then construct

        .. MATH::

            \Phi^{\ast}_n(f)(x,y) = \sum_{d \mid n} (yF_d(x,y) - xG_d(x,y))^{\mu(n/d)}

        where `\mu` is the Möbius function.

        For a pair `[m,n]`, let `f^m = [F_m,G_m]`. Compute

        .. MATH::

            \Phi^{\ast}_{m,n}(f)(x,y) = \Phi^{\ast}_n(f)(F_m,G_m)/\Phi^{\ast}_n(f)(F_{m-1},G_{m-1})

        REFERENCES:

        .. [Hutz] \B. Hutz. Determination of all rational preperiodic points
           for morphisms of PN. Mathematics of Computation, 84:291 (2015), 289-308.

        .. [MoPa] \P. Morton and P. Patel. The Galois theory of periodic points
           of polynomial maps. Proc. London Math. Soc., 68 (1994), 225-263.

        INPUT:

        - ``period`` -- a positive integer or a list/tuple `[m,n]` where
          `m` is the preperiod and `n` is the period.

        OUTPUT:

        - If possible, a two variable polynomial in the coordinate ring of this map.
          Otherwise a fraction field element of the coordinate ring of this map. Or,
          a Symbolic Ring element.

        .. TODO::

            - Do the division when the base ring is p-adic so that the output is a polynomial.

            - Convert back to a polynomial when the base ring is a function field (not over `\QQ` or `F_p`)

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = Hom(P,P)
            sage: f = H([x^2+y^2, y^2])
            sage: f.dynatomic_polynomial(2)
            x^2 + x*y + 2*y^2

        ::

            sage: P.<x,y> = ProjectiveSpace(ZZ,1)
            sage: H = Hom(P,P)
            sage: f = H([x^2+y^2, x*y])
            sage: f.dynatomic_polynomial(4)
            2*x^12 + 18*x^10*y^2 + 57*x^8*y^4 + 79*x^6*y^6 + 48*x^4*y^8 + 12*x^2*y^10 + y^12

        ::

            sage: P.<x,y> = ProjectiveSpace(CC,1)
            sage: H = Hom(P,P)
            sage: f = H([x^2+y^2, 3*x*y])
            sage: f.dynatomic_polynomial(3)
            13.0000000000000*x^6 + 117.000000000000*x^4*y^2 +
            78.0000000000000*x^2*y^4 + y^6

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = Hom(P,P)
            sage: f = H([x^2-10/9*y^2, y^2])
            sage: f.dynatomic_polynomial([2,1])
            x^4*y^2 - 11/9*x^2*y^4 - 80/81*y^6

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = Hom(P,P)
            sage: f = H([x^2-29/16*y^2, y^2])
            sage: f.dynatomic_polynomial([2,3])
            x^12 - 95/8*x^10*y^2 + 13799/256*x^8*y^4 - 119953/1024*x^6*y^6 +
            8198847/65536*x^4*y^8 - 31492431/524288*x^2*y^10 +
            172692729/16777216*y^12

        ::

            sage: P.<x,y> = ProjectiveSpace(ZZ,1)
            sage: H = Hom(P,P)
            sage: f = H([x^2-y^2, y^2])
            sage: f.dynatomic_polynomial([1,2])
            x^2 - x*y

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = Hom(P,P)
            sage: f = H([x^3-y^3, 3*x*y^2])
            sage: f.dynatomic_polynomial([0,4])==f.dynatomic_polynomial(4)
            True

        ::

            sage: P.<x,y,z> = ProjectiveSpace(QQ,2)
            sage: H = Hom(P,P)
            sage: f = H([x^2+y^2, x*y, z^2])
            sage: f.dynatomic_polynomial(2)
            Traceback (most recent call last):
            ...
            TypeError: does not make sense in dimension >1

        ::

            sage: P.<x,y> = ProjectiveSpace(Qp(5),1)
            sage: H = Hom(P,P)
            sage: f = H([x^2+y^2, y^2])
            sage: f.dynatomic_polynomial(2)
            (x^4*y + (2 + O(5^20))*x^2*y^3 - x*y^4 + (2 + O(5^20))*y^5)/(x^2*y -
            x*y^2 + y^3)

        ::

            sage: L.<t> = PolynomialRing(QQ)
            sage: P.<x,y> = ProjectiveSpace(L,1)
            sage: H = Hom(P,P)
            sage: f = H([x^2+t*y^2, y^2])
            sage: f.dynatomic_polynomial(2)
            x^2 + x*y + (t + 1)*y^2

        ::

            sage: K.<c> = PolynomialRing(ZZ)
            sage: P.<x,y> = ProjectiveSpace(K,1)
            sage: H = Hom(P,P)
            sage: f = H([x^2+ c*y^2, y^2])
            sage: f.dynatomic_polynomial([1, 2])
            x^2 - x*y + (c + 1)*y^2

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = Hom(P,P)
            sage: f = H([x^2+y^2, y^2])
            sage: f.dynatomic_polynomial(2)
            x^2 + x*y + 2*y^2
            sage: R.<X> = PolynomialRing(QQ)
            sage: K.<c> = NumberField(X^2 + X + 2)
            sage: PP = P.change_ring(K)
            sage: ff = f.change_ring(K)
            sage: p = PP((c, 1))
            sage: ff(ff(p)) == p
            True

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = Hom(P,P)
            sage: f = H([x^2+y^2, x*y])
            sage: f.dynatomic_polynomial([2, 2])
            x^4 + 4*x^2*y^2 + y^4
            sage: R.<X> = PolynomialRing(QQ)
            sage: K.<c> = NumberField(X^4 + 4*X^2 + 1)
            sage: PP = P.change_ring(K)
            sage: ff = f.change_ring(K)
            sage: p = PP((c, 1))
            sage: ff.nth_iterate(p, 4) == ff.nth_iterate(p, 2)
            True

        ::

            sage: P.<x,y> = ProjectiveSpace(CC, 1)
            sage: H = Hom(P,P)
            sage: f = H([x^2-CC.0/3*y^2, y^2])
            sage: f.dynatomic_polynomial(2)
            0.666666666666667*x^2 + 0.333333333333333*y^2

        ::

            sage: L.<t> = PolynomialRing(QuadraticField(2).maximal_order())
            sage: P.<x, y> = ProjectiveSpace(L.fraction_field() , 1)
            sage: H = Hom(P,P)
            sage: f = H([x^2 + (t^2 + 1)*y^2 , y^2])
            sage: f.dynatomic_polynomial(2)
            x^2 + x*y + (t^2 + 2)*y^2

        TESTS:

        We check that the dynatomic polynomial has the right parent (see :trac:`18409`)::

            sage: P.<x,y> = ProjectiveSpace(QQbar,1)
            sage: H = End(P)
            sage: R = P.coordinate_ring()
            sage: f = H([x^2-1/3*y^2, y^2])
            sage: f.dynatomic_polynomial(2).parent()
            Multivariate Polynomial Ring in x, y over Algebraic Field

        ::

            sage: T.<v> = QuadraticField(33)
            sage: S.<t> = PolynomialRing(T)
            sage: P.<x,y> = ProjectiveSpace(FractionField(S),1)
            sage: H = End(P)
            sage: f = H([t*x^2-1/t*y^2, y^2])
            sage: f.dynatomic_polynomial([1, 2]).parent()
            Multivariate Polynomial Ring in x, y over Fraction Field of Univariate Polynomial
            Ring in t over Number Field in v with defining polynomial x^2 - 33

        ::

            sage: P.<x, y> = ProjectiveSpace(QQ, 1)
            sage: H = End(P)
            sage: f = H([x^3-y^3*2, y^3])
            sage: f.dynatomic_polynomial(1).parent()
            Multivariate Polynomial Ring in x, y over Rational Field

        ::

            sage: P.<x, y> = ProjectiveSpace(QQ, 1)
            sage: H = End(P)
            sage: f = H([x^2 + y^2, y^2])
            sage: f.dynatomic_polynomial(0)
            0
            sage: f.dynatomic_polynomial([0,0])
            0

        Some rings still return Symoblic Ring elements::

            sage: S.<t> = FunctionField(CC)
            sage: P.<x,y> = ProjectiveSpace(S,1)
            sage: H = End(P)
            sage: R = P.coordinate_ring()
            sage: f = H([t*x^2-1*y^2, t*y^2])
            sage: f.dynatomic_polynomial([1, 2]).parent()
            Symbolic Ring

        ::

            sage: R.<x,y> = PolynomialRing(QQ)
            sage: S = R.quo(R.ideal(y^2-x+1))
            sage: P.<u,v> = ProjectiveSpace(FractionField(S),1)
            sage: H = End(P)
            sage: f = H([u^2 + S(x^2)*v^2, v^2])
            sage: dyn = f.dynatomic_polynomial([1,1]); dyn
            v^3*xbar^2 + u^2*v + u*v^2
            sage: dyn.parent()
            Symbolic Ring
       """
        if self.domain().ngens() > 2:
            raise TypeError("does not make sense in dimension >1")
        if not isinstance(period, (list, tuple)):
            period = [0, period]
        x = self.domain().gen(0)
        y = self.domain().gen(1)
        f0, f1 = F0, F1 = self._polys
        PHI = self.base_ring().one()
        m = period[0]
        n = period[1]
        if m == 0:
            if n == 0:
                return self[0].parent().zero()
            elif n == 1:
                return y*F0 - x*F1
        for d in range(1, n):
            if n % d == 0:
                PHI = PHI * ((y*F0 - x*F1)**moebius(n//d))
            F0, F1 = f0(F0, F1), f1(F0, F1)
        PHI = PHI * (y*F0 - x*F1)
        if m != 0:
            fm = self.nth_iterate_map(m)
            fm1 = self.nth_iterate_map(m - 1)
            PHI = PHI(fm._polys)/ PHI(fm1._polys)
        try:
            QR = PHI.numerator().quo_rem(PHI.denominator())
            if not QR[1]:
                return(QR[0])
        except (TypeError, NotImplementedError): # something Singular can't handle
            pass
        #even when the ring can be passed to singular in quo_rem,
        #it can't always do the division, so we call Maxima
        from sage.rings.padics.generic_nodes import is_pAdicField, is_pAdicRing
        if period != [0,1]: #period==[0,1] we don't need to do any division
            BR = self.domain().base_ring().base_ring()
            if not (is_pAdicRing(BR) or is_pAdicField(BR)):
                try:
                    PHI = PHI.numerator()._maxima_().divide(PHI.denominator())[0].sage()
                    # do it again to divide out by denominators of coefficients
                    PHI = PHI.numerator()._maxima_().divide(PHI.denominator())[0].sage()
                    if not is_FractionFieldElement(PHI):
                        from sage.symbolic.expression_conversions import polynomial
                        PHI = polynomial(PHI, ring=self.coordinate_ring())
                except (TypeError, NotImplementedError): #something Maxima, or the conversion, can't handle
                    pass
        return PHI

    def nth_iterate_map(self, n):
        r"""
        Returns the ``n``-th iterate of this map as a new map.

        ALGORITHM:

        Uses a form of successive squaring to reducing computations.


        .. TODO:: This could be improved.

        INPUT:

        - ``n`` -- a positive integer.

        OUTPUT:

        - A map between projective spaces.

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = Hom(P,P)
            sage: f = H([x^2+y^2, y^2])
            sage: f.nth_iterate_map(2)
            Scheme endomorphism of Projective Space of dimension 1 over Rational
            Field
              Defn: Defined on coordinates by sending (x : y) to
                    (x^4 + 2*x^2*y^2 + 2*y^4 : y^4)

        ::

            sage: P.<x,y> = ProjectiveSpace(CC,1)
            sage: H = Hom(P,P)
            sage: f = H([x^2-y^2, x*y])
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
            sage: f = H([x^2-y^2, x*y, z^2+x^2])
            sage: f.nth_iterate_map(2)
            Scheme endomorphism of Projective Space of dimension 2 over Integer Ring
              Defn: Defined on coordinates by sending (x : y : z) to
                    (x^4 - 3*x^2*y^2 + y^4 : x^3*y - x*y^3 : 2*x^4 - 2*x^2*y^2 + y^4
            + 2*x^2*z^2 + z^4)

        ::

            sage: P.<x,y,z> = ProjectiveSpace(QQ,2)
            sage: X = P.subscheme(x*z-y^2)
            sage: H = Hom(X,X)
            sage: f = H([x^2, x*z, z^2])
            sage: f.nth_iterate_map(2)
            Scheme endomorphism of Closed subscheme of Projective Space of dimension
            2 over Rational Field defined by:
              -y^2 + x*z
              Defn: Defined on coordinates by sending (x : y : z) to
                    (x^4 : x^2*z^2 : z^4)

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: P2.<u,v,w> = ProjectiveSpace(QQ, 2)
            sage: H = Hom(P, P2)
            sage: f = H([x^2, y^2, x^2-y^2])
            sage: f.nth_iterate_map(1)
            Scheme morphism:
              From: Projective Space of dimension 1 over Rational Field
              To:   Projective Space of dimension 2 over Rational Field
              Defn: Defined on coordinates by sending (x : y) to
                    (x^2 : y^2 : x^2 - y^2)
        """
        D = int(n)
        if D < 0:
            raise TypeError("iterate number must be a positive integer")
        if D == 1:
            return self
        if not self.is_endomorphism():
            raise TypeError("map is not an endomorphism")

        N = self.codomain().ambient_space().dimension_relative() + 1
        F = list(self._polys)
        Coord_ring = self.codomain().coordinate_ring()
        if isinstance(Coord_ring, QuotientRing_generic):
            PHI = [Coord_ring.gen(i).lift() for i in range(N)]
        else:
            PHI = [Coord_ring.gen(i) for i in range(N)]

        while D:
            if D&1:
                PHI = [PHI[j](*F) for j in range(N)]
            if D > 1: #avoid extra iterate
                F = [F[j](*F) for j in range(N)] #'square'
            D >>= 1
        return End(self.domain())(PHI)

    def nth_iterate(self, P, n, **kwds):
        r"""
        Returns the ``n``-th iterate of the point ``P`` by this map.

        If ``normalize`` is ``True``, then the coordinates are
        automatically normalized.

        .. TODO:: Is there a more efficient way to do this?

        INPUT:

        - ``P`` -- a point in this map's domain.

        - ``n`` -- a positive integer.

        kwds:

        - ``normalize`` - Boolean (optional Default: ``False``).

        OUTPUT:

        - A point in this map's codomain.

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(ZZ,1)
            sage: H = Hom(P,P)
            sage: f = H([x^2+y^2, 2*y^2])
            sage: Q = P(1,1)
            sage: f.nth_iterate(Q,4)
            (32768 : 32768)

        ::

            sage: P.<x,y> = ProjectiveSpace(ZZ,1)
            sage: H = Hom(P,P)
            sage: f = H([x^2+y^2, 2*y^2])
            sage: Q = P(1,1)
            sage: f.nth_iterate(Q, 4, normalize=True)
            (1 : 1)

        ::

            sage: P.<x,y,z> = ProjectiveSpace(QQ,2)
            sage: H = Hom(P,P)
            sage: f = H([x^2, 2*y^2, z^2-x^2])
            sage: Q = P(2,7,1)
            sage: f.nth_iterate(Q,2)
            (-16/7 : -2744 : 1)

        ::

            sage: R.<t> = PolynomialRing(QQ)
            sage: P.<x,y,z> = ProjectiveSpace(R,2)
            sage: H = Hom(P,P)
            sage: f = H([x^2+t*y^2, (2-t)*y^2, z^2])
            sage: Q = P(2+t,7,t)
            sage: f.nth_iterate(Q,2)
            (t^4 + 2507*t^3 - 6787*t^2 + 10028*t + 16 : -2401*t^3 + 14406*t^2 -
            28812*t + 19208 : t^4)

        ::

            sage: P.<x,y,z> = ProjectiveSpace(ZZ,2)
            sage: X = P.subscheme(x^2-y^2)
            sage: H = Hom(X,X)
            sage: f = H([x^2, y^2, z^2])
            sage: f.nth_iterate(X(2,2,3), 3)
            (256 : 256 : 6561)

        ::

            sage: K.<c> = FunctionField(QQ)
            sage: P.<x,y> = ProjectiveSpace(K,1)
            sage: H = Hom(P,P)
            sage: f = H([x^3 - 2*x*y^2 - c*y^3, x*y^2])
            sage: f.nth_iterate(P(c,1), 2)
            ((c^6 - 9*c^4 + 25*c^2 - c - 21)/(c^2 - 3) : 1)
        """
        return(P.nth_iterate(self, n, **kwds))

    def degree(self):
        r"""
        Return the degree of this map.

        The degree is defined as the degree of the homogeneous
        polynomials that are the coordinates of this map.

        OUTPUT:

        - A positive integer

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = Hom(P,P)
            sage: f = H([x^2+y^2, y^2])
            sage: f.degree()
            2

        ::

            sage: P.<x,y,z> = ProjectiveSpace(CC,2)
            sage: H = Hom(P,P)
            sage: f = H([x^3+y^3, y^2*z, z*x*y])
            sage: f.degree()
            3

        ::

            sage: R.<t> = PolynomialRing(QQ)
            sage: P.<x,y,z> = ProjectiveSpace(R,2)
            sage: H = Hom(P,P)
            sage: f = H([x^2+t*y^2, (2-t)*y^2, z^2])
            sage: f.degree()
            2

        ::

            sage: P.<x,y,z> = ProjectiveSpace(ZZ,2)
            sage: X = P.subscheme(x^2-y^2)
            sage: H = Hom(X,X)
            sage: f = H([x^2, y^2, z^2])
            sage: f.degree()
            2
        """
        return(self._polys[0].degree())

    def dehomogenize(self, n):
        r"""
        Returns the standard dehomogenization at the ``n[0]`` coordinate for the domain
        and the ``n[1]`` coordinate for the codomain.

        Note that the new function is defined over the fraction field
        of the base ring of this map.

        INPUT:

        - ``n`` -- a tuple of nonnegative integers.  If ``n`` is an integer, then the two values of
            the tuple are assumed to be the same.

        OUTPUT:

        - :class:`SchemeMorphism_polynomial_affine_space`.

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(ZZ,1)
            sage: H = Hom(P,P)
            sage: f = H([x^2+y^2, y^2])
            sage: f.dehomogenize(0)
            Scheme endomorphism of Affine Space of dimension 1 over Integer Ring
              Defn: Defined on coordinates by sending (x) to
                    (x^2/(x^2 + 1))

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = Hom(P,P)
            sage: f = H([x^2-y^2, y^2])
            sage: f.dehomogenize((0,1))
            Scheme endomorphism of Affine Space of dimension 1 over Rational Field
              Defn: Defined on coordinates by sending (x) to
                    ((-x^2 + 1)/x^2)

        ::

            sage: P.<x,y,z> = ProjectiveSpace(QQ,2)
            sage: H = Hom(P,P)
            sage: f = H([x^2+y^2, y^2-z^2, 2*z^2])
            sage: f.dehomogenize(2)
            Scheme endomorphism of Affine Space of dimension 2 over Rational Field
              Defn: Defined on coordinates by sending (x0, x1) to
                    (1/2*x0^2 + 1/2*x1^2, 1/2*x1^2 - 1/2)

        ::

            sage: R.<t> = PolynomialRing(QQ)
            sage: P.<x,y,z> = ProjectiveSpace(FractionField(R),2)
            sage: H = Hom(P,P)
            sage: f = H([x^2+t*y^2, t*y^2-z^2, t*z^2])
            sage: f.dehomogenize(2)
            Scheme endomorphism of Affine Space of dimension 2 over Fraction Field
            of Univariate Polynomial Ring in t over Rational Field
              Defn: Defined on coordinates by sending (x0, x1) to
                    (1/t*x0^2 + x1^2, x1^2 - 1/t)

        ::

            sage: P.<x,y,z> = ProjectiveSpace(ZZ,2)
            sage: X = P.subscheme(x^2-y^2)
            sage: H = Hom(X,X)
            sage: f = H([x^2, y^2, x*z])
            sage: f.dehomogenize(2)
            Scheme endomorphism of Closed subscheme of Affine Space of dimension 2 over Integer Ring defined by:
              x0^2 - x1^2
              Defn: Defined on coordinates by sending (x0, x1) to
                    (x0, x1^2/x0)

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = End(P)
            sage: f = H([x^2 - 2*x*y, y^2])
            sage: f.dehomogenize(0).homogenize(0) == f
            True
        """
        #the dehomogenizations are stored for future use.
        try:
            return self.__dehomogenization[n]
        except AttributeError:
            self.__dehomogenization = {}
        except KeyError:
            pass
        #it is possible to dehomogenize the domain and codomain at different coordinates
        if isinstance(n,(tuple,list)):
            ind=tuple(n)
        else:
            ind=(n,n)
        PS_domain = self.domain()
        A_domain = PS_domain.ambient_space()
        if self._polys[ind[1]].substitute({A_domain.gen(ind[0]):1}) == 0:
            raise ValueError("can't dehomogenize at 0 coordinate")
        else:
            Aff_domain = PS_domain.affine_patch(ind[0])
            S = Aff_domain.ambient_space().coordinate_ring()
            N = A_domain.dimension_relative()
            R = A_domain.coordinate_ring()
            phi = R.hom([S.gen(j) for j in range(0, ind[0])] + [1] + [S.gen(j) for j in range(ind[0], N)], S)
            F = []
            G = phi(self._polys[ind[1]])
            for i in range(0, N + 1):
                if i != ind[1]:
                    F.append(phi(self._polys[i]) / G)
            H = Hom(Aff_domain, self.codomain().affine_patch(ind[1]))
            #since often you dehomogenize at the same coordinate in domain
            #and codomain it should be stored appropriately.
            if ind == (n,n):
                self.__dehomogenization[ind]=H(F)
                return self.__dehomogenization[ind]
            else:
                self.__dehomogenization[n]=H(F)
                return self.__dehomogenization[n]

    def orbit(self, P, N, **kwds):
        r"""
        Return the orbit of the point ``P`` by this map.

        If ``N`` is an integer it returns `[P,self(P),\ldots,self^N(P)]`.
        If ``N`` is a list or tuple `N=[m,k]` it returns `[self^m(P),\ldots,self^k(P)]`.
        Automatically normalize the points if ``normalize=True``. Perform the checks on point initialize if
        ``check=True``

        INPUT:

        - ``P`` -- a point in this map's domain.

        - ``n`` -- a non-negative integer or list or tuple of two non-negative integers.

        kwds:

        - ``check`` -- boolean (optional - default: ``True``).

        - ``normalize`` -- boolean (optional - default: ``False``).


        OUTPUT:

        - a list of points in this map's codomain.

        EXAMPLES::

            sage: P.<x,y,z> = ProjectiveSpace(ZZ,2)
            sage: H = Hom(P,P)
            sage: f = H([x^2+y^2, y^2-z^2, 2*z^2])
            sage: f.orbit(P(1,2,1), 3)
            [(1 : 2 : 1), (5 : 3 : 2), (34 : 5 : 8), (1181 : -39 : 128)]

        ::

            sage: P.<x,y,z> = ProjectiveSpace(ZZ,2)
            sage: H = Hom(P,P)
            sage: f = H([x^2+y^2, y^2-z^2, 2*z^2])
            sage: f.orbit(P(1,2,1), [2,4])
            [(34 : 5 : 8), (1181 : -39 : 128), (1396282 : -14863 : 32768)]

        ::

            sage: P.<x,y,z> = ProjectiveSpace(ZZ,2)
            sage: X = P.subscheme(x^2-y^2)
            sage: H = Hom(X,X)
            sage: f = H([x^2, y^2, x*z])
            sage: f.orbit(X(2,2,3), 3, normalize=True)
            [(2 : 2 : 3), (2 : 2 : 3), (2 : 2 : 3), (2 : 2 : 3)]

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = Hom(P,P)
            sage: f = H([x^2+y^2, y^2])
            sage: f.orbit(P.point([1,2],False), 4, check=False)
            [(1 : 2), (5 : 4), (41 : 16), (1937 : 256), (3817505 : 65536)]

        ::

            sage: K.<c> = FunctionField(QQ)
            sage: P.<x,y> = ProjectiveSpace(K,1)
            sage: H = Hom(P,P)
            sage: f = H([x^2+c*y^2, y^2])
            sage: f.orbit(P(0,1), 3)
            [(0 : 1), (c : 1), (c^2 + c : 1), (c^4 + 2*c^3 + c^2 + c : 1)]
        """
        return(P.orbit(self, N, **kwds))

    @cached_method
    def is_morphism(self):
        r"""
        returns ``True`` if this map is a morphism.

        The map is a morphism if and only if the ideal generated by
        the defining polynomials is the unit ideal
        (no common zeros of the defining polynomials).

        OUTPUT:

        - Boolean

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = Hom(P,P)
            sage: f = H([x^2+y^2, y^2])
            sage: f.is_morphism()
            True

        ::

            sage: P.<x,y,z> = ProjectiveSpace(RR,2)
            sage: H = Hom(P,P)
            sage: f = H([x*z-y*z, x^2-y^2, z^2])
            sage: f.is_morphism()
            False

        ::

            sage: R.<t> = PolynomialRing(GF(5))
            sage: P.<x,y,z> = ProjectiveSpace(R,2)
            sage: H = Hom(P,P)
            sage: f = H([x*z-t*y^2, x^2-y^2, t*z^2])
            sage: f.is_morphism()
            True

        Map that is not morphism on projective space, but is over a subscheme::

            sage: P.<x,y,z> = ProjectiveSpace(RR,2)
            sage: X = P.subscheme([x*y + y*z])
            sage: H = Hom(X,X)
            sage: f = H([x*z-y*z, x^2-y^2, z^2])
            sage: f.is_morphism()
            True
        """

        R = self.coordinate_ring()
        F = self._polys
        defpolys = list(self.domain().defining_polynomials())
        if R.base_ring().is_field():
            F.extend(defpolys)
            J = R.ideal(F)
        else:
            S = PolynomialRing(R.base_ring().fraction_field(), R.gens(), R.ngens())
            L = [S(f) for f in F] + [S(f) for f in defpolys]
            J = S.ideal(L)
        if J.dimension() > 0:
            return False
        else:
            return True

    def resultant(self, normalize=False):
        r"""
        Computes the resultant of the defining polynomials of this map.

        If ``normalize`` is ``True``, then first normalize the coordinate
        functions with :meth:`normalize_coordinates`.

        INPUT:

        - ``normalize`` -- Boolean (optional - default: ``False``).

        OUTPUT:

        - an element of the base ring of this map.

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = Hom(P,P)
            sage: f = H([x^2+y^2, 6*y^2])
            sage: f.resultant()
            36

        ::

            sage: R.<t> = PolynomialRing(GF(17))
            sage: P.<x,y> = ProjectiveSpace(R,1)
            sage: H = Hom(P,P)
            sage: f = H([t*x^2+t*y^2, 6*y^2])
            sage: f.resultant()
            2*t^2

        ::

            sage: R.<t> = PolynomialRing(GF(17))
            sage: P.<x,y,z> = ProjectiveSpace(R,2)
            sage: H = Hom(P,P)
            sage: f = H([t*x^2+t*y^2, 6*y^2, 2*t*z^2])
            sage: f.resultant()
            13*t^8

        ::

            sage: P.<x,y,z> = ProjectiveSpace(QQ,2)
            sage: H = Hom(P,P)
            sage: F = H([x^2+y^2,6*y^2,10*x*z+z^2+y^2])
            sage: F.resultant()
            1296

        ::

            sage: R.<t>=PolynomialRing(QQ)
            sage: s = (t^3+t+1).roots(QQbar)[0][0]
            sage: P.<x,y>=ProjectiveSpace(QQbar,1)
            sage: H = Hom(P,P)
            sage: f = H([s*x^3-13*y^3, y^3-15*y^3])
            sage: f.resultant()
            871.6925062959149?
            """

        if self.domain().dimension_relative() != self.codomain().dimension_relative():
            raise ValueError("domain and codomain should be of same dimension")
        if normalize:
            F = copy(self)
            F.normalize_coordinates()
        else:
            F = self

        if self.domain().dimension_relative() == 1:
            x = self.domain().gen(0)
            y = self.domain().gen(1)
            d = self.degree()
            f = F[0].substitute({y:1})
            g = F[1].substitute({y:1})
            #Try to use pari first, as it is faster for one dimensional case
            #however the coercion from a Pari object to a sage object breaks
            #in the case of QQbar, so we just pass it into the macaulay resultant
            try:
                res = (f.lc() ** (d - g.degree()) * g.lc() ** (d - f.degree()) * f._pari_().polresultant(g, x))
                return(self.domain().base_ring()(res))
            except (TypeError, PariError):
                pass
        #Otherwise, use Macaulay
        R = F[0].parent()
        res = R.macaulay_resultant(F._polys)
        return res #Coercion here is not necessary as it is already done in Macaulay Resultant

    @cached_method
    def primes_of_bad_reduction(self, check=True):
        r"""
        Determines the primes of bad reduction for an endomorphism
        defined over number fields.

        If ``check`` is ``True``, each prime is verified to be of bad reduction.

        ALGORITHM:

        `p` is a prime of bad reduction if and only if the defining
        polynomials of self have a common zero. Or stated another way,
        `p` is a prime of bad reduction if and only if the radical of
        the ideal defined by the defining polynomials of self is not
        `(x_0,x_1,\ldots,x_N)`.  This happens if and only if some
        power of each `x_i` is not in the ideal defined by the
        defining polynomials of self. This last condition is what is
        checked. The lcm of the coefficients of the monomials `x_i` in
        a Groebner basis is computed. This may return extra primes.

        INPUT:

        - ``check`` -- Boolean (optional - default: ``True``).

        OUTPUT:

        - a list of integer primes.

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = Hom(P,P)
            sage: f = H([1/3*x^2+1/2*y^2, y^2])
            sage: f.primes_of_bad_reduction()
            [2, 3]

        ::

            sage: P.<x,y,z,w> = ProjectiveSpace(QQ,3)
            sage: H = Hom(P,P)
            sage: f = H([12*x*z-7*y^2, 31*x^2-y^2, 26*z^2, 3*w^2-z*w])
            sage: f.primes_of_bad_reduction()
            [2, 3, 7, 13, 31]

        A number field example ::

            sage: R.<z> = QQ[]
            sage: K.<a> = NumberField(z^2 - 2)
            sage: P.<x,y> = ProjectiveSpace(K,1)
            sage: H = Hom(P,P)
            sage: f = H([1/3*x^2+1/a*y^2, y^2])
            sage: f.primes_of_bad_reduction()
            [Fractional ideal (a), Fractional ideal (3)]

        This is an example where check = False returns extra primes::

            sage: P.<x,y,z> = ProjectiveSpace(ZZ,2)
            sage: H = Hom(P,P)
            sage: f = H([3*x*y^2 + 7*y^3 - 4*y^2*z + 5*z^3, -5*x^3 + x^2*y + y^3 + 2*x^2*z,\
                -2*x^2*y + x*y^2 + y^3 - 4*y^2*z + x*z^2])
            sage: f.primes_of_bad_reduction(False)
            [2, 5, 37, 2239, 304432717]
            sage: f.primes_of_bad_reduction()
            [5, 37, 2239, 304432717]
        """
        from sage.schemes.projective.projective_space import is_ProjectiveSpace
        if is_ProjectiveSpace(self.domain()) is False or is_ProjectiveSpace(self.codomain()) is False:
            raise NotImplementedError
        K = FractionField(self.codomain().base_ring())
        #The primes of bad reduction are the support of the resultant for number fields

        if K in NumberFields():
            if K != QQ:
                F = copy(self)
                F.normalize_coordinates()
                return (K(F.resultant()).support())
            else:
                #For the rationals, we can use groebner basis, as it is quicker in practice
                R = self.coordinate_ring()
                F = self._polys

                if R.base_ring().is_field():
                    J = R.ideal(F)
                else:
                    S = PolynomialRing(R.base_ring().fraction_field(), R.gens(), R.ngens())
                    J = S.ideal([S.coerce(F[i]) for i in range(R.ngens())])
                if J.dimension() > 0:
                    raise TypeError("not a morphism")
                #normalize to coefficients in the ring not the fraction field.
                F = [F[i] * lcm([F[j].denominator() for j in range(len(F))]) for i in range(len(F))]

                #move the ideal to the ring of integers
                if R.base_ring().is_field():
                    S = PolynomialRing(R.base_ring().ring_of_integers(), R.gens(), R.ngens())
                    F = [F[i].change_ring(R.base_ring().ring_of_integers()) for i in range(len(F))]
                    J = S.ideal(F)
                else:
                    J = R.ideal(F)
                GB = J.groebner_basis()
                badprimes = []

                #get the primes dividing the coefficients of the monomials x_i^k_i
                for i in range(len(GB)):
                    LT = GB[i].lt().degrees()
                    power = 0
                    for j in range(R.ngens()):
                        if LT[j] != 0:
                            power += 1
                    if power == 1:
                        badprimes = badprimes + GB[i].lt().coefficients()[0].support()
                badprimes = sorted(set(badprimes))

                #check to return only the truly bad primes
                if check:
                    index = 0
                    while index < len(badprimes):  #figure out which primes are really bad primes...
                        S = PolynomialRing(GF(badprimes[index]), R.gens(), R.ngens())
                        J = S.ideal([S.coerce(F[j]) for j in range(R.ngens())])
                        if J.dimension() == 0:
                            badprimes.pop(index)
                        else:
                            index += 1
                return(badprimes)
        else:
            raise TypeError("base ring must be number field or number field ring")

    def conjugate(self, M):
        r"""
        Conjugates this map by ``M``, i.e. `M^{-1} \circ f \circ M`.

        If possible the new map will be defined over the same space as
        this map. Otherwise, will try to coerce to the base ring of
        ``M``.

        INPUT:

        - ``M`` -- a square invertible matrix.

        OUTPUT:

        - a map from the domain to the codomain of this map.

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(ZZ,1)
            sage: H = Hom(P,P)
            sage: f = H([x^2+y^2, y^2])
            sage: f.conjugate(matrix([[1,2], [0,1]]))
            Scheme endomorphism of Projective Space of dimension 1 over Integer Ring
              Defn: Defined on coordinates by sending (x : y) to
                    (x^2 + 4*x*y + 3*y^2 : y^2)

        ::

            sage: R.<x> = PolynomialRing(QQ)
            sage: K.<i> = NumberField(x^2+1)
            sage: P.<x,y> = ProjectiveSpace(ZZ,1)
            sage: H = Hom(P,P)
            sage: f = H([x^3+y^3, y^3])
            sage: f.conjugate(matrix([[i,0], [0,-i]]))
            Scheme endomorphism of Projective Space of dimension 1 over Integer Ring
              Defn: Defined on coordinates by sending (x : y) to
                    (-x^3 + y^3 : -y^3)

        ::

            sage: P.<x,y,z> = ProjectiveSpace(ZZ,2)
            sage: H = Hom(P,P)
            sage: f = H([x^2+y^2 ,y^2, y*z])
            sage: f.conjugate(matrix([[1,2,3], [0,1,2], [0,0,1]]))
            Scheme endomorphism of Projective Space of dimension 2 over Integer Ring
              Defn: Defined on coordinates by sending (x : y : z) to
                    (x^2 + 4*x*y + 3*y^2 + 6*x*z + 9*y*z + 7*z^2 : y^2 + 2*y*z : y*z + 2*z^2)

        ::

            sage: P.<x,y> = ProjectiveSpace(ZZ,1)
            sage: H = Hom(P,P)
            sage: f = H([x^2+y^2, y^2])
            sage: f.conjugate(matrix([[2,0], [0,1/2]]))
            Scheme endomorphism of Projective Space of dimension 1 over Rational Field
              Defn: Defined on coordinates by sending (x : y) to
                    (2*x^2 + 1/8*y^2 : 1/2*y^2)

        ::

            sage: R.<x> = PolynomialRing(QQ)
            sage: K.<i> = NumberField(x^2+1)
            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = Hom(P,P)
            sage: f = H([1/3*x^2+1/2*y^2, y^2])
            sage: f.conjugate(matrix([[i,0], [0,-i]]))
            Scheme endomorphism of Projective Space of dimension 1 over Number Field in i with defining polynomial x^2 + 1
              Defn: Defined on coordinates by sending (x : y) to
                    ((1/3*i)*x^2 + (1/2*i)*y^2 : (-i)*y^2)
        """

        if M.is_square() == 1 and M.determinant() != 0 and  M.ncols() == self.domain().ambient_space().dimension_relative() + 1:
            X = M * vector(self[0].parent().gens())
            F = vector(self._polys)
            F = F(list(X))
            N = M.inverse()
            F = N * F
            R = self.codomain().ambient_space().coordinate_ring()
            try:
                F = [R(f) for f in F]
                PS = self.codomain()
            except TypeError: #no longer defined over same ring
                R = R.change_ring(M.base_ring())
                F = [R(f) for f in F]
                PS = self.codomain().change_ring(M.base_ring())
            H = Hom(PS, PS)
            return(H(F))
        else:
            raise TypeError("matrix must be invertible and size dimension +1 ")

    def green_function(self, P, v, **kwds):
        r"""
        Evaluates the local Green's function at the place ``v`` for ``P`` with ``N`` terms of the
        series or to within a given error bound.

        Must be over a number field or order of a number field. Note that this is
        the absolute local Green's function so is scaled by the degree of the base field.

        Use ``v=0`` for the archimedean place over `\QQ` or field embedding. Non-archimedean
        places are prime ideals for number fields or primes over `\QQ`.

        ALGORITHM:

        See Exercise 5.29 and Figure 5.6 of [Silverman-ADS]_.

        REFERENCES:

        .. [Silverman-ADS] Joseph H. Silverman. The Arithmetic of Dynamics Systems. Springer, GTM 241, 2007.

        INPUT:

        - ``P`` - a projective point.

        - ``v`` - non-negative integer. a place, use v=0 for the archimedean place.

        kwds:

        - ``N`` - positive integer. number of terms of the series to use, (optional - default: 10).

        - ``prec`` - positive integer, float point or p-adic precision, default: 100.

        - ``error_bound`` - a positive real number (optional).

        OUTPUT:

        - a real number.

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = Hom(P,P)
            sage: f = H([x^2+y^2, x*y])
            sage: f.green_function(P.point([5,2], False), 0, N=30)
            1.7315451844777407992085512000
            sage: f.green_function(P.point([2,1], False), 0, N=30)
            0.86577259223181088325226209926
            sage: f.green_function(P.point([1,1], False), 0, N=30)
            0.43288629610862338612700146098
        """
        return(P.green_function(self, v, **kwds))

    def canonical_height(self, P, **kwds):
        r"""
        Evaluates the (absolute) canonical height of ``P`` with respect to this map.

        Must be over number field or order of a number field. Specify either
        the number of terms of the series to evaluate or the error bound required.

        ALGORITHM:

        The sum of the Green's function at the archimedean places and the places of bad reduction.

        INPUT:

        - ``P`` -- a projective point.

        kwds:

        - ``badprimes`` - a list of primes of bad reduction (optional).

        - ``N`` - positive integer. number of terms of the series to use in the local green functions
          (optional - default: 10).

        - ``prec`` - positive integer, float point or p-adic precision, default: 100.

        - ``error_bound`` - a positive real number (optional).

        OUTPUT:

        - a real number

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(ZZ,1)
            sage: H = Hom(P,P)
            sage: f = H([x^2+y^2, 2*x*y]);
            sage: f.canonical_height(P.point([5,4]), error_bound=0.001)
            2.1970553519503404898926835324
            sage: f.canonical_height(P.point([2,1]), error_bound=0.001)
            1.0984430632822307984974382955

        Notice that preperiodic points may not be exactly 0::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = Hom(P,P)
            sage: f = H([x^2-29/16*y^2, y^2]);
            sage: f.canonical_height(P.point([1,4]), error_bound=0.000001)
            2.9868196689972114460185071428e-7

        ::

            sage: P.<x,y,z> = ProjectiveSpace(QQ,2)
            sage: X = P.subscheme(x^2-y^2);
            sage: H = Hom(X,X)
            sage: f = H([x^2,y^2, 4*z^2]);
            sage: Q = X([4,4,1])
            sage: f.canonical_height(Q, badprimes=[2])
            0.0013538030870311431824555314882
        """
        return(P.canonical_height(self, **kwds))

    def global_height(self, prec=None):
        r"""
        Returns the maximum of the absolute logarithmic heights of the coefficients
        in any of the coordinate functions of this map.

        INPUT:

        - ``prec`` -- desired floating point precision (default:
          default RealField precision).

        OUTPUT:

        - a real number.

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = Hom(P,P)
            sage: f = H([1/1331*x^2+1/4000*y^2, 210*x*y]);
            sage: f.global_height()
            8.29404964010203

        This function does not automatically normalize::

            sage: P.<x,y,z> = ProjectiveSpace(ZZ,2)
            sage: H = Hom(P,P)
            sage: f = H([4*x^2+100*y^2, 210*x*y, 10000*z^2]);
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
            sage: f = H([2*x^2 + 3*O(w)*y^2, O(w)*y^2])
            sage: f.global_height()
            1.44518587894808

        ::

            sage: P.<x,y> = ProjectiveSpace(QQbar,1)
            sage: P2.<u,v,w> = ProjectiveSpace(QQbar,2)
            sage: H = Hom(P,P2)
            sage: f = H([x^2 + QQbar(I)*x*y + 3*y^2, y^2, QQbar(sqrt(5))*x*y])
            sage: f.global_height()
            1.09861228866811
        """
        K = self.domain().base_ring()
        if K in _NumberFields or is_NumberFieldOrder(K):
            f = self
        elif K is QQbar:
            f = self._number_field_from_algebraics()
        else:
            raise TypeError("Must be over a Numberfield or a Numberfield Order or QQbar")
        H = 0
        for i in range(self.domain().ambient_space().dimension_relative() + 1):
            C = f[i].coefficients()
            h = max([c.global_height(prec) for c in C])
            H = max(H, h)
        return(H)

    def local_height(self, v, prec=None):
        r"""
        Returns the maximum of the local height of the coefficients in any
        of the coordinate functions of this map.

        INPUT:

        - ``v`` -- a prime or prime ideal of the base ring.

        - ``prec`` -- desired floating point precision (default:
          default RealField precision).

        OUTPUT:

        - a real number.

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = Hom(P,P)
            sage: f = H([1/1331*x^2+1/4000*y^2, 210*x*y]);
            sage: f.local_height(1331)
            7.19368581839511

        This function does not automatically normalize::

            sage: P.<x,y,z> = ProjectiveSpace(QQ,2)
            sage: H = Hom(P,P)
            sage: f = H([4*x^2+3/100*y^2, 8/210*x*y, 1/10000*z^2]);
            sage: f.local_height(2)
            2.77258872223978
            sage: f.normalize_coordinates()
            sage: f.local_height(2)
            0.000000000000000

        ::

            sage: R.<z> = PolynomialRing(QQ)
            sage: K.<w> = NumberField(z^2-2)
            sage: P.<x,y> = ProjectiveSpace(K,1)
            sage: H = Hom(P,P)
            sage: f = H([2*x^2 + w/3*y^2, 1/w*y^2])
            sage: f.local_height(K.ideal(3))
            1.09861228866811
        """
        K = FractionField(self.domain().base_ring())
        if K not in _NumberFields:
            raise TypeError("must be over a number field or a number field order")
        return max([K(c).local_height(v, prec) for f in self for c in f.coefficients()])

    def local_height_arch(self, i, prec=None):
        r"""
        Returns the maximum of the local height at the ``i``-th infinite place of the coefficients in any
        of the coordinate functions of this map.

        INPUT:

        - ``i`` -- an integer.

        - ``prec`` -- desired floating point precision (default:
          default RealField precision).

        OUTPUT:

        - a real number.

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = Hom(P,P)
            sage: f = H([1/1331*x^2+1/4000*y^2, 210*x*y]);
            sage: f.local_height_arch(0)
            5.34710753071747

        ::

            sage: R.<z> = PolynomialRing(QQ)
            sage: K.<w> = NumberField(z^2-2)
            sage: P.<x,y> = ProjectiveSpace(K,1)
            sage: H = Hom(P,P)
            sage: f = H([2*x^2 + w/3*y^2, 1/w*y^2])
            sage: f.local_height_arch(1)
            0.6931471805599453094172321214582
        """
        K = FractionField(self.domain().base_ring())
        if K not in _NumberFields:
            raise TypeError("must be over a number field or a number field order")
        if K == QQ:
            return max([K(c).local_height_arch(prec=prec) for f in self for c in f.coefficients()])
        else:
            return max([K(c).local_height_arch(i, prec=prec) for f in self for c in f.coefficients()])


    def height_difference_bound(self, prec=None):
        r"""
        Return an upper bound on the different between the canonical
        height of a point with respect to this map and the absolute
        height of the point.

        This map must be a morphism.

        ALGORITHM:

            Uses a Nullstellensatz argument to compute the constant.
            For details: see [Hutz]_.

        INPUT:

        - ``prec`` - positive integer, float point, default: RealField default.

        OUTPUT:

        - a real number.

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = End(P)
            sage: f = H([x^2+y^2, x*y]);
            sage: f.height_difference_bound()
            1.38629436111989

        This function does not automatically normalize. ::

            sage: P.<x,y,z> = ProjectiveSpace(ZZ,2)
            sage: H = End(P)
            sage: f = H([4*x^2+100*y^2, 210*x*y, 10000*z^2]);
            sage: f.height_difference_bound()
            11.0020998412042
            sage: f.normalize_coordinates()
            sage: f.height_difference_bound()
            10.3089526606443

       A number field example::

            sage: R.<x> = QQ[]
            sage: K.<c> = NumberField(x^3 - 2)
            sage: P.<x,y,z> = ProjectiveSpace(K,2)
            sage: H = End(P)
            sage: f = H([1/(c+1)*x^2+c*y^2, 210*x*y, 10000*z^2])
            sage: f.height_difference_bound()
            11.0020998412042

        ::

            sage: P.<x,y,z> = ProjectiveSpace(QQbar,2)
            sage: H = End(P)
            sage: f = H([x^2, QQbar(sqrt(-1))*y^2, QQbar(sqrt(3))*z^2])
            sage: f.height_difference_bound()
            3.43967790223022
        """
        FF = FractionField(self.domain().base_ring()) #lift will only work over fields, so coercing into FF
        if not FF in NumberFields():
            if FF == QQbar:
                #since this is absolute height, we can choose any number field over which the
                #function is defined.
                f = self._number_field_from_algebraics()
            else:
                raise NotImplementedError("fraction field of the base ring must be a number field or QQbar")
        else:
            f = self.change_ring(FF)
        if not self.is_endomorphism():
            raise NotImplementedError("must be an endomorphism of projective space")
        if prec is None:
            R = RealField()
        else:
            R = RealField(prec)
        N = f.domain().dimension_relative()
        d = f.degree()
        D = (N + 1) * (d - 1) + 1
        #compute upper bound
        U = f.global_height(prec) + R(binomial(N + d, d)).log()
        #compute lower bound - from explicit polynomials of Nullstellensatz
        CR = f.domain().coordinate_ring()
        I = CR.ideal(f.defining_polynomials())
        MCP = []
        for k in range(N + 1):
            CoeffPolys = (CR.gen(k) ** D).lift(I)
            Res = 1
            for j in range(len(CoeffPolys)):
                if CoeffPolys[j] != 0:
                    for i in range(len(CoeffPolys[j].coefficients())):
                        Res = lcm(Res, abs(CoeffPolys[j].coefficients()[i].denominator()))
            h = max([c.global_height() for g in CoeffPolys for c in (Res*g).coefficients()])
            MCP.append([Res, h]) #since we need to clear denominators
        maxh = 0
        gcdRes = 0
        for k in range(len(MCP)):
            gcdRes = gcd(gcdRes, MCP[k][0])
            maxh = max(maxh, MCP[k][1])
        L = abs( R(gcdRes).log() - R((N + 1) * binomial(N + D - d, D - d)).log() - maxh)
        C = max(U, L) #height difference dh(P) - L <= h(f(P)) <= dh(P) +U
        return(C / (d - 1))

    def multiplier(self, P, n, check=True):
        r"""
        Returns the multiplier of the point ``P`` of period ``n`` with respect to this map.

        This map must have same domain and codomain.

        INPUT:

        - ``P`` - a point on domain of this map.

        - ``n`` - a positive integer, the period of ``P``.

        - ``check`` -- verify that ``P`` has period ``n``, Default:True.

        OUTPUT:

        - a square matrix of size ``self.codomain().dimension_relative()``
          in the ``base_ring`` of this map.

        EXAMPLES::

            sage: P.<x,y,z> = ProjectiveSpace(QQ,2)
            sage: H = End(P)
            sage: f = H([x^2,y^2, 4*z^2]);
            sage: Q = P.point([4,4,1], False);
            sage: f.multiplier(Q,1)
            [2 0]
            [0 2]

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = End(P)
            sage: f = H([7*x^2 - 28*y^2, 24*x*y])
            sage: f.multiplier(P(2,5), 4)
            [231361/20736]

        ::

            sage: P.<x,y> = ProjectiveSpace(CC,1)
            sage: H = End(P)
            sage: f = H([x^3 - 25*x*y^2 + 12*y^3, 12*y^3])
            sage: f.multiplier(P(1,1), 5)
            [0.389017489711935]

        ::

            sage: P.<x,y> = ProjectiveSpace(RR,1)
            sage: H = End(P)
            sage: f = H([x^2-2*y^2, y^2])
            sage: f.multiplier(P(2,1), 1)
            [4.00000000000000]

        ::

            sage: P.<x,y> = ProjectiveSpace(Qp(13),1)
            sage: H = End(P)
            sage: f = H([x^2-29/16*y^2, y^2])
            sage: f.multiplier(P(5,4), 3)
            [6 + 8*13 + 13^2 + 8*13^3 + 13^4 + 8*13^5 + 13^6 + 8*13^7 + 13^8 +
            8*13^9 + 13^10 + 8*13^11 + 13^12 + 8*13^13 + 13^14 + 8*13^15 + 13^16 +
            8*13^17 + 13^18 + 8*13^19 + O(13^20)]

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = End(P)
            sage: f = H([x^2-y^2, y^2])
            sage: f.multiplier(P(0,1), 1)
            Traceback (most recent call last):
            ...
            ValueError: (0 : 1) is not periodic of period 1
        """
        if not self.is_endomorphism():
            raise TypeError("must be an endomorphism")
        if check:
            if self.nth_iterate(P, n) != P:
                raise ValueError("%s is not periodic of period %s"%(P, n))
            if n < 1:
                raise ValueError("period must be a positive integer")
        N = self.domain().ambient_space().dimension_relative()
        l = identity_matrix(FractionField(self.codomain().base_ring()), N, N)
        Q = P
        Q.normalize_coordinates()
        index = N
        indexlist = [] #keep track of which dehomogenizations are needed
        while Q[index] == 0:
            index -= 1
        indexlist.append(index)
        for i in range(0, n):
            F = []
            R = self(Q)
            R.normalize_coordinates()
            index = N
            while R[index] == 0:
                index -= 1
            indexlist.append(index)
            #dehomogenize and compute multiplier
            F = self.dehomogenize((indexlist[i],indexlist[i+1]))
            #get the correct order for chain rule matrix multiplication
            l = F.jacobian()(tuple(Q.dehomogenize(indexlist[i])))*l
            Q = R
        return l

    def _multipliermod(self, P, n, p, k):
        r"""
        Returns the multiplier of the point ``P`` of period ``n`` with respect to
        this map modulo `p^k`.

        This map must be an endomorphism of projective space defined over `\QQ` or '\ZZ'.
        This function should not be used at the top level as it does not perform input checks.
        It is used primarily for the rational preperiodic and periodic point algorithms.

        INPUT:

        - ``P`` - a point on domain of this map.

        - ``n`` - a positive integer, the period of ``P``.

        - ``p`` - a positive integer.

        - ``k`` - a positive integer.

        OUTPUT:

        - a square matrix of size ``self.codomain().dimension_relative()`` in `Zmod(p^k)`.

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = End(P)
            sage: f = H([x^2-29/16*y^2, y^2])
            sage: f._multipliermod(P(5,4), 3, 11, 1)
            [3]

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = End(P)
            sage: f = H([x^2-29/16*y^2, y^2])
            sage: f._multipliermod(P(5,4), 3, 11, 2)
            [80]
        """
        N = self.domain().dimension_relative()
        BR = FractionField(self.codomain().base_ring())
        l = identity_matrix(BR, N, N)
        Q = copy(P)
        g = gcd(Q._coords) #we can't use normalize_coordinates since it can cause denominators
        Q.scale_by(1 / g)
        index = N
        indexlist = [] #keep track of which dehomogenizations are needed
        while Q[index] % p == 0:
            index -= 1
        indexlist.append(index)
        for i in range(0, n):
            F = []
            R = self(Q, False)
            g = gcd(R._coords)
            R.scale_by(1 / g)
            for index in range(N + 1):
                R._coords[index] = R._coords[index] % (p ** k)
            index = N
            while R[index] % p == 0:
                index -= 1
            indexlist.append(index)
            #dehomogenize and compute multiplier
            F = self.dehomogenize((indexlist[i],indexlist[i+1]))
            l = (F.jacobian()(tuple(Q.dehomogenize(indexlist[i])))*l) % (p ** k)
            Q = R
        return(l)

    def possible_periods(self, **kwds):
        r"""
        Returns the set of possible periods for rational periodic points of this map.

        Must be defined over `\ZZ` or `\QQ`.

        ALGORITHM:
            Calls ``self.possible_periods()`` modulo all primes of good reduction in range ``prime_bound``.
            Returns the intersection of those lists.

        INPUT:

        kwds:

        - ``prime_bound`` - a list or tuple of two positive integers. Or an integer for the upper bound. (optional)
            default: [1,20].

        - ``bad_primes`` - a list or tuple of integer primes, the primes of bad reduction.  (optional)

        - ``ncpus`` - number of cpus to use in parallel.  (optional)
            default: all available cpus.

        OUTPUT:

        - a list of positive integers.

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = End(P)
            sage: f = H([x^2-29/16*y^2, y^2])
            sage: f.possible_periods(ncpus=1)
            [1, 3]

        ::

            sage: PS.<x,y> = ProjectiveSpace(1,QQ)
            sage: H = End(PS)
            sage: f = H([5*x^3 - 53*x*y^2 + 24*y^3, 24*y^3])
            sage: f.possible_periods(prime_bound=[1,5])
            Traceback (most recent call last):
            ...
            ValueError: no primes of good reduction in that range
            sage: f.possible_periods(prime_bound=[1,10])
            [1, 4, 12]
            sage: f.possible_periods(prime_bound=[1,20])
            [1, 4]

        ::

            sage: P.<x,y,z> = ProjectiveSpace(ZZ,2)
            sage: H = End(P)
            sage: f = H([2*x^3 - 50*x*z^2 + 24*z^3, 5*y^3 - 53*y*z^2 + 24*z^3, 24*z^3])
            sage: f.possible_periods(prime_bound=10)
            [1, 2, 6, 20, 42, 60, 140, 420]
            sage: f.possible_periods(prime_bound=20) # long time
            [1, 20]
        """
        if not self.is_endomorphism():
            raise NotImplementedError("must be an endomorphism of projective space")
        if self.domain().base_ring() != ZZ and self.domain().base_ring() != QQ:
            raise NotImplementedError("must be ZZ or QQ")


        primebound = kwds.pop("prime_bound", [1, 20])
        badprimes = kwds.pop("bad_primes", None)
        num_cpus = kwds.pop("ncpus", ncpus())

        if (isinstance(primebound, (list, tuple)) == False):
            try:
                primebound = [1, ZZ(primebound)]
            except TypeError:
                raise TypeError("prime bound must be an integer")
        else:
            try:
                primebound[0] = ZZ(primebound[0])
                primebound[1] = ZZ(primebound[1])
            except TypeError:
                raise TypeError("prime bounds must be integers")

        if badprimes is None:
            badprimes = self.primes_of_bad_reduction()

        firstgood = 0

        def parallel_function(morphism):
            return morphism.possible_periods()

        # Calling possible_periods for each prime in parallel
        parallel_data = []
        for q in primes(primebound[0], primebound[1] + 1):
            if not (q in badprimes):
                F = self.change_ring(GF(q))
                parallel_data.append(((F,), {}))

        parallel_iter = p_iter_fork(num_cpus, 0)
        parallel_results = list(parallel_iter(parallel_function, parallel_data))

        for result in parallel_results:
            possible_periods = result[1]
            if firstgood == 0:
                periods = set(possible_periods)
                firstgood = 1
            else:
                periodsq = set(possible_periods)
                periods = periods.intersection(periodsq)

        if firstgood == 0:
            raise ValueError("no primes of good reduction in that range")
        else:
            return(sorted(periods))

    def _preperiodic_points_to_cyclegraph(self, preper):
        r"""
        Given the complete set of periodic or preperiodic points return the
        digraph representing the orbit.

        If ``preper`` is not the complete set, this function will not fill in the gaps.


        INPUT:

        - ``preper`` - a list or tuple of projective points. The complete set
          of rational periodic or preperiodic points.

        OUTPUT:

        - a digraph representing the orbit the rational preperiodic points
          ``preper`` in projective space.

        Examples::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = End(P)
            sage: f = H([x^2-2*y^2, y^2])
            sage: preper = [P(-2, 1), P(1, 0), P(0, 1), P(1, 1), P(2, 1), P(-1, 1)]
            sage: f._preperiodic_points_to_cyclegraph(preper)
            Looped digraph on 6 vertices
        """
        V = []
        E = []
        #We store the points we encounter is a list, D. Each new point is checked to
        #see if it is in that list (which uses ==) so that equal points with different
        #representations only appear once in the graph.
        D=[]
        for i in range(0, len(preper)):
            try:
                V.append(D[D.index(preper[i])])
            except ValueError:
                D.append(preper[i])
                V.append(preper[i])
            Q = self(preper[i])
            Q.normalize_coordinates()
            try:
                E.append([D[D.index(Q)]])
            except ValueError:
                D.append(Q)
                E.append([Q])
        from sage.graphs.digraph import DiGraph
        g = DiGraph(dict(zip(V, E)), loops=True)
        return(g)

    def is_PGL_minimal(self, prime_list=None):
        r"""
        Checks if this map is a minimal model in its conjugacy class.

        See [Bruin-Molnar]_ and [Molnar]_ for a description of the algorithm.

        INPUT:

        - ``prime_list`` -- list of primes to check minimality, if None, check all places.

        OUTPUT:

        - Boolean - True if this map is minimal, False otherwise.

        EXAMPLES::

            sage: PS.<X,Y> = ProjectiveSpace(QQ,1)
            sage: H = End(PS)
            sage: f = H([X^2+3*Y^2, X*Y])
            sage: f.is_PGL_minimal()
            True

        ::

            sage: PS.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = End(PS)
            sage: f = H([6*x^2+12*x*y+7*y^2, 12*x*y])
            sage: f.is_PGL_minimal()
            False

        ::

            sage: PS.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = End(PS)
            sage: f = H([6*x^2+12*x*y+7*y^2, y^2])
            sage: f.is_PGL_minimal()
            Traceback (most recent call last):
            ...
            TypeError: affine minimality is only considered for maps not of the form
            f or 1/f for a polynomial f
        """
        if self.base_ring() != QQ and self.base_ring() != ZZ:
            raise NotImplementedError("minimal models only implemented over ZZ or QQ")
        if not self.is_morphism():
            raise TypeError("the function is not a morphism")
        if self.degree() == 1:
            raise NotImplementedError("minimality is only for degree 2 or higher")

        from .endPN_minimal_model import affine_minimal
        return(affine_minimal(self, False , prime_list , True))

    def minimal_model(self, return_transformation=False, prime_list=None):
        r"""
        Determine if this morphisms is minimal.

        This map must be defined over the projective line over the rationals.
        In particular, determine if this map is affine minimal, which is enough
        to decide if it is minimal or not. See Proposition 2.10 in [Bruin-Molnar].

        REFERENCES:

        [Bruin-Molnar]_, [Molnar]_

        INPUT:

        - ``return_transformation`` -- a boolean value, default value True. This
                                    signals a return of the `PGL_2` transformation
                                    to conjugate this map to the calculated minimal
                                    model. default: False.

        - ``prime_list`` -- a list of primes, in case one only wants to determine minimality
                   at those specific primes.

        OUTPUT:

        - a scheme morphism on the projective line which is a minimal model of this map.

        - a `PGL(2,\QQ)` element which conjugates this map to a minimal model.

        EXAMPLES::

            sage: PS.<X,Y> = ProjectiveSpace(QQ,1)
            sage: H = End(PS)
            sage: f = H([X^2+3*Y^2, X*Y])
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
            sage: f = H([7365/2*X^4 + 6282*X^3*Y + 4023*X^2*Y^2 + 1146*X*Y^3 + 245/2*Y^4,\
             -12329/2*X^4 - 10506*X^3*Y - 6723*X^2*Y^2 - 1914*X*Y^3 - 409/2*Y^4])
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
            sage: f = H([6*x^2+12*x*y+7*y^2, 12*x*y])
            sage: f.minimal_model()
            Scheme endomorphism of Projective Space of dimension 1 over Rational
            Field
              Defn: Defined on coordinates by sending (x : y) to
                    (x^2 + 12*x*y + 42*y^2 : 2*x*y)

        ::

            sage: PS.<x,y> = ProjectiveSpace(ZZ,1)
            sage: H = End(PS)
            sage: f = H([6*x^2+12*x*y+7*y^2, 12*x*y + 42*y^2])
            sage: g,M=f.minimal_model(return_transformation=True)
            sage: f.conjugate(M) == g
            True

        ::

            sage: PS.<X,Y> = ProjectiveSpace(QQ,1)
            sage: H = End(PS)
            sage: f = H([X+Y, X-3*Y])
            sage: f.minimal_model()
            Traceback (most recent call last):
            ...
            NotImplementedError: minimality is only for degree 2 or higher

        ::

            sage: PS.<X,Y> = ProjectiveSpace(QQ,1)
            sage: H = End(PS)
            sage: f = H([X^2-Y^2, X^2+X*Y])
            sage: f.minimal_model()
            Traceback (most recent call last):
            ...
            TypeError: the function is not a morphism

        """
        if self.base_ring() != QQ and self.base_ring() != ZZ:
            raise NotImplementedError("minimal models only implemented over ZZ or QQ")
        if not self.is_morphism():
            raise TypeError("the function is not a morphism")
        if self.degree() == 1:
            raise NotImplementedError("minimality is only for degree 2 or higher")

        from .endPN_minimal_model import affine_minimal
        return(affine_minimal(self, return_transformation, prime_list, False))

    def automorphism_group(self, **kwds):
        r"""
        Calculates the subgroup of `PGL2` that is the automorphism group of this map.

        The automorphism group is the set of `PGL(2)` elements that fix this map
        under conjugation.

        INPUT:

        keywords:

        - ``starting_prime`` -- The first prime to use for CRT. default: 5.(optional)

        - ``algorithm``-- Choose ``CRT``-Chinese Remainder Theorem- or ``fixed_points`` algorithm.
            default: depends on this map. (optional)

        - ``return_functions``-- Boolean - True returns elements as linear fractional transformations.
            False returns elements as `PGL2` matrices. default: False. (optional)

        - ``iso_type`` -- Boolean - True returns the isomorphism type of the automorphism group.
            default: False (optional)

        OUTPUT:

        - list - elements of automorphism group.

        AUTHORS:

        - Original algorithm written by Xander Faber, Michelle Manes, Bianca Viray

        - Modified by Joao Alberto de Faria, Ben Hutz, Bianca Thompson

        REFERENCES:

        .. [FMV] Xander Faber, Michelle Manes, and Bianca Viray. Computing Conjugating Sets
            and Automorphism Groups of Rational Functions. Journal of Algebra, 423 (2014), 1161-1190.

        EXAMPLES::

            sage: R.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = End(R)
            sage: f = H([x^2-y^2, x*y])
            sage: f.automorphism_group(return_functions=True)
            [x, -x]

        ::

            sage: R.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = End(R)
            sage: f = H([x^2 + 5*x*y + 5*y^2, 5*x^2 + 5*x*y + y^2])
            sage: f.automorphism_group()
            [
            [1 0]  [0 2]
            [0 1], [2 0]
            ]

        ::

            sage: R.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = End(R)
            sage: f=H([x^2-2*x*y-2*y^2, -2*x^2-2*x*y+y^2])
            sage: f.automorphism_group(return_functions=True)
            [x, 2/(2*x), -x - 1, -2*x/(2*x + 2), (-x - 1)/x, -1/(x + 1)]

        ::

            sage: R.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = End(R)
            sage: f = H([3*x^2*y - y^3, x^3 - 3*x*y^2])
            sage: f.automorphism_group(algorithm='CRT', return_functions=True, iso_type=True)
            ([x, (x + 1)/(x - 1), (-x + 1)/(x + 1), -x, 1/x, -1/x,
            (x - 1)/(x + 1), (-x - 1)/(x - 1)], 'Dihedral of order 8')

        ::

            sage: A.<z> = AffineSpace(QQ,1)
            sage: H = End(A)
            sage: f = H([1/z^3])
            sage: F = f.homogenize(1)
            sage: F.automorphism_group()
            [
            [1 0]  [0 2]  [-1  0]  [ 0 -2]
            [0 1], [2 0], [ 0  1], [ 2  0]
            ]

        ::

            sage: P.<x,y,z> = ProjectiveSpace(QQ,2)
            sage: H = End(P)
            sage: f = H([x**2 + x*z, y**2, z**2])
            sage: f.automorphism_group() # long test
            [
            [1 0 0]
            [0 1 0]
            [0 0 1]
            ]

        ::

            sage: K.<w> = CyclotomicField(3)
            sage: P.<x,y> = ProjectiveSpace(K, 1)
            sage: H = End(P)
            sage: D6 = H([y^2,x^2])
            sage: D6.automorphism_group()
            [
            [1 0]  [0 w]  [0 1]  [w 0]  [-w - 1      0]  [     0 -w - 1]
            [0 1], [1 0], [1 0], [0 1], [     0      1], [     1      0]
            ]
        """

        alg = kwds.get('algorithm', None)
        p = kwds.get('starting_prime', 5)
        return_functions = kwds.get('return_functions', False)
        iso_type = kwds.get('iso_type', False)
        if self.domain().dimension_relative() != 1:
            return self.conjugating_set(self)
        if self.base_ring() != QQ  and self.base_ring != ZZ:
            return self.conjugating_set(self)
        f = self.dehomogenize(1)
        R = PolynomialRing(f.base_ring(),'x')
        if is_FractionFieldElement(f[0]):
            F = (f[0].numerator().univariate_polynomial(R))/f[0].denominator().univariate_polynomial(R)
        else:
            F = f[0].univariate_polynomial(R)
        from .endPN_automorphism_group import automorphism_group_QQ_CRT, automorphism_group_QQ_fixedpoints
        if alg is None:
            if self.degree() <= 12:
                return(automorphism_group_QQ_fixedpoints(F, return_functions, iso_type))
            return(automorphism_group_QQ_CRT(F, p, return_functions, iso_type))
        elif alg == 'CRT':
            return(automorphism_group_QQ_CRT(F, p, return_functions, iso_type))

        return(automorphism_group_QQ_fixedpoints(F, return_functions, iso_type))

    def wronskian_ideal(self):
        r"""
        Returns the ideal generated by the critical point locus.

        This is the vanishing of the maximal minors of the Jacobian matrix.
        Not implemented for subvarieties.

        OUTPUT: an ideal in the coordinate ring of the domain of this map.

        Examples::

            sage: R.<x> = PolynomialRing(QQ)
            sage: K.<w> = NumberField(x^2+11)
            sage: P.<x,y> = ProjectiveSpace(K,1)
            sage: H = End(P)
            sage: f = H([x^2-w*y^2, w*y^2])
            sage: f.wronskian_ideal()
            Ideal ((4*w)*x*y) of Multivariate Polynomial Ring in x, y over Number
            Field in w with defining polynomial x^2 + 11

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: P2.<u,v,t> = ProjectiveSpace(K,2)
            sage: H = Hom(P,P2)
            sage: f = H([x^2-2*y^2, y^2, x*y])
            sage: f.wronskian_ideal()
            Ideal (4*x*y, 2*x^2 + 4*y^2, -2*y^2) of Multivariate Polynomial Ring in
            x, y over Rational Field
        """
        dom = self.domain()
        from sage.schemes.projective.projective_space import is_ProjectiveSpace
        if not (is_ProjectiveSpace(dom) and is_ProjectiveSpace(self.codomain())):
            raise NotImplementedError("not implemented for subschemes")
        N = dom.dimension_relative()+1
        R = dom.coordinate_ring()
        J = jacobian(self.defining_polynomials(),dom.gens())
        return(R.ideal(J.minors(N)))

    def critical_subscheme(self):
        r"""
        Returns the critical subscheme of this endomorphism defined over the base ring of this map.


        OUTPUT: the critical subscheme of the endomorphism.

        Examples::

            sage: set_verbose(None)
            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = End(P)
            sage: f = H([x^3-2*x*y^2 + 2*y^3, y^3])
            sage: f.critical_subscheme()
            Closed subscheme of Projective Space of dimension 1 over Rational Field
            defined by:
            9*x^2*y^2 - 6*y^4

        ::

            sage: set_verbose(None)
            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = End(P)
            sage: f = H([2*x^2-y^2, x*y])
            sage: f.critical_subscheme()
            Closed subscheme of Projective Space of dimension 1 over Rational Field
            defined by:
            4*x^2 + 2*y^2

        ::

            sage: P.<x,y,z> = ProjectiveSpace(QQ,2)
            sage: H = End(P)
            sage: f = H([2*x^2-y^2, x*y, z^2])
            sage: f.critical_subscheme()
            Closed subscheme of Projective Space of dimension 2 over Rational Field
            defined by:
            8*x^2*z + 4*y^2*z

        ::

            sage: P.<x,y,z,w> = ProjectiveSpace(GF(81),3)
            sage: H = End(P)
            sage: g = H([x^3+y^3, y^3+z^3, z^3+x^3, w^3])
            sage: g.critical_subscheme()
            Closed subscheme of Projective Space of dimension 3 over Finite Field in
            z4 of size 3^4 defined by:
              0

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = End(P)
            sage: f = H([x^2,x*y])
            sage: f.critical_subscheme()
            Traceback (most recent call last):
            ...
            TypeError: the function is not a morphism
        """
        from sage.schemes.projective.projective_space import is_ProjectiveSpace
        PS = self.domain()
        if not is_ProjectiveSpace(PS):
            raise NotImplementedError("not implemented for subschemes")
        if not self.is_endomorphism():
            raise NotImplementedError("must be an endomorphism")
        if not self.is_morphism():
            raise TypeError("the function is not a morphism")
        wr = self.wronskian_ideal()
        crit_subscheme = self.codomain().subscheme(wr)
        return crit_subscheme

    def critical_points(self, R=None):
        r"""
        Returns the critical points of this endomorphism defined over the ring ``R``
        or the base ring of this map.

        Must be dimension 1.

        INPUT:

            - ``R`` - a ring (optional).

        OUTPUT: a list of projective space points defined over ``R``.

        Examples::

            sage: set_verbose(None)
            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = End(P)
            sage: f = H([x^3-2*x*y^2 + 2*y^3, y^3])
            sage: f.critical_points()
            [(1 : 0)]
            sage: K.<w> = QuadraticField(6)
            sage: f.critical_points(K)
            [(-1/3*w : 1), (1/3*w : 1), (1 : 0)]

        ::

            sage: set_verbose(None)
            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = End(P)
            sage: f = H([2*x^2-y^2, x*y])
            sage: f.critical_points(QQbar)
            [(-0.7071067811865475?*I : 1), (0.7071067811865475?*I : 1)]
        """
        from sage.schemes.projective.projective_space import is_ProjectiveSpace
        PS = self.domain()
        if PS.dimension_relative() > 1:
            raise NotImplementedError("use .wronskian_ideal() for dimension > 1")
        if R is None:
            F = self
        else:
            F = self.change_ring(R)
        P = F.codomain()
        X = F.critical_subscheme()
        crit_points = [P(Q) for Q in X.rational_points()]
        return crit_points

    def is_postcritically_finite(self, err=0.01, embedding=None):
        r"""
        Determine if this map is post-critically finite.

        Only for endomorphisms of `\mathbb{P}^1`. It checks if each critical point
        is preperiodic. The optional parameter ``err`` is passed into
        ``is_preperiodic()`` as part of the preperiodic check.

        INPUT:

            - ``err`` - positive real number (optional, Default: 0.01).

            - ``embedding`` - embedding of base ring into `\QQbar`

        OUTPUT: Boolean

        Examples::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = End(P)
            sage: f = H([x^2 - y^2, y^2])
            sage: f.is_postcritically_finite()
            True

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = End(P)
            sage: f = H([x^3- y^3, y^3])
            sage: f.is_postcritically_finite()
            False

        ::

            sage: R.<z> = QQ[]
            sage: K.<v> = NumberField(z^8 + 3*z^6 + 3*z^4 + z^2 + 1)
            sage: PS.<x,y> = ProjectiveSpace(K,1)
            sage: H = End(PS)
            sage: f = H([x^3+v*y^3, y^3])
            sage: f.is_postcritically_finite(embedding=K.embeddings(QQbar)[0]) # long time
            True

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = End(P)
            sage: f = H([6*x^2+16*x*y+16*y^2, -3*x^2-4*x*y-4*y^2])
            sage: f.is_postcritically_finite()
            True
        """
        if not self.is_endomorphism():
            raise TypeError("must be an endomorphism")

        #iteration of subschemes not yet implemented
        if self.domain().dimension_relative() > 1:
            raise NotImplementedError("only implemented in dimension 1")

        #Since is_preperiodic uses heights we need to be over a numberfield
        K = FractionField(self.codomain().base_ring())
        if not K in _NumberFields and not K is QQbar:
            raise NotImplementedError("must be over a number field or a number field order or QQbar")

        if embedding is None:
            F = self.change_ring(QQbar)
        else:
            F = self.change_ring(embedding)
        crit_points = F.critical_points()
        pcf = True
        i = 0
        while pcf and i < len(crit_points):
            if crit_points[i].is_preperiodic(F, err) == False:
                pcf = False
            i += 1
        return(pcf)

    def critical_point_portrait(self, check=True, embedding=None):
        r"""
        If this map is post-critically finite, return its critical point portrait.

        This is the directed graph of iterates starting with the critical points. Must be
        dimension 1. If ``check`` is True, then the map is first checked to see if it is
        postcritically finite.

        INPUT:

            - check - Boolean.

            - ``embedding`` - embedding of base ring into `\QQbar`

        OUTPUT: a digraph.

        Examples::

            sage: R.<z> = QQ[]
            sage: K.<v> = NumberField(z^6 + 2*z^5 + 2*z^4 + 2*z^3 + z^2 + 1)
            sage: PS.<x,y> = ProjectiveSpace(K,1)
            sage: H = End(PS)
            sage: f = H([x^2+v*y^2, y^2])
            sage: f.critical_point_portrait(check=False, embedding=K.embeddings(QQbar)[0]) # long time
            Looped digraph on 6 vertices

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = End(P)
            sage: f = H([x^5 + 5/4*x*y^4, y^5])
            sage: f.critical_point_portrait(check=False)
            Looped digraph on 5 vertices

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = End(P)
            sage: f = H([x^2 + 2*y^2, y^2])
            sage: f.critical_point_portrait()
            Traceback (most recent call last):
            ...
            TypeError: map be be post-critically finite
        """
        #input checking done in is_postcritically_finite
        if check:
            if not self.is_postcritically_finite():
                raise TypeError("map be be post-critically finite")
        if embedding is None:
            F = self.change_ring(QQbar)
        else:
            F = self.change_ring(embedding)
        crit_points = F.critical_points()
        N = len(crit_points)
        for i in range(N):
            done = False
            Q= F(crit_points[i])
            while not done:
                if Q in crit_points:
                    done = True
                else:
                    crit_points.append(Q)
                Q = F(Q)
        return(F._preperiodic_points_to_cyclegraph(crit_points))

    def critical_height(self, **kwds):
        r"""
        Compute the critical height of this map.

        The critical height is defined by J. Silverman as
        the sum of the canonical heights of the critical points. This must be dimension 1 and
        defined over a number field or number field order.

        INPUT:

        kwds:

        - ``badprimes`` - a list of primes of bad reduction. (optional)

        - ``N`` - positive integer. number of terms of the series to use in the local green functions.
          (optional - Default: 10)

        - ``prec`` - positive integer, float point or p-adic precision, Default: 100.

        - ``error_bound`` - a positive real number. (optional)

        - ``embedding`` - the embedding of the base field to `\QQbar` (optional)

        OUTPUT: Real number

        Examples::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = End(P)
            sage: f = H([x^3+7*y^3, 11*y^3])
            sage: f.critical_height()
            1.1989273321156851418802151128

        ::

            sage: K.<w> = QuadraticField(2)
            sage: P.<x,y> = ProjectiveSpace(K,1)
            sage: H = Hom(P,P)
            sage: f = H([x^2+w*y^2, y^2])
            sage: f.critical_height()
            0.16090842452312941163719755472

        Postcritically finite maps have critical height 0::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = End(P)
            sage: f = H([x^3-3/4*x*y^2 + 3/4*y^3, y^3])
            sage: f.critical_height(error_bound=0.0001)
            0.000011738508366948556443245983996
        """
        if not self.is_endomorphism():
            raise TypeError("must be an endomorphism")
        PS = self.codomain()
        if PS.dimension_relative() > 1:
            raise NotImplementedError("only implemented in dimension 1")

        K = FractionField(PS.base_ring())
        if not K in _NumberFields and not K is QQbar:
            raise NotImplementedError("must be over a number field or a number field order or QQbar")
        #doesn't really matter which we choose as Galois conjugates have the same height
        emb = kwds.get("embedding", K.embeddings(QQbar)[0])
        F = self.change_ring(K).change_ring(emb)
        crit_points = F.critical_points()
        n = len(crit_points)
        err_bound = kwds.get("error_bound", None)
        if not err_bound is None:
            kwds.update({"error_bound": err_bound/n})
        ch = 0
        for P in crit_points:
            ch += P.canonical_height(F, **kwds)
        return(ch)

    def periodic_points(self, n, minimal=True, R=None, algorithm='variety'):
        r"""
        Computes the periodic points of period ``n`` of this map
        defined over the ring ``R`` or the base ring of the map.

        This can be done either by finding the rational points on the variety
        defining the points of period ``n``, or, for finite fields,
        finding the cycle of appropriate length in the cyclegraph. For small
        cardinality fields, the cyclegraph algorithm is effective for any
        map and length cycle, but is slow when the cyclegraph is large.
        The variety algorithm is good for small period, degree, and dimension,
        but is slow as the defining equations of the variety get more
        complicated.

        The map must be a projective morphism.

        INPUT:

        - ``n`` - a positive integer.

        - ``minimal`` - Boolean. True specifies to find only the periodic points of minimal period ``n``.
            False specifies to find all periodic points of period ``n``. Default: True.

        - ``R`` a commutative ring.

        - ``algorithm`` - a string. Either ``variety`` to find the rational points
          on the appropriate variety or ``cyclegraph`` to find the cycles from the
          cycle graph. Default: ``variety``.

        OUTPUT:

        - a list of periodic points of this map.

        EXAMPLES::

            sage: set_verbose(None)
            sage: P.<x,y> = ProjectiveSpace(QQbar,1)
            sage: H = Hom(P,P)
            sage: f = H([x^2-x*y+y^2, x^2-y^2+x*y])
            sage: f.periodic_points(1)
            [(-0.500000000000000? - 0.866025403784439?*I : 1), (-0.500000000000000? + 0.866025403784439?*I : 1),
            (1 : 1)]

        ::

            sage: P.<x,y,z> = ProjectiveSpace(QuadraticField(5,'t'),2)
            sage: H = Hom(P,P)
            sage: f = H([x^2 - 21/16*z^2, y^2-z^2, z^2])
            sage: f.periodic_points(2)
            [(-5/4 : -1 : 1), (-5/4 : -1/2*t + 1/2 : 1), (-5/4 : 0 : 1), (-5/4 : 1/2*t + 1/2 : 1), (-3/4 : -1 : 1),
            (-3/4 : 0 : 1), (1/4 : -1 : 1), (1/4 : -1/2*t + 1/2 : 1), (1/4 : 0 : 1), (1/4 : 1/2*t + 1/2 : 1),
            (7/4 : -1 : 1), (7/4 : 0 : 1)]

        ::

            sage: w = QQ['w'].0
            sage: K = NumberField(w^6 - 3*w^5 + 5*w^4 - 5*w^3 + 5*w^2 - 3*w + 1,'s')
            sage: P.<x,y,z> = ProjectiveSpace(K,2)
            sage: H = Hom(P,P)
            sage: f = H([x^2+z^2, y^2+x^2, z^2+y^2])
            sage: f.periodic_points(1)
            [(-2*s^5 + 4*s^4 - 5*s^3 + 3*s^2 - 4*s : -2*s^5 + 5*s^4 - 7*s^3 + 6*s^2 - 7*s + 3 : 1),
            (-s^5 + 3*s^4 - 4*s^3 + 4*s^2 - 4*s + 2 : -s^5 + 2*s^4 - 2*s^3 + s^2 - s : 1),
            (2*s^5 - 6*s^4 + 9*s^3 - 8*s^2 + 7*s - 4 : 2*s^5 - 5*s^4 + 7*s^3 - 5*s^2 + 6*s - 2 : 1),
            (s^5 - 2*s^4 + 2*s^3 + s : s^5 - 3*s^4 + 4*s^3 - 3*s^2 + 2*s - 1 : 1),
            (s^5 - 2*s^4 + 3*s^3 - 3*s^2 + 3*s - 1 : -s^5 + 3*s^4 - 5*s^3 + 4*s^2 - 4*s + 2 : 1), (1 : 1 : 1),
            (-s^5 + 3*s^4 - 5*s^3 + 4*s^2 - 3*s + 1 : s^5 - 2*s^4 + 3*s^3 - 3*s^2 + 4*s - 1 : 1)]

        ::

            sage: P.<x,y,z> = ProjectiveSpace(QQ,2)
            sage: H = Hom(P,P)
            sage: f = H([x^2 - 21/16*z^2, y^2-2*z^2, z^2])
            sage: f.periodic_points(2, False)
            [(-5/4 : -1 : 1), (-5/4 : 2 : 1), (-3/4 : -1 : 1), (-3/4 : 2 : 1), (0 : 1 : 0), (1/4 : -1 : 1),
            (1/4 : 2 : 1), (1 : 0 : 0), (1 : 1 : 0), (7/4 : -1 : 1), (7/4 : 2 : 1)]

        ::

            sage: P.<x,y,z> = ProjectiveSpace(QQ,2)
            sage: H = Hom(P,P)
            sage: f = H([x^2 - 21/16*z^2, y^2-2*z^2, z^2])
            sage: f.periodic_points(2)
            [(-5/4 : -1 : 1), (-5/4 : 2 : 1), (1/4 : -1 : 1), (1/4 : 2 : 1)]

        ::

            sage: set_verbose(None)
            sage: P.<x,y> = ProjectiveSpace(ZZ, 1)
            sage: H = End(P)
            sage: f = H([x^2+y^2,y^2])
            sage: f.periodic_points(2, R=QQbar, minimal=False)
            [(-0.500000000000000? - 1.322875655532296?*I : 1),
             (-0.500000000000000? + 1.322875655532296?*I : 1),
             (0.500000000000000? - 0.866025403784439?*I : 1),
             (0.500000000000000? + 0.866025403784439?*I : 1),
             (1 : 0)]

        ::

            sage: P.<x,y> = ProjectiveSpace(GF(307), 1)
            sage: H = End(P)
            sage: f = H([x^10+y^10, y^10])
            sage: f.periodic_points(16, minimal=True, algorithm='cyclegraph')
            [(69 : 1), (185 : 1), (120 : 1), (136 : 1), (97 : 1), (183 : 1),
             (170 : 1), (105 : 1), (274 : 1), (275 : 1), (154 : 1), (156 : 1),
             (87 : 1), (95 : 1), (161 : 1), (128 : 1)]

        ::

            sage: P.<x,y> = ProjectiveSpace(GF(13^2,'t'),1)
            sage: H = End(P)
            sage: f = H([x^3 + 3*y^3, x^2*y])
            sage: f.periodic_points(30, minimal=True, algorithm='cyclegraph')
            [(t + 3 : 1), (6*t + 6 : 1), (7*t + 1 : 1), (2*t + 8 : 1),
             (3*t + 4 : 1), (10*t + 12 : 1), (8*t + 10 : 1), (5*t + 11 : 1),
             (7*t + 4 : 1), (4*t + 8 : 1), (9*t + 1 : 1), (2*t + 2 : 1),
             (11*t + 9 : 1), (5*t + 7 : 1), (t + 10 : 1), (12*t + 4 : 1),
             (7*t + 12 : 1), (6*t + 8 : 1), (11*t + 10 : 1), (10*t + 7 : 1),
             (3*t + 9 : 1), (5*t + 5 : 1), (8*t + 3 : 1), (6*t + 11 : 1),
             (9*t + 12 : 1), (4*t + 10 : 1), (11*t + 4 : 1), (2*t + 7 : 1),
             (8*t + 12 : 1), (12*t + 11 : 1)]

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = End(P)
            sage: f = H([3*x^2+5*y^2,y^2])
            sage: f.periodic_points(2, R=GF(3), minimal=False)
            Traceback (most recent call last):
            ...
            NotImplementedError: must be a projective morphism
        """
        if n <= 0:
            raise ValueError("a positive integer period must be specified")
        if not self.is_endomorphism():
            raise TypeError("must be an endomorphism")
        if R is None:
            f = self
            R = self.base_ring()
        else:
            f = self.change_ring(R)
        if not f.is_morphism():
            # if not, the variety is not dimension 0 and
            # can cannot construct the cyclegraph due to
            #indeterminacies
            raise NotImplementedError("must be a projective morphism")
        PS = f.codomain()
        if algorithm == 'variety':
            if R in NumberFields() or R is QQbar or R in FiniteFields():
                N = PS.dimension_relative() + 1
                R = PS.coordinate_ring()
                F = f.nth_iterate_map(n)
                L = [F[i]*R.gen(j) - F[j]*R.gen(i) for i in range(0,N) for j in range(i+1, N)]
                X = PS.subscheme(L)
                points = [PS(Q) for Q in X.rational_points()]

                if not minimal:
                    return points
                else:
                    #we want only the points with minimal period n
                    #so we go through the list and remove any that
                    #have smaller period by checking the iterates
                    rem_indices = []
                    for i in range(len(points)-1,-1,-1):
                        # iterate points to check if minimal
                        P = points[i]
                        for j in range(1,n):
                            P = f(P)
                            if P == points[i]:
                                points.pop(i)
                                break
                    return points
            else:
                raise NotImplementedError("ring must a number field or finite field")
        elif algorithm == 'cyclegraph':
            if R in FiniteFields():
                g = f.cyclegraph()
                points = []
                for cycle in g.all_simple_cycles():
                    m = len(cycle)-1
                    if minimal:
                        if m == n:
                            points = points + cycle[:-1]
                    else:
                        if n % m == 0:
                            points = points + cycle[:-1]
                return(points)
            else:
                raise TypeError("ring must be finite to generate cyclegraph")
        else:
            raise ValueError("algorithm must be either 'variety' or 'cyclegraph'")

    def multiplier_spectra(self, n, formal=True, embedding=None):
        r"""
        Computes the formal ``n`` multiplier spectra of this map.

        This is the set of multipliers of the periodic points of formal period
        ``n`` included with the appropriate multiplicity.
        User can also specify to compute the ``n`` multiplier spectra instead which includes the
        multipliers of all periodic points of period ``n``.The map must be defined over
        projective space of dimension 1 over a number field.

        INPUT:

        - ``n`` - a positive integer, the period.

        - ``formal`` - a Boolean. True specifies to find the formal ``n`` multiplier spectra
            of this map. False specifies to find the ``n`` multiplier spectra
            of this map. Default: True

        - ``embedding`` - embedding of the base field into `\QQbar`

        OUTPUT:

        - a list of `\QQbar` elements.

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = End(P)
            sage: f = H([4608*x^10 - 2910096*x^9*y + 325988068*x^8*y^2 + 31825198932*x^7*y^3 - 4139806626613*x^6*y^4\
            - 44439736715486*x^5*y^5 + 2317935971590902*x^4*y^6 - 15344764859590852*x^3*y^7 + 2561851642765275*x^2*y^8\
            + 113578270285012470*x*y^9 - 150049940203963800*y^10, 4608*y^10])
            sage: f.multiplier_spectra(1)
            [0, -7198147681176255644585/256, 848446157556848459363/19683, -3323781962860268721722583135/35184372088832,
            529278480109921/256, -4290991994944936653/2097152, 1061953534167447403/19683, -3086380435599991/9,
            82911372672808161930567/8192, -119820502365680843999, 3553497751559301575157261317/8192]

        ::

            sage: set_verbose(None)
            sage: z = QQ['z'].0
            sage: K.<w> = NumberField(z^4 - 4*z^2 + 1,'z')
            sage: P.<x,y> = ProjectiveSpace(K,1)
            sage: H = End(P)
            sage: f = H([x^2 - w/4*y^2, y^2])
            sage: f.multiplier_spectra(2, False, embedding=K.embeddings(QQbar)[0])
            [0,
             5.931851652578137? + 0.?e-49*I,
             0.0681483474218635? - 1.930649271699173?*I,
             0.0681483474218635? + 1.930649271699173?*I]

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = End(P)
            sage: f = H([x^2 - 3/4*y^2, y^2])
            sage: f.multiplier_spectra(2)
            [1]

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = End(P)
            sage: f = H([x^2 - 7/4*y^2, y^2])
            sage: f.multiplier_spectra(3)
            [1, 1]

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = End(P)
            sage: f = H([x^2 + y^2, x*y])
            sage: f.sigma_invariants(1)
            [3, 3, 1]
        """
        PS = self.domain()
        n = Integer(n)

        if (n < 1):
            raise ValueError("period must be a positive integer")
        from sage.schemes.projective.projective_space import is_ProjectiveSpace
        if not is_ProjectiveSpace(PS):
            raise NotImplementedError("not implemented for subschemes")
        if (PS.dimension_relative() > 1):
            raise NotImplementedError("only implemented for dimension 1")
        if not self.is_endomorphism():
            raise TypeError("self must be an endomorphism")
        if not PS.base_ring() in NumberFields() and not PS.base_ring() is QQbar:
            raise NotImplementedError("self must be a map over a number field")

        if embedding is None:
            f = self.change_ring(QQbar)
        else:
            f = self.change_ring(embedding)

        PS = f.domain()

        if not formal:
            G = f.nth_iterate_map(n)
            F = G[0]*PS.gens()[1] - G[1]*PS.gens()[0]
        else:
            # periodic points of formal period n are the roots of the nth dynatomic polynomial
            K = f._number_field_from_algebraics()
            F = K.dynatomic_polynomial(n)
            if K.domain().base_ring() != QQ: # need to coerce F to poly over QQbar. This is only needed if base ring is not QQ
                abspoly = K.domain().base_ring().absolute_polynomial()
                phi = K.domain().base_ring().hom(QQbar.polynomial_root(abspoly,abspoly.any_root(CIF)),QQbar)
                Kx = K.coordinate_ring()
                QQbarx = QQbar[Kx.variable_names()]
                phix = Kx.hom(phi,QQbarx)
                F = phix(F)

        other_roots = F.parent()(F([(f.domain().gens()[0]),1])).univariate_polynomial().roots(ring=QQbar)

        points = []

        minfty = min([e[1] for e in F.exponents()]) # include the point at infinity with the right multiplicity
        for i in range(minfty):
            points.append(PS([1,0]))

        for R in other_roots:
            for i in range(R[1]):
                points.append(PS([R[0],1])) # include copies of higher multiplicity roots

        newpoints = [] # should include one representative point per cycle, included with the right multiplicity

        while(points != []):
            P = points[0]
            newpoints.append(P)
            points.pop(0)
            Q = P
            for i in range(1,n):
                try:
                    points.remove(f(Q))
                except ValueError:
                    pass
                Q = f(Q)

        multipliers = [f.multiplier(P,n)[0,0] for P in newpoints]

        return multipliers

    def sigma_invariants(self, n, formal=True, embedding=None):
        r"""
        Computes the values of the elementary symmetric polynomials of the formal ``n`` multilpier spectra
        of this map.

        Can specify to instead compute the values corresponding to the elementary symmetric
        polynomials of the ``n`` multiplier spectra, which includes the multipliers of all periodic
        points of period ``n``. The map must be defined over projective space of dimension 1 over
        a number field.

        INPUT:

        - ``n`` - a positive integer, the period.

        - ``formal`` - a Boolean. True specifies to find the values of the elementary symmetric polynomials
            corresponding to the formal ``n`` multiplier spectra of this map. False specifies to instead find
            the values corresponding to the ``n`` multiplier spectra of this map, which includes the multipliers
            of all periodic points of period ``n`` of this map. Default: True

        - ``embedding`` - embedding of the base field into `\QQbar`

        OUTPUT: a list of elements in the base ring.

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = End(P)
            sage: f = H([512*x^5 - 378128*x^4*y + 76594292*x^3*y^2 - 4570550136*x^2*y^3 - 2630045017*x*y^4\
            + 28193217129*y^5, 512*y^5])
            sage: f.sigma_invariants(1)
            [19575526074450617/1048576, -9078122048145044298567432325/2147483648,
            -2622661114909099878224381377917540931367/1099511627776,
            -2622661107937102104196133701280271632423/549755813888,
            338523204830161116503153209450763500631714178825448006778305/72057594037927936, 0]

        ::

            sage: set_verbose(None)
            sage: z = QQ['z'].0
            sage: K = NumberField(z^4 - 4*z^2 + 1,'z')
            sage: P.<x,y> = ProjectiveSpace(K,1)
            sage: H = End(P)
            sage: f = H([x^2 - 5/4*y^2, y^2])
            sage: f.sigma_invariants(2, False, embedding=K.embeddings(QQbar)[0])
            [13, 11, -25, 0]

        ::

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: H = End(P)
            sage: f = H([x^2, y^2, z^2])
            sage: f.sigma_invariants(1)
            Traceback (most recent call last):
            ...
            NotImplementedError: only implemented for dimension 1
        """
        polys = []

        multipliers = self.multiplier_spectra(n, formal, embedding=embedding)

        e = SymmetricFunctions(QQbar).e()

        N = len(multipliers)
        R = self.base_ring()
        for i in range(0,N):
            polys.append(R(e([i+1]).expand(N)(multipliers)))
        return polys

    def reduced_form(self, prec=300, return_conjugation=True, error_limit=0.000001):
        r"""
        Returns reduced form of this projective morphism.

        The reduced form is the `SL(2, \ZZ)` equivalent morphism obtained by applying
        the binary form reduction algorithm from Stoll and Cremona [SC]_
        to the homogeneous polynomial defining the periodic points (the dynatomic polynomial).
        The smallest period `n` with enough periodic points is used.

        This should also minimize the sum of the squares of the coefficients,
        but this is not always the case.

        See :meth:`sage.rings.polynomial.multi_polynomial.reduced_form` for the information
        on binary form reduction.

        Implemented by Rebecca Lauren Miller as part of GSOC 2016.

        INPUT:

        - ``prec`` -- integer, desired precision (default: 300)

        - ``return_conjuagtion`` -- A Boolean. Returns element of `SL(2, \ZZ)`. (default: True)

        - ``error_limit`` -- sets the error tolerance (default:0.000001)

        OUTPUT:

        - a projective morphism

        - a matrix

        EXAMPLES::

            sage: PS.<x,y> = ProjectiveSpace(QQ, 1)
            sage: H = End(PS)
            sage: f = H([x^3 + x*y^2, y^3])
            sage: m = matrix(QQ, 2, 2, [-221, -1, 1, 0])
            sage: f = f.conjugate(m)
            sage: f.reduced_form(prec=100) #needs 2 periodic
            Traceback (most recent call last):
            ...
            ValueError: not enough precision
            sage: f.reduced_form()
            (
            Scheme endomorphism of Projective Space of dimension 1 over Rational Field
            Defn: Defined on coordinates by sending (x : y) to
            (x^3 + x*y^2 : y^3)
            ,
            [  0  -1]
            [  1 221]
            )

        ::

            sage: PS.<x,y> = ProjectiveSpace(ZZ, 1)
            sage: H = End(PS)
            sage: f = H([x^2+ x*y, y^2]) #needs 3 periodic
            sage: m = matrix(QQ, 2, 2, [-221, -1, 1, 0])
            sage: f = f.conjugate(m)
            sage: f.reduced_form(prec=200)
            (
            Scheme endomorphism of Projective Space of dimension 1 over Integer Ring
            Defn: Defined on coordinates by sending (x : y) to
            (-x^2 + x*y - y^2 : -y^2)
            ,
            [  0  -1]
            [  1 220]
            )

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: H = End(P)
            sage: f = H([x^3, y^3])
            sage: f.reduced_form()
            (
            Scheme endomorphism of Projective Space of dimension 1 over Rational Field
            Defn: Defined on coordinates by sending (x : y) to
            (x^3 : y^3)
            ,
            [-1  0]
            [ 0 -1]
            )

        ::

            sage: PS.<X,Y> = ProjectiveSpace(QQ,1)
            sage: H = End(PS)
            sage: f = H([7365*X^4 + 12564*X^3*Y + 8046*X^2*Y^2 + 2292*X*Y^3 + 245*Y^4,\
            -12329*X^4 - 21012*X^3*Y - 13446*X^2*Y^2 - 3828*X*Y^3 - 409*Y^4])
            sage: f.reduced_form(prec=30)
            Traceback (most recent call last):
            ...
            ValueError: accuracy of Newton's root not within tolerance(1.2519607 > 1e-06), increase precision
            sage: f.reduced_form()
            (
            Scheme endomorphism of Projective Space of dimension 1 over Rational Field
            Defn: Defined on coordinates by sending (X : Y) to
            (-7*X^4 - 12*X^3*Y - 42*X^2*Y^2 - 12*X*Y^3 - 7*Y^4 : -X^4 - 4*X^3*Y - 6*X^2*Y^2 - 4*X*Y^3 - Y^4),
            <BLANKLINE>
            [-1  2]
            [ 2 -5]
            )

        ::

            sage: P.<x,y> = ProjectiveSpace(RR, 1)
            sage: H = End(P)
            sage: f = H([x^4, RR(sqrt(2))*y^4])
            sage: m = matrix(RR, 2, 2, [1,12,0,1])
            sage: f = f.conjugate(m)
            sage: f.reduced_form()
            (
            Scheme endomorphism of Projective Space of dimension 1 over Real Field with 53 bits of precision
              Defn: Defined on coordinates by sending (x : y) to
                    (-x^4 + 2.86348722511320e-12*y^4 : -1.41421356237310*y^4)
            ,
            [-1 12]
            [ 0 -1]
            )

        ::

            sage: P.<x,y> = ProjectiveSpace(CC, 1)
            sage: H = End(P)
            sage: f = H([x^4, CC(sqrt(-2))*y^4])
            sage: m = matrix(CC, 2, 2, [1,12,0,1])
            sage: f = f.conjugate(m)
            sage: f.reduced_form()
            (
            Scheme endomorphism of Projective Space of dimension 1 over Complex Field with 53 bits of precision
              Defn: Defined on coordinates by sending (x : y) to
                    (-x^4 + (-1.03914726748259e-15)*y^4 : (-8.65956056235493e-17 - 1.41421356237309*I)*y^4)
            ,
            [-1 12]
            [ 0 -1]
            )

        ::

            sage: K.<w> = QuadraticField(2)
            sage: P.<x,y> = ProjectiveSpace(K, 1)
            sage: H = End(P)
            sage: f = H([x^3, w*y^3])
            sage: m = matrix(K, 2, 2, [1,12,0,1])
            sage: f = f.conjugate(m)
            sage: f.reduced_form()
            (
             Scheme endomorphism of Projective Space of dimension 1 over Number Field in w
             with defining polynomial x^2 - 2
              Defn: Defined on coordinates by sending (x : y) to
                    (x^3 : (w)*y^3)
            ,
            [-1 12]
            [ 0 -1]
            )

        ::

            sage: R.<x> = QQ[]
            sage: K.<w> = NumberField(x^5+x-3, embedding=(x^5+x-3).roots(ring=CC)[0][0])
            sage: P.<x,y> = ProjectiveSpace(K, 1)
            sage: H = End(P)
            sage: f = H([12*x^3, 2334*w*y^3])
            sage: m = matrix(K, 2, 2, [-12,1,1,0])
            sage: f = f.conjugate(m)
            sage: f.reduced_form()
            (
            Scheme endomorphism of Projective Space of dimension 1 over Number Field
            in w with defining polynomial x^5 + x - 3
              Defn: Defined on coordinates by sending (x : y) to
                    (12*x^3 : (2334*w)*y^3)
            ,
            [  0  -1]
            [  1 -12]
            )
        """
        R = self.coordinate_ring()
        F = R(self.dynatomic_polynomial(1))
        x,y = R.gens()
        d = gcd(F, F.derivative(x)).degree() #counts multiple roots
        n = 2
        # Checks to make sure there are enough distinct, roots we need 3
        # if there are not it finds the nth periodic points until there are enough
        while F.degree()-d <= 2:
            F = self.dynatomic_polynomial(n) # finds n periodic points
            d = gcd(F, F.derivative(x)).degree() #counts multiple roots
            n += 1
        G,m = F.reduced_form(prec=prec, return_conjugation=return_conjugation)
        if return_conjugation:
            return (self.conjugate(m), m)
        return self.conjugate(m)

class SchemeMorphism_polynomial_projective_space_field(SchemeMorphism_polynomial_projective_space):

    def lift_to_rational_periodic(self, points_modp, B=None):
        r"""
        Given a list of points in projective space over `GF(p)`, determine if they lift to
        `\QQ`-rational periodic points.

        The map must be an endomorphism of projective space defined over `\QQ`

        ALGORITHM:
            Use Hensel lifting to find a `p`-adic approximation for that rational point. The accuracy needed
            is determined by the height bound ``B``. Then apply the LLL algorithm to determine if the lift
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
            sage: f = H([x^2 - y^2, y^2])
            sage: f.lift_to_rational_periodic([[P(0,1).change_ring(GF(7)), 4]])
            [[(0 : 1), 2]]

        ::

            There may be multiple points in the lift.
            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = End(P)
            sage: f = H([-5*x^2 + 4*y^2, 4*x*y])
            sage: f.lift_to_rational_periodic([[P(1,0).change_ring(GF(3)), 1]]) # long time
            [[(1 : 0), 1], [(2/3 : 1), 1], [(-2/3 : 1), 1]]

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = End(P)
            sage: f = H([16*x^2 - 29*y^2, 16*y^2])
            sage: f.lift_to_rational_periodic([[P(3,1).change_ring(GF(13)), 3]])
            [[(-1/4 : 1), 3]]

        ::

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: H = End(P)
            sage: f = H([76*x^2 - 180*x*y + 45*y^2 + 14*x*z + 45*y*z - 90*z^2, 67*x^2 - 180*x*y - 157*x*z + 90*y*z, -90*z^2])
            sage: f.lift_to_rational_periodic([[P(14,19,1).change_ring(GF(23)), 9]]) # long time
            [[(-9 : -4 : 1), 9]]
        """
        if points_modp == []:
            return([])
        else:
            if B is None:
                B = e ** self.height_difference_bound()

            p = points_modp[0][0].codomain().base_ring().characteristic()
            if p == 0:
                raise TypeError("must be positive characteristic")
            PS = self.domain()
            N = PS.dimension_relative()
            R = RealField()
            #compute the maximum p-adic precision needed to conclusively determine
            #if the rational point exists
            L = R((R(2 ** (N / 2 + 1) * sqrt(N + 1) * B ** 2).log()) / R(p).log() + 1).trunc()

            points = []
            for i in range(len(points_modp)):
                #[point mod p, period, current p-adic precision]
                points.append([points_modp[i][0].change_ring(QQ, check = False), points_modp[i][1], 1])
            good_points = []
            #shifts is used in non-Hensel lifting
            shifts = None
            #While there are still points to consider try to lift to next precision
            while points != []:
                q = points.pop()
                qindex = N
                #Find the last non-zero coordinate to use for normalizations
                while q[0][qindex] % p == 0:
                    qindex -= 1
                T = q[0]
                n = q[1]
                k = q[2]
                T.scale_by(1 / T[qindex]) #normalize
                bad = 0
                #stop where we reach the needed precision or the point is bad
                while k < L and bad == 0:
                    l = self._multipliermod(T, n, p, 2 * k)
                    l -= l.parent().one() #f^n(x) - x
                    lp = l.change_ring(Zmod(p ** k))
                    ldet = lp.determinant()
                    # if the matrix is invertible then we can Hensel lift
                    if ldet % p != 0:
                        RQ = ZZ.quo(p ** (2 * k))
                        T.clear_denominators()
                        newT = T.change_ring(RQ, check = False)
                        fp = self.change_ring(RQ, check = False)
                        S = newT.nth_iterate(fp, n, normalize = False).change_ring(QQ, check = False)
                        T.scale_by(1 / T[qindex])
                        S.scale_by(1 / S[qindex])
                        for i in range(N + 1):
                            S._coords[i] = S._coords[i] - T._coords[i]
                            if S[i] % (p ** k) != 0 and i != N:
                                bad = 1
                                break
                        if bad == 1:
                            break
                        S.scale_by(-1 / p ** k)
                        vecs = [Zmod(p ** k)(S._coords[iS]) for iS in range(N + 1)]
                        vecs.pop(qindex)
                        newvecs = list((lp.inverse()) * vector(vecs)) #l.inverse should be mod p^k!!
                        newS = []
                        [newS.append(QQ(newvecs[i])) for i in range(qindex)]
                        newS.append(0)
                        [newS.append(QQ(newvecs[i])) for i in range(qindex, N)]
                        S = PS.point(newS, False) #don't check for [0,...,0]
                        for i in range(N + 1):
                            S._coords[i] = S._coords[i] % (p ** k)
                        for i in range(N + 1):
                            T._coords[i] += S._coords[i] * (p ** k)
                        T.normalize_coordinates()
                        #Hensel gives us 2k for the newprecision
                        k = min(2 * k, L)
                    else:
                        #we are unable to Hensel Lift so must try all possible lifts
                        #to the next precision (k+1)
                        first = 0
                        newq = []
                        RQ = Zmod(p ** (k + 1))
                        fp = self.change_ring(RQ, check = False)
                        if shifts is None:
                            shifts = xmrange([p for i in range(N)])
                        for shift in shifts:
                            newT = T.change_ring(RQ, check = False)
                            shiftindex = 0
                            for i in range(N + 1):
                                if i != qindex:
                                    newT._coords[i] = newT[i] + shift[shiftindex] * p ** k
                                    shiftindex += 1
                            TT = fp.nth_iterate(newT, n, normalize = False)
                            if TT == newT:
                                if first == 0:
                                    newq.append(newT.change_ring(QQ, check = False))
                                    newq.append(n)
                                    newq.append(k + 1)
                                    first = 1
                                else:
                                    points.append([newT.change_ring(QQ, check = False), n, k + 1])
                        if newq == []:
                            bad = 1
                            break
                        else:
                            T = newq[0]
                            k += 1
                #given a p-adic lift of appropriate precision
                #perform LLL to find the "smallest" rational approximation
                #If this height is small enough, then it is a valid rational point
                if bad == 0:
                    M = matrix(N + 2, N + 1)
                    T.clear_denominators()
                    for i in range(N + 1):
                        M[0, i] = T[i]
                        M[i + 1, i] = p ** L
                    M[N + 1, N] = p ** L
                    M = M.LLL()
                    Q = []
                    [Q.append(M[1, i]) for i in range(N + 1)]
                    g = gcd(Q)
                    #remove gcds since this is a projective point
                    newB = B * g
                    for i in range(N + 1):
                        if abs(Q[i]) > newB:
                            #height too big, so not a valid point
                            bad = 1
                            break
                    if bad == 0:
                        P = PS.point(Q, False)
                        #check that it is actually periodic
                        newP = copy(P)
                        k = 1
                        done = False
                        while done == False and k <= n:
                              newP = self(newP)
                              if newP == P:
                                  if ([P, k] in good_points) == False:
                                      good_points.append([newP, k])
                                  done = True
                              k += 1

            return(good_points)

    def rational_periodic_points(self, **kwds):
        r"""
        Determine the set of rational periodic points for an endomorphism of projective space.

        The map must be defined over `\QQ` and be an endomorphism of projective space.
        If the map is a polynomial endomorphism of `\mathbb{P}^1`, i.e. has a totally
        ramified fixed point, then the base ring can be an absolute number field.
        This is done by passing to the Weil restriction.

        The default parameter values are typically good choices for `\mathbb{P}^1`. If you are having
        trouble getting a particular map to finish, try first computing the possible periods, then
        try various different ``lifting_prime`` values.

        ALGORITHM:
            Modulo each prime of good reduction `p` determine the set of periodic points modulo `p`.
            For each cycle modulo `p` compute the set of possible periods (`mrp^e`). Take the intersection
            of the list of possible periods modulo several primes of good reduction to get a possible list
            of minimal periods of rational periodic points. Take each point modulo `p` associated to each
            of these possible periods and try to lift it to a rational point with a combination of
            `p`-adic approximation and the LLL basis reduction algorithm.

            See [Hutz]_.

        INPUT:

        kwds:

        - ``prime_bound`` - a pair (list or tuple) of positive integers that represent the
            limits of primes to use in the reduction step. Or an integer that represents the upper bound. (optional)
            default: [1,20]

        -  ``lifting_prime`` - a prime integer. (optional) argument that specifies modulo which prime to try and perform the
            lifting. default: 23

        - ``periods`` - a list of positive integers which is the list of possible periods. (optional)

        - ``bad_primes`` - a list or tuple of integer primes, the primes of bad reduction.  (optional)

        - ``ncpus`` - number of cpus to use in parallel.  (optional)
            default: all available cpus.

        OUTPUT:

        - a list of rational points in projective space.

        Examples::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = End(P)
            sage: f = H([x^2-3/4*y^2, y^2])
            sage: sorted(f.rational_periodic_points(prime_bound=20, lifting_prime=7)) # long time
            [(-1/2 : 1), (1 : 0), (3/2 : 1)]

        ::

            sage: P.<x,y,z> = ProjectiveSpace(QQ,2)
            sage: H = End(P)
            sage: f = H([2*x^3 - 50*x*z^2 + 24*z^3, 5*y^3 - 53*y*z^2 + 24*z^3, 24*z^3])
            sage: sorted(f.rational_periodic_points(prime_bound=[1,20])) # long time
            [(-3 : -1 : 1), (-3 : 0 : 1), (-3 : 1 : 1), (-3 : 3 : 1), (-1 : -1 : 1),
            (-1 : 0 : 1), (-1 : 1 : 1), (-1 : 3 : 1), (0 : 1 : 0), (1 : -1 : 1), (1
            : 0 : 0), (1 : 0 : 1), (1 : 1 : 1), (1 : 3 : 1), (3 : -1 : 1), (3 : 0 :
            1), (3 : 1 : 1), (3 : 3 : 1), (5 : -1 : 1), (5 : 0 : 1), (5 : 1 : 1), (5
            : 3 : 1)]

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = End(P)
            sage: f = H([-5*x^2 + 4*y^2, 4*x*y])
            sage: sorted(f.rational_periodic_points()) # long time
            [(-2 : 1), (-2/3 : 1), (2/3 : 1), (1 : 0), (2 : 1)]

        ::

            sage: R.<x> = QQ[]
            sage: K.<w> = NumberField(x^2-x+1)
            sage: P.<u,v> = ProjectiveSpace(K,1)
            sage: H = End(P)
            sage: f = H([u^2 + v^2,v^2])
            sage: f.rational_periodic_points()
            [(w : 1), (1 : 0), (-w + 1 : 1)]

        ::

            sage: R.<x> = QQ[]
            sage: K.<w> = NumberField(x^2-x+1)
            sage: P.<u,v> = ProjectiveSpace(K,1)
            sage: H = End(P)
            sage: f = H([u^2+v^2,u*v])
            sage: f.rational_periodic_points()
            Traceback (most recent call last):
            ...
            NotImplementedError: rational periodic points for number fields only implemented for polynomials
        """
        if not self.is_endomorphism():
            raise NotImplementedError("must be an endomorphism of projective space")
        PS = self.domain()
        K = PS.base_ring()
        if K in _NumberFields:
            if not K.is_absolute():
                raise TypeError("base field must be an absolute field")
            d = K.absolute_degree()
            #check that we are not over QQ
            if d > 1:
                if PS.dimension_relative() != 1:
                    raise NotImplementedError("rational periodic points for number fields only implemented in dimension 1")
                w = K.absolute_generator()
                #we need to dehomogenize for the Weil restriction and will check that point at infty
                #separately. We also check here that we are working with a polynomial. If the map
                #is not a polynomial, the Weil restriction will not be a morphism and we cannot
                #apply this algorithm.
                g = self.dehomogenize(1)
                inf = PS([1,0])
                k = 1
                if isinstance(g[0], FractionFieldElement):
                    g = self.dehomogenize(0)
                    inf = PS([0,1])
                    k = 0
                    if isinstance(g[0], FractionFieldElement):
                        raise NotImplementedError("rational periodic points for number fields only implemented for polynomials")
                #determine rational periodic points
                #infinity is a totally ramified fixed point for a polynomial
                periodic_points = set([inf])
                #compute the weil restriction
                G = g.weil_restriction()
                F = G.homogenize(d)
                #find the QQ rational periodic points for the weil restriction
                Fper = F.rational_periodic_points(**kwds)
                for P in Fper:
                    #take the 'good' points in the weil restriction and find the
                    #associated number field points.
                    if P[d] == 1:
                        pt = [sum([P[i]*w**i for i in range(d)])]
                        pt.insert(k,1)
                        Q = PS(pt)
                        #for each periodic point get the entire cycle
                        if not Q in periodic_points:
                            #check periodic not preperiodic and add all points in cycle
                            orb = set([Q])
                            Q2 = self(Q)
                            while Q2 not in orb:
                                orb.add(Q2)
                                Q2 = self(Q2)
                            if Q2 == Q:
                                periodic_points = periodic_points.union(orb)
                return list(periodic_points)
            else:
                primebound = kwds.pop("prime_bound", [1, 20])
                p = kwds.pop("lifting_prime", 23)
                periods = kwds.pop("periods", None)
                badprimes = kwds.pop("bad_primes", None)
                num_cpus = kwds.pop("ncpus", ncpus())

                if (isinstance(primebound, (list, tuple)) == False):
                    try:
                        primebound = [1, ZZ(primebound)]
                    except TypeError:
                        raise TypeError("bound on primes must be an integer")
                else:
                    try:
                        primebound[0] = ZZ(primebound[0])
                        primebound[1] = ZZ(primebound[1])
                    except TypeError:
                        raise TypeError("prime bounds must be integers")

                if badprimes is None:
                    badprimes = self.primes_of_bad_reduction()
                if periods is None:
                    periods = self.possible_periods(prime_bound=primebound, bad_primes=badprimes, ncpus=num_cpus)
                PS = self.domain()
                R = PS.base_ring()
                periodic = set()
                while p in badprimes:
                    p = next_prime(p + 1)
                B = e ** self.height_difference_bound()

                f = self.change_ring(GF(p))
                all_points = f.possible_periods(True) #return the list of points and their periods.
                pos_points = []
                for i in range(len(all_points)):
                    if all_points[i][1] in periods and  (all_points[i] in pos_points) == False:  #check period, remove duplicates
                        pos_points.append(all_points[i])
                periodic_points = self.lift_to_rational_periodic(pos_points,B)
                for p,n in periodic_points:
                    for k in range(n):
                        p.normalize_coordinates()
                        periodic.add(p)
                        p = self(p)
                return list(periodic)
        else:
            raise TypeError("base field must be an absolute number field")

    def rational_preimages(self, Q, k=1):
        r"""
        Determine all of the rational `k`-th preimages of ``Q`` by this map.

        Given a rational point ``Q`` in the domain of this map, return all the rational points ``P``
        in the domain with `f^k(P)==Q`. In other words, the set of `k`-th preimages of ``Q``.
        The map must be defined over a number field and be an endomorphism for `k > 1`.

        If ``Q`` is a subscheme, then return the subscheme that maps to ``Q`` by this map.
        In particular, `f^{-k}(V(h_1,\ldots,h_t)) = V(h_1 \circ f^k, \ldots, h_t \circ f^k)`.

        INPUT:

        - ``Q`` - a rational point or subscheme in the domain of this map.

        - ``k`` - positive integer.

        OUTPUT:

        - a list of rational points or a subscheme in the domain of this map.

        Examples::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: H = End(P)
            sage: f = H([16*x^2 - 29*y^2, 16*y^2])
            sage: f.rational_preimages(P(-1, 4))
            [(-5/4 : 1), (5/4 : 1)]

        ::

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: H = End(P)
            sage: f = H([76*x^2 - 180*x*y + 45*y^2 + 14*x*z + 45*y*z\
            - 90*z^2, 67*x^2 - 180*x*y - 157*x*z + 90*y*z, -90*z^2])
            sage: f.rational_preimages(P(-9, -4, 1))
            [(0 : 4 : 1)]

        A non-periodic example ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: H = End(P)
            sage: f = H([x^2 + y^2, 2*x*y])
            sage: f.rational_preimages(P(17, 15))
            [(3/5 : 1), (5/3 : 1)]

        ::

            sage: P.<x,y,z,w> = ProjectiveSpace(QQ, 3)
            sage: H = End(P)
            sage: f = H([x^2 - 2*y*w - 3*w^2, -2*x^2 + y^2 - 2*x*z\
            + 4*y*w + 3*w^2, x^2 - y^2 + 2*x*z + z^2 - 2*y*w - w^2, w^2])
            sage: f.rational_preimages(P(0, -1, 0, 1))
            []

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: H = End(P)
            sage: f = H([x^2 + y^2, 2*x*y])
            sage: f.rational_preimages([CC.0, 1])
            Traceback (most recent call last):
            ...
            TypeError: point must be in codomain of self

        A number field example ::

            sage: z = QQ['z'].0
            sage: K.<a> = NumberField(z^2 - 2);
            sage: P.<x,y> = ProjectiveSpace(K, 1)
            sage: H = End(P)
            sage: f = H([x^2 + y^2, y^2])
            sage: f.rational_preimages(P(3, 1))
            [(-a : 1), (a : 1)]

        ::

            sage: z = QQ['z'].0
            sage: K.<a> = NumberField(z^2 - 2);
            sage: P.<x,y,z> = ProjectiveSpace(K, 2)
            sage: X = P.subscheme([x^2 - z^2])
            sage: H = End(X)
            sage: f= H([x^2 - z^2, a*y^2, z^2 - x^2])
            sage: f.rational_preimages(X([1, 2, -1]))
            []

        ::

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: X = P.subscheme([x^2 - z^2])
            sage: H = End(X)
            sage: f= H([x^2-z^2, y^2, z^2-x^2])
            sage: f.rational_preimages(X([0, 1, 0]))
            Traceback (most recent call last):
            ...
            NotImplementedError: subschemes as preimages not implemented

        ::

            sage: P.<x, y> = ProjectiveSpace(QQ, 1)
            sage: H = End(P)
            sage: f = H([x^2-y^2, y^2])
            sage: f.rational_preimages(P.subscheme([x]))
            Closed subscheme of Projective Space of dimension 1 over Rational Field
            defined by:
              x^2 - y^2

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: H = End(P)
            sage: f = H([x^2 - 29/16*y^2, y^2])
            sage: f.rational_preimages(P(5/4, 1), k=4)
            [(-3/4 : 1), (3/4 : 1), (-7/4 : 1), (7/4 : 1)]

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: P2.<u,v,w> = ProjectiveSpace(QQ, 2)
            sage: H = Hom(P, P2)
            sage: f = H([x^2, y^2, x^2-y^2])
            sage: f.rational_preimages(P2(1, 1, 0))
            [(-1 : 1), (1 : 1)]
        """
        k = ZZ(k)
        if k <= 0:
            raise ValueError("k (=%s) must be a positive integer"%(k))
        #first check if subscheme
        from sage.schemes.generic.algebraic_scheme import AlgebraicScheme_subscheme_projective
        if isinstance(Q, AlgebraicScheme_subscheme_projective):
            return(Q.preimage(self, k))

        #else assume a point
        BR = self.base_ring()
        if k > 1 and not self.is_endomorphism():
            raise TypeError("must be an endomorphism of projective space")
        if not Q in self.codomain():
            raise TypeError("point must be in codomain of self")
        if isinstance(BR.base_ring(),(ComplexField_class, RealField_class,RealIntervalField_class, ComplexIntervalField_class)):
            raise NotImplementedError("not implemented over precision fields")
        PS = self.domain().ambient_space()
        N = PS.dimension_relative()
        L = [Q]
        for n in range(k):
            L2 = []
            for P in L:
                I = list(self.domain().defining_polynomials())
                for i in range(N+1):
                    for j in range(i+1, N+1):
                        I.append(P[i]*self[j] - P[j]*self[i])
                X = PS.subscheme(I)
                if X.dimension() > 0:
                    raise NotImplementedError("subschemes as preimages not implemented")
                preimages = []
                for T in X.rational_points():
                    if not all([g(tuple(T)) == 0 for g in self]):
                        preimages.append(PS(T))
                L2 = L2 + preimages
            L = L2
        return L

    def all_rational_preimages(self, points):
        r"""
        Given a set of rational points in the domain of this map, return all the rational
        preimages of those points.

        In others words, all the rational points which have some
        iterate in the set points. This function repeatedly calls ``rational_preimages``.
        If the degree is at least two, by Northocott, this is always a finite set.
        The map must be defined over number fields and be an endomorphism.

        INPUT:

        - ``points`` - a list of rational points in the domain of this map.

        OUTPUT:

        - a list of rational points in the domain of this map.

        Examples::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = End(P)
            sage: f = H([16*x^2 - 29*y^2, 16*y^2])
            sage: sorted(f.all_rational_preimages([P(-1,4)]))
            [(-7/4 : 1), (-5/4 : 1), (-3/4 : 1), (-1/4 : 1), (1/4 : 1), (3/4 : 1),
            (5/4 : 1), (7/4 : 1)]

        ::

            sage: P.<x,y,z> = ProjectiveSpace(QQ,2)
            sage: H = End(P)
            sage: f = H([76*x^2 - 180*x*y + 45*y^2 + 14*x*z + 45*y*z - 90*z^2, 67*x^2 - 180*x*y - 157*x*z + 90*y*z, -90*z^2])
            sage: sorted(f.all_rational_preimages([P(-9,-4,1)]))
            [(-9 : -4 : 1), (0 : -1 : 1), (0 : 0 : 1), (0 : 1 : 1), (0 : 4 : 1), (1
            : 0 : 1), (1 : 1 : 1), (1 : 2 : 1), (1 : 3 : 1)]

        A non-periodic example ::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = End(P)
            sage: f = H([x^2 + y^2, 2*x*y])
            sage: sorted(f.all_rational_preimages([P(17,15)]))
            [(1/3 : 1), (3/5 : 1), (5/3 : 1), (3 : 1)]

        A number field example.::

            sage: z = QQ['z'].0
            sage: K.<w> = NumberField(z^3 + (z^2)/4 - (41/16)*z + 23/64);
            sage: P.<x,y> = ProjectiveSpace(K,1)
            sage: H = End(P)
            sage: f = H([16*x^2 - 29*y^2, 16*y^2])
            sage: f.all_rational_preimages([P(16*w^2 - 29,16)])
            [(w^2 + w - 25/16 : 1), (-w - 1/2 : 1), (w^2 - 29/16 : 1), (-w : 1), (w + 1/2 : 1), (w^2 - 21/16 : 1),
            (w : 1), (-w^2 + 21/16 : 1), (-w^2 - w + 25/16 : 1), (w^2 + w - 33/16 : 1), (-w^2 + 29/16 : 1),
            (-w^2 - w + 33/16 : 1)]

        ::

            sage: K.<w> = QuadraticField(3)
            sage: P.<u,v> = ProjectiveSpace(K,1)
            sage: H = End(P)
            sage: f = H([u^2+v^2, v^2])
            sage: f.all_rational_preimages(P(4))
            [(-w : 1), (w : 1)]
        """
        if not self.is_endomorphism():
            raise NotImplementedError("must be an endomorphism of projective space")
        if self.domain().base_ring() not in NumberFields():
            raise TypeError("field won't return finite list of elements")
        if not isinstance(points, (list, tuple)):
            points = [points]

        PS = self.domain()
        RPS = PS.base_ring()
        preperiodic = set()
        while points != []:
            P = points.pop()
            preimages = self.rational_preimages(P)
            for i in range(len(preimages)):
                if not preimages[i] in preperiodic:
                    points.append(preimages[i])
                    preperiodic.add(preimages[i])
        return(list(preperiodic))

    def rational_preperiodic_points(self, **kwds):
        r"""
        Determine the set of rational preperiodic points for this map.

        The map must be defined over `\QQ` and be an endomorphism of projective space.
        If the map is a polynomial endomorphism of `\mathbb{P}^1`, i.e. has a totally
        ramified fixed point, then the base ring can be an absolute number field.
        This is done by passing to the Weil restriction.

        The default parameter values are typically good choices for `\mathbb{P}^1`. If you are having
        trouble getting a particular map to finish, try first computing the possible periods, then
        try various different values for ``lifting_prime``.

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

        - ``ncpus`` - number of cpus to use in parallel.  (optional)
            default: all available cpus.

        OUTPUT:

        - a list of rational points in projective space.

        Examples::

            sage: PS.<x,y> = ProjectiveSpace(1,QQ)
            sage: H = End(PS)
            sage: f = H([x^2 -y^2, 3*x*y])
            sage: sorted(f.rational_preperiodic_points())
            [(-2 : 1), (-1 : 1), (-1/2 : 1), (0 : 1), (1/2 : 1), (1 : 0), (1 : 1),
            (2 : 1)]

        ::

            sage: PS.<x,y> = ProjectiveSpace(1,QQ)
            sage: H = End(PS)
            sage: f = H([5*x^3 - 53*x*y^2 + 24*y^3, 24*y^3])
            sage: sorted(f.rational_preperiodic_points(prime_bound=10))
            [(-1 : 1), (0 : 1), (1 : 0), (1 : 1), (3 : 1)]

        ::

            sage: PS.<x,y,z> = ProjectiveSpace(2,QQ)
            sage: H = End(PS)
            sage: f = H([x^2 - 21/16*z^2, y^2-2*z^2, z^2])
            sage: sorted(f.rational_preperiodic_points(prime_bound=[1,8], lifting_prime=7, periods=[2])) # long time
            [(-5/4 : -2 : 1), (-5/4 : -1 : 1), (-5/4 : 0 : 1), (-5/4 : 1 : 1), (-5/4
            : 2 : 1), (-1/4 : -2 : 1), (-1/4 : -1 : 1), (-1/4 : 0 : 1), (-1/4 : 1 :
            1), (-1/4 : 2 : 1), (1/4 : -2 : 1), (1/4 : -1 : 1), (1/4 : 0 : 1), (1/4
            : 1 : 1), (1/4 : 2 : 1), (5/4 : -2 : 1), (5/4 : -1 : 1), (5/4 : 0 : 1),
            (5/4 : 1 : 1), (5/4 : 2 : 1)]

        ::

            sage: K.<w> = QuadraticField(33)
            sage: PS.<x,y> = ProjectiveSpace(K,1)
            sage: H = End(PS)
            sage: f = H([x^2-71/48*y^2, y^2])
            sage: sorted(f.rational_preperiodic_points()) # long time
            [(-1/12*w - 1 : 1),
             (-1/6*w - 1/4 : 1),
             (-1/12*w - 1/2 : 1),
             (-1/6*w + 1/4 : 1),
             (1/12*w - 1 : 1),
             (1/12*w - 1/2 : 1),
             (-1/12*w + 1/2 : 1),
             (-1/12*w + 1 : 1),
             (1/6*w - 1/4 : 1),
             (1/12*w + 1/2 : 1),
             (1 : 0),
             (1/6*w + 1/4 : 1),
             (1/12*w + 1 : 1)]
        """
        PS = self.domain()
        K = PS.base_ring()
        if K in _NumberFields:
            if not K.is_absolute():
                raise TypeError("base field must be an absolute field")
            d = K.absolute_degree()
            #check that we are not over QQ
            if d > 1:
                if PS.dimension_relative() != 1:
                    raise NotImplementedError("rational preperiodic points for number fields only implemented in dimension 1")
                w = K.absolute_generator()
                #we need to dehomogenize for the Weil restriction and will check that point at infty
                #separately. We also check here that we are working with a polynomial. If the map
                #is not a polynomial, the Weil restriction will not be a morphism and we cannot
                #apply this algorithm.
                g = self.dehomogenize(1)
                inf = PS([1,0])
                k = 1
                if isinstance(g[0], FractionFieldElement):
                    g = self.dehomogenize(0)
                    inf = PS([0,1])
                    k = 0
                    if isinstance(g[0], FractionFieldElement):
                        raise NotImplementedError("rational preperiodic points for number fields only implemented for polynomials")
                #determine rational preperiodic points
                #infinity is a totally ramified fixed point for a polynomial
                preper = set([inf])
                #compute the weil restriction
                G = g.weil_restriction()
                F = G.homogenize(d)
                #find the QQ rational preperiodic points for the weil restriction
                Fpre = F.rational_preperiodic_points(**kwds)
                for P in Fpre:
                    #take the 'good' points in the weil restriction and find the
                    #associated number field points.
                    if P[d] == 1:
                        pt = [sum([P[i]*w**i for i in range(d)])]
                        pt.insert(k,1)
                        Q = PS(pt)
                        #for each preperiodic point get the entire connected component
                        if not Q in preper:
                            for t in self.connected_rational_component(Q):
                                preper.add(t)
                preper = list(preper)
            else:
                #input error checking done in possible_periods and rational_periodic_points
                badprimes = kwds.pop("bad_primes", None)
                periods = kwds.pop("periods", None)
                primebound = kwds.pop("prime_bound", [1, 20])
                num_cpus = kwds.pop("ncpus", ncpus())
                if badprimes is None:
                    badprimes = self.primes_of_bad_reduction()
                if periods is None:
                    periods = self.possible_periods(prime_bound=primebound, bad_primes=badprimes, ncpus=num_cpus) #determine the set of possible periods
                if periods == []:
                    return([]) #no rational preperiodic points
                else:
                    p = kwds.pop("lifting_prime", 23)
                    T = self.rational_periodic_points(prime_bound=primebound, lifting_prime=p, periods=periods, bad_primes=badprimes, ncpus=num_cpus) #find the rational preperiodic points
                    preper = self.all_rational_preimages(T) #find the preperiodic points
                    preper = list(preper)
            return(preper)
        else:
            raise TypeError("base field must be an absolute number field")

    def rational_preperiodic_graph(self, **kwds):
        r"""
        Determine the directed graph of the rational preperiodic points for this map.

        The map must be defined over `\QQ` and be an endomorphism of projective space.
        If this map is a polynomial endomorphism of `\mathbb{P}^1`, i.e. has a totally
        ramified fixed point, then the base ring can be an absolute number field.
        This is done by passing to the Weil restriction.

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

        - ``ncpus`` - number of cpus to use in parallel.  (optional)
            default: all available cpus.

        OUTPUT:

        - a digraph representing the orbits of the rational preperiodic points in projective space.

        Examples::

            sage: PS.<x,y> = ProjectiveSpace(1,QQ)
            sage: H = End(PS)
            sage: f = H([7*x^2 - 28*y^2, 24*x*y])
            sage: f.rational_preperiodic_graph()
            Looped digraph on 12 vertices

        ::

            sage: PS.<x,y> = ProjectiveSpace(1,QQ)
            sage: H = End(PS)
            sage: f = H([-3/2*x^3 +19/6*x*y^2, y^3])
            sage: f.rational_preperiodic_graph(prime_bound=[1,8])
            Looped digraph on 12 vertices

        ::

            sage: PS.<x,y,z> = ProjectiveSpace(2,QQ)
            sage: H = End(PS)
            sage: f = H([2*x^3 - 50*x*z^2 + 24*z^3, 5*y^3 - 53*y*z^2 + 24*z^3, 24*z^3])
            sage: f.rational_preperiodic_graph(prime_bound=[1,11], lifting_prime=13) # long time
            Looped digraph on 30 vertices

        ::

            sage: K.<w> = QuadraticField(-3)
            sage: P.<x,y> = ProjectiveSpace(K,1)
            sage: H = End(P)
            sage: f = H([x^2+y^2, y^2])
            sage: f.rational_preperiodic_graph() # long time
            Looped digraph on 5 vertices
        """
        #input checking done in .rational_preperiodic_points()
        preper = self.rational_preperiodic_points(**kwds)
        g = self._preperiodic_points_to_cyclegraph(preper)
        return(g)

    def connected_rational_component(self, P, n=0):
        r"""
        Computes the connected component of a rational preperiodic point ``P`` by this map.

        Will work for non-preperiodic points if ``n`` is positive.
        Otherwise this will not terminate.

        INPUT:

        - ``P`` - A rational preperiodic point of this map.

        - ``n`` - Maximum distance from ``P`` to branch out. A value of 0 indicates no bound. Default: 0

        OUTPUT:

        - a list of points connected to ``P`` up to the specified distance.

        Examples::

            sage: R.<x> = PolynomialRing(QQ)
            sage: K.<w> = NumberField(x^3+1/4*x^2-41/16*x+23/64)
            sage: PS.<x,y> = ProjectiveSpace(1,K)
            sage: H = End(PS)
            sage: f = H([x^2 - 29/16*y^2, y^2])
            sage: P = PS([w,1])
            sage: f.connected_rational_component(P)
            [(w : 1),
             (w^2 - 29/16 : 1),
             (w^2 + w - 25/16 : 1),
             (-w^2 - w + 25/16 : 1),
             (-w : 1),
             (w + 1/2 : 1),
             (-w - 1/2 : 1),
             (-w^2 + 29/16 : 1),
             (-w^2 + 21/16 : 1),
             (w^2 - 21/16 : 1),
             (w^2 + w - 33/16 : 1),
             (-w^2 - w + 33/16 : 1)]

        ::

            sage: PS.<x,y,z> = ProjectiveSpace(2,QQ)
            sage: H = End(PS)
            sage: f = H([x^2 - 21/16*z^2, y^2-2*z^2, z^2])
            sage: P = PS([17/16,7/4,1])
            sage: f.connected_rational_component(P,3)
            [(17/16 : 7/4 : 1),
             (-47/256 : 17/16 : 1),
             (-83807/65536 : -223/256 : 1),
             (-17/16 : -7/4 : 1),
             (-17/16 : 7/4 : 1),
             (17/16 : -7/4 : 1),
             (1386468673/4294967296 : -81343/65536 : 1),
             (-47/256 : -17/16 : 1),
             (47/256 : -17/16 : 1),
             (47/256 : 17/16 : 1),
             (-1/2 : -1/2 : 1),
             (-1/2 : 1/2 : 1),
             (1/2 : -1/2 : 1),
             (1/2 : 1/2 : 1)]

        """
        points = [[],[]] # list of points and a list of their corresponding levels
        points[0].append(P)
        points[1].append(0) # P is treated as level 0

        nextpoints = []
        nextpoints.append(P)

        level = 1
        foundall = False # whether done or not
        while not foundall:
            newpoints = []
            for Q in nextpoints:
                # forward image
                newpoints.append(self(Q))
                # preimages
                newpoints.extend(self.rational_preimages(Q))
            del nextpoints[:] # empty list
            # add any points that are not already in the connected component
            for Q in newpoints:
                if (Q not in points[0]):
                    points[0].append(Q)
                    points[1].append(level)
                    nextpoints.append(Q)
            # done if max level was achieved or if there were no more points to add
            if ((level + 1 > n and n != 0) or len(nextpoints) == 0):
                foundall = True
            level = level + 1

        return points[0]

    def _number_field_from_algebraics(self):
        r"""
        Given a projective map defined over `\QQbar`, return the same map, but defined
        over a number field.

        This is only implemented for maps of projective space.

        OUTPUT: scheme morphism

        EXAMPLES::

            sage: R.<x> = PolynomialRing(QQ)
            sage: P.<x,y> = ProjectiveSpace(QQbar,1)
            sage: H = End(P)
            sage: f = H([QQbar(3^(1/3))*x^2 + QQbar(sqrt(-2))*y^2, y^2])
            sage: f._number_field_from_algebraics()
            Scheme endomorphism of Projective Space of dimension 1 over Number Field
            in a with defining polynomial y^6 + 6*y^4 + 6*y^3 + 12*y^2 - 36*y + 17
              Defn: Defined on coordinates by sending (z0 : z1) to
                    ((48/269*a^5 + 27/269*a^4 + 320/269*a^3 + 468/269*a^2 + 772/269*a
                    - 1092/269)*z0^2 + (48/269*a^5 + 27/269*a^4 + 320/269*a^3 + 468/269*a^2
                    + 1041/269*a - 1092/269)*z1^2 : z1^2)

        ::

            sage: R.<x> = PolynomialRing(QQ)
            sage: K.<a> = NumberField(x^3-x+1)
            sage: P.<x,y> = ProjectiveSpace(QQbar,1)
            sage: P2.<u,v,w> = ProjectiveSpace(QQbar,2)
            sage: H = Hom(P, P2)
            sage: f = H([x^2 + QQbar(I)*x*y + 3*y^2, y^2, QQbar(sqrt(5))*x*y])
            sage: f._number_field_from_algebraics()
            Scheme morphism:
              From: Projective Space of dimension 1 over Number Field in a with defining polynomial y^4 + 3*y^2 + 1
              To:   Projective Space of dimension 2 over Number Field in a with defining polynomial y^4 + 3*y^2 + 1
              Defn: Defined on coordinates by sending (z0 : z1) to
                    (z0^2 + (a^3 + 2*a)*z0*z1 + 3*z1^2 : z1^2 : (2*a^2 + 3)*z0*z1)
        """
        from sage.schemes.projective.projective_space import is_ProjectiveSpace
        if not (is_ProjectiveSpace(self.domain()) and is_ProjectiveSpace(self.domain())):
            raise NotImplementedError("not implemented for subschemes")

        K,C,phi = number_field_elements_from_algebraics([c for f in self for c in f.coefficients()])
        from sage.schemes.projective.projective_space import ProjectiveSpace
        N = self.domain().dimension_relative()
        PS = ProjectiveSpace(K,N,'z')
        if self.is_endomorphism():
            H = End(PS)
        else:
            PS2 = ProjectiveSpace(K,self.codomain().dimension_relative(),'w')
            H = Hom(PS,PS2)
        R = PS.coordinate_ring()
        exps = [f.exponents() for f in self]
        F = []
        j = 0
        for t in exps:
            G = 0
            for e in t:
                G += C[j]*prod([R.gen(i)**e[i] for i in range(N+1)])
                j += 1
            F.append(G)
        return(H(F))

    def conjugating_set(self, other):
        r"""
        Returns the set of elements in PGL that conjugates one map to the other.

        Given two nonconstant rational functions of equal degree determine to see if there is an element of PGL that
        conjugates one rational function to another. It does this by taking the fixed points of one map and mapping
        them to all unique permutations of the fixed points of the other map. If there are not enough fixed points the
        function compares the mapping between rational preimages of fixed points and the rational preimages of the preimages of
        fixed points until there are enough points; such that there are `n+2` points with all `n+1` subsets linearly independent.

        ALGORITHM:
        
        Implementing invariant set algorithim from the paper [FMV]_. Given that the set of  `n` th preimages of fixed points is
        invariant under conjugation find all elements of PGL that take one set to another.

        INPUT: Two nonconstant rational functions of same degree.

        OUTPUT: Set of conjugating `n+1` by `n+1` matrices.

        AUTHORS:
        
        - Original algorithm written by Xander Faber, Michelle Manes, Bianca Viray [FMV]_.
        
        - Implimented by Rebecca Lauren Miller, as part of GSOC 2016.

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = End(P)
            sage: f = H([x^2 - 2*y^2, y^2])
            sage: m = matrix(QQbar, 2, 2, [-1, 3, 2, 1])
            sage: g = f.conjugate(m)
            sage: f.conjugating_set(g)
            [
            [-1  3]
            [ 2  1]
            ]

        ::

            sage: K.<w> = QuadraticField(-1)
            sage: P.<x,y> = ProjectiveSpace(K,1)
            sage: H = End(P)
            sage: f = H([x^2 + y^2, x*y])
            sage: m = matrix(K, 2, 2, [1, 1, 2, 1])
            sage: g = f.conjugate(m)
            sage: f.conjugating_set(g) # long test
            [
            [1 1]  [-1 -1]
            [2 1], [ 2  1]
            ]

        ::

            sage: K.<i> = QuadraticField(-1)
            sage: P.<x,y> = ProjectiveSpace(K,1)
            sage: H = End(P)
            sage: D8 = H([y^3, x^3])
            sage: D8.conjugating_set(D8) # long test
            [
            [1 0]  [0 1]  [ 0 -i]  [i 0]  [ 0 -1]  [-1  0]  [-i  0]  [0 i]
            [0 1], [1 0], [ 1  0], [0 1], [ 1  0], [ 0  1], [ 0  1], [1 0]
            ]

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = End(P)
            sage: D8 = H([y^2, x^2])
            sage: D8.conjugating_set(D8)
            Traceback (most recent call last):
            ...
            ValueError: not enough rational preimages

        ::
        
            sage: P.<x,y> = ProjectiveSpace(GF(7),1)
            sage: H = End(P)
            sage: D6 = H([y^2, x^2])
            sage: D6.conjugating_set(D6)
            [
            [1 0]  [0 1]  [0 2]  [4 0]  [2 0]  [0 4]
            [0 1], [1 0], [1 0], [0 1], [0 1], [1 0]
            ]

        ::

            sage: P.<x,y,z> = ProjectiveSpace(QQ,2)
            sage: H = End(P)
            sage: f = H([x^2 + x*z, y^2, z^2])
            sage: f.conjugating_set(f) # long test
            [
            [1 0 0]
            [0 1 0]
            [0 0 1]
            ]
        """
        f = copy(self)
        g = copy(other)
        try:
            f.normalize_coordinates()
            g.normalize_coordinates()
        except (ValueError):
            pass
        if f.degree() != g.degree():# checks that maps are of equal degree
            return []
        n = f.domain().dimension_relative()
        L = Set(f.periodic_points(1))
        K = Set(g.periodic_points(1))
        if len(L) != len(K):  # checks maps have the same number of fixed points
            return []
        d = len(L)
        r = f.domain().base_ring()
        more = True
        if d >= n+2: # need at least n+2 points
            for i in Subsets(L, n+2):# makes sure all n+1 subsets are linearly independent
                Ml = matrix(r, [list(s) for s in i])
                if not any([j == 0 for j in Ml.minors(n+1)]):
                    Tf = list(i)
                    more= False
                    break
        while more:
            Tl = [Q for i in L for Q in f.rational_preimages(i)] #  finds preimages of fixed points
            Tk = [Q for i in K for Q in g.rational_preimages(i)]
            if len(Tl) != len(Tk):
                return []
            L = L.union(Set(Tl))
            K = K.union(Set(Tk))
            if d == len(L): # if no new preimages then not enough points
                raise ValueError("not enough rational preimages")
            d = len(L)
            if d >= n+2: # makes sure all n+1 subsets are linearly independent
                for i in Subsets(L, n+2):
                    Ml = matrix(r, [list(s) for s in i])
                    if not any([j == 0 for j in Ml.minors(n+1)]):
                        more = False
                        Tf = list(i)
                        break
        Conj = []
        for i in Arrangements(K,(n+2)): # try all possible conjugations between invariant sets
            try: # need all n+1 subsets linearly independenet
                s = f.domain().point_transformation_matrix(i,Tf)# finds elements of PGL that maps one map to another
                if self.conjugate(s) == other:
                    Conj.append(s)
            except (ValueError):
                pass
        return Conj

    def is_conjugate(self, other):
        r"""
        Returns whether or not two maps are conjugate.

        ALGORITHM:

        Implementing invariant set algorithim from the paper [FMV]_. Given that the set of `n` th preimages is
        invariant under conjugation this function finds whether two maps are conjugate.

        INPUT: Two nonconstant rational functions of same degree.

        OUTPUT: Boolean.

        AUTHORS:

        - Original algorithm written by Xander Faber, Michelle Manes, Bianca Viray [FMV]_.

        - Implimented by Rebecca Lauren Miller as part of GSOC 2016.

        EXAMPLES::

            sage: K.<w> = CyclotomicField(3)
            sage: P.<x,y> = ProjectiveSpace(K,1)
            sage: H = End(P)
            sage: D8 = H([y^2, x^2])
            sage: D8.is_conjugate(D8)
            True

        ::

            sage: set_verbose(None)
            sage: P.<x,y> = ProjectiveSpace(QQbar,1)
            sage: H = End(P)
            sage: f = H([x^2 + x*y,y^2])
            sage: m = matrix(QQbar, 2, 2, [1, 1, 2, 1])
            sage: g = f.conjugate(m)
            sage: f.is_conjugate(g) # long test
            True
            
        ::
        
            sage: P.<x,y> = ProjectiveSpace(GF(5),1)
            sage: H = End(P)
            sage: f = H([x^3 + x*y^2,y^3])
            sage: m = matrix(GF(5), 2, 2, [1, 3, 2, 9])
            sage: g = f.conjugate(m)
            sage: f.is_conjugate(g)
            True
            
        ::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = End(P)
            sage: f = H([x^2 + x*y,y^2])
            sage: g = H([x^3 + x^2*y, y^3])
            sage: f.is_conjugate(g) 
            False

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = End(P)
            sage: f = H([x^2 + x*y, y^2])
            sage: g = H([x^2 - 2*y^2, y^2])
            sage: f.is_conjugate(g)
            False
        """
        f = copy(self)
        g = copy(other)
        try:
            f.normalize_coordinates()
            g.normalize_coordinates()
        except (ValueError):
            pass
        if f.degree() != g.degree(): # checks that maps are of equal degree
            return False
        n = f.domain().dimension_relative()
        L = Set(f.periodic_points(1))
        K = Set(g.periodic_points(1))
        if len(L) != len(K): # checks maps have the same number of fixed points
            return False
        d = len(L)
        r = f.domain().base_ring()
        more = True
        if d >= n+2: # need at least n+2 points
            for i in Subsets(L, n+2): # makes sure all n+1 subsets are linearly independent
                Ml = matrix(r, [list(s) for s in i])
                if not any([j == 0 for j in Ml.minors(n+1)]):
                    Tf = list(i)
                    more = False
                    break
        while more:
            Tl = [Q for i in L for Q in f.rational_preimages(i)] # finds preimages of fixed points
            Tk = [Q for i in K for Q in g.rational_preimages(i)]
            if len(Tl) != len(Tk):
                return False
            L = L.union(Set(Tl))
            K = K.union(Set(Tk))
            if d == len(L):# if no new preimages then not enough points
                raise ValueError("not enough rational preimages")
            d = len(L)
            if d >= n+2: # makes sure all n+1 subsets are linearly independent
                for i in Subsets(L, n+2): # checks at least n+1 are linearly independent
                    Ml = matrix(r, [list(s) for s in i])
                    if not any([j == 0 for j in Ml.minors(n+1)]):
                        more = False
                        Tf = list(i)
                        break
        Conj = []
        for i in Arrangements(K,n+2):# try all possible conjugations between invariant sets
            try: # need all n+1 subsets linearly independenet
                s = f.domain().point_transformation_matrix(i,Tf) # finds elements of PGL that maps one map to another
                if self.conjugate(s) == other:
                    return True
            except (ValueError):
                pass
        return False

    def is_polynomial(self):
        r"""
        Checks to see if the function has a totally ramified fixed point.

        The function must be defined over an absolute number field or a
        finite field.

        OUTPUT: Boolean

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: K.<w> = QuadraticField(7)
            sage: P.<x,y> = ProjectiveSpace(K, 1)
            sage: H = End(P)
            sage: f = H([x**2 + 2*x*y - 5*y**2, 2*x*y])
            sage: f.is_polynomial()
            False

        ::

            sage: R.<x> = QQ[]
            sage: K.<w> = QuadraticField(7)
            sage: P.<x,y> = ProjectiveSpace(K, 1)
            sage: H = End(P)
            sage: f = H([x**2 - 7*x*y, 2*y**2])
            sage: m = matrix(K, 2, 2, [w, 1, 0, 1])
            sage: f = f.conjugate(m)
            sage: f.is_polynomial()
            True

        ::

            sage: K.<w> = QuadraticField(4/27)
            sage: P.<x,y> = ProjectiveSpace(K,1)
            sage: H = End(P)
            sage: S = P.coordinate_ring()
            sage: f = H([x**3 + w*y^3,x*y**2])
            sage: f.is_polynomial()
            False

        ::

            sage: K = GF(3**2, prefix='w')
            sage: P.<x,y> = ProjectiveSpace(K,1)
            sage: H = End(P)
            sage: f = H([x**2 + K.gen()*y**2, x*y])
            sage: f.is_polynomial()
            False

        ::

            sage: PS.<x,y> = ProjectiveSpace(QQ, 1)
            sage: H = End(PS)
            sage: f = H([6*x^2+12*x*y+7*y^2, 12*x*y + 42*y^2])
            sage: f.is_polynomial()
            False
        """
        if self.codomain().dimension_relative() != 1:
            raise NotImplementedError("space must have dimension equal to 1")
        K = self.base_ring()
        if not K in FiniteFields() and (not K in _NumberFields or not K.is_absolute()):
            raise NotImplementedError("must be over an absolute number field or finite field")
        if K in FiniteFields():
            q = K.characteristic()
            deg = K.degree()
            var = K.variable_name()
        g = self
        #get polynomial defining fixed points
        G = self.dehomogenize(1).dynatomic_polynomial(1)
        # see if infty = (1,0) is fixed
        if G.degree() <= g.degree():
            #check if infty is totally ramified
            if len((g[1]).factor()) == 1:
                return True
        #otherwise we need to create the tower of extensions
        #which contain the fixed points. We do
        #this successively so we can exit early if
        #we find one and not go all the way to the splitting field
        i = 0 #field index
        if G.degree() != 0:
            G = G.polynomial(G.variable(0))
        while G.degree() != 0:
            Y = G.factor()
            R = G.parent()
            u = G
            for p,e in Y:
                if p.degree() == 1:
                    if len((g[0]*p[1] + g[1]*p[0]).factor()) == 1:
                        return True
                    G = R(G/p) # we already checked this root
                else:
                    u = p #need to extend to get these roots
            if G.degree() != 0:
                #create the next extension
                if K == QQ:
                    L = NumberField(u, 't'+str(i))
                    i += 1
                    phi = K.embeddings(L)[0]
                    K = L
                elif K in FiniteFields():
                    deg = deg*G.degree()
                    K = GF(q**(deg), prefix=var)
                else:
                    L = K.extension(u, 't'+str(i))
                    i += 1
                    phi1 = K.embeddings(L)[0]
                    K = L
                    L = K.absolute_field('t'+str(i))
                    i += 1
                    phi = K.embeddings(L)[0]*phi1
                    K = L
                if K in FiniteFields():
                    G = G.change_ring(K)
                    g = g.change_ring(K)
                else:
                    G = G.change_ring(phi)
                    g = g.change_ring(phi)
        return False

    def normal_form(self, return_conjugation=False):
        r"""
        Returns a normal form for the map in the moduli space of dynamical systems.

        Currently implemented only for polynomials. The totally ramified fixed point is
        moved to infinity and the map is conjugated to the form
        `x^n + a_{n-2}x^{n-2} + \cdots + a_{0}`. Note that for finite fields
        we can only remove the `(n-1)`-st term when the characteristic
        does not divide `n`.

        INPUT:

        - ``return_conjugation`` -- Boolean - True returns conjugatation element of PGL.
          along with the embedding into the new field. Default: False. (optional)

        OUTPUT:

        - :class:`SchemeMorphism_polynomial`

        - Element of PGL as a matrix. (optional)

        - Field embedding. (option)

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: H = End(P)
            sage: f = H([x^2 + 2*x*y - 5*x^2, 2*x*y])
            sage: f.normal_form()
            Traceback (most recent call last):
            ...
            NotImplementedError: map is not a polynomial

        ::

            sage: R.<x> = QQ[]
            sage: K.<w> = NumberField(x^2 - 5)
            sage: P.<x,y> = ProjectiveSpace(K,1)
            sage: H = End(P)
            sage: f = H([x^2 + w*x*y, y^2])
            sage: g,m,psi = f.normal_form(return_conjugation = True);m
            [     1 -1/2*w]
            [     0      1]
            sage: f.change_ring(psi).conjugate(m) == g
            True

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: H = End(P)
            sage: f = H([13*x^2 + 4*x*y + 3*y^2, 5*y^2])
            sage: f.normal_form()
            Scheme endomorphism of Projective Space of dimension 1 over Rational Field
              Defn: Defined on coordinates by sending (x : y) to
                    (5*x^2 + 9*y^2 : 5*y^2)

        ::

            sage: K = GF(3^3, prefix = 'w')
            sage: P.<x,y> = ProjectiveSpace(K,1)
            sage: H = End(P)
            sage: f = H([x^3 + 2*x^2*y + 2*x*y^2 + K.gen()*y^3, y^3])
            sage: f.normal_form()
            Scheme endomorphism of Projective Space of dimension 1 over Finite Field
            in w3 of size 3^3
              Defn: Defined on coordinates by sending (x : y) to
                    (x^3 + x^2*y + x*y^2 + (-w3)*y^3 : y^3)
        """
        #defines the field of fixed points
        if self.codomain().dimension_relative() != 1:
            raise NotImplementedError("space must have dimension equal to 1")
        K = self.base_ring()
        if not K in FiniteFields() and (not K in _NumberFields or not K.is_absolute()):
            raise NotImplementedError("must be over an absolute number field or finite field")
        if K in FiniteFields():
            q = K.characteristic()
            deg = K.degree()
            var = K.variable_name()
        else:
            psi = K.hom([K.gen()]) #identity hom for return_embedding
        g = self
        G = self.dehomogenize(1).dynatomic_polynomial(1)
        done = False
        bad = True
        #check infty = (1,0) is fixed
        if G.degree() <= g.degree():
            #check infty totally ramified
            if len((g[1]).factor()) == 1:
                T = self.domain()(1,0)
                bad = False
                done = True
                m = matrix(K, 2, 2, [1,0,0,1])
        #otherwise we need to create the tower of extensions
        #which contain the fixed points. We do
        #this successively so we can early exit if
        #we find one and not go all the way to the splitting field
        i = 0
        if G.degree() != 0:
            G = G.polynomial(G.variable(0))
        else:
            #no other fixed points
            raise NotImplementedError("map is not a polynomial")
        #check other fixed points
        while not done:
            Y = G.factor()
            R = G.parent()
            done = True
            for p,e in Y:
                if p.degree() == 1:
                    if len((g[0]*p[1] + g[1]*p[0]).factor()) == 1:
                        T = self.domain()(-p[0], p[1])
                        bad = False
                        done = True
                        break # bc only 1 totally ramified fixed pt
                    G = R(G/p)
                else:
                    done = False
                    u = p
            if not done:
                #extend
                if K == QQ:
                    L = NumberField(u, 't'+str(i))
                    i += 1
                    phi = K.embeddings(L)[0]
                    psi = phi*psi
                    K = L
                elif K in FiniteFields():
                    deg = deg*G.degree()
                    K = GF(q**(deg), prefix=var)
                else:
                    L = K.extension(u, 't'+str(i))
                    i += 1
                    phi1 = K.embeddings(L)[0]
                    K = L
                    L = K.absolute_field('t'+str(i))
                    i += 1
                    phi = K.embeddings(L)[0]*phi1
                    psi = phi*psi
                    K = L
                #switch to the new field
                if K in FiniteFields():
                    G = G.change_ring(K)
                    g = g.change_ring(K)
                else:
                    G = G.change_ring(phi)
                    g = g.change_ring(phi)
        if bad:
            raise NotImplementedError("map is not a polynomial")
        #conjugate to normal form
        Q = T.codomain()
        #moved totally ramified fixed point to infty
        target = [T, Q(T[0]+1, 1), Q(T[0]+2, 1)]
        source = [Q(1, 0), Q(0, 1), Q(1, 1)]
        m = Q.point_transformation_matrix(source, target)
        N = g.base_ring()
        d = g.degree()
        gc = g.conjugate(m)
        #make monic
        R = PolynomialRing(N, 'z')
        v = N(gc[1].coefficient([0,d])/gc[0].coefficient([d,0]))
        #need a (d-1)-st root to make monic
        u = R.gen(0)**(d-1) - v
        if d != 2 and u.is_irreducible():
            #we need to extend again
            if N in FiniteFields():
                deg = deg*(d-1)
                M = GF(q**(deg), prefix=var)
            else:
                L = N.extension(u,'t'+str(i))
                i += 1
                phi1 = N.embeddings(L)[0]
                M = L.absolute_field('t'+str(i))
                phi = L.embeddings(M)[0]*phi1
                psi = phi*psi
            if M in FiniteFields():
                gc = gc.change_ring(M)
            else:
                gc = gc.change_ring(phi)
            m = matrix(M, 2, 2, [phi(s) for t in list(m) for s in t])
            rv = phi(v).nth_root(d-1)
        else: #root is already in the field
            M = N
            rv = v.nth_root(d-1)
        mc = matrix(M, 2, 2, [rv,0,0,1])
        gcc = gc.conjugate(mc)
        if not (M in FiniteFields() and q.divides(d)):
            #remove 2nd order term
            mc2 = matrix(M, 2, 2, [1, M((-gcc[0].coefficient([d-1, 1]) \
                /(d*gcc[1].coefficient([0, d]))).constant_coefficient()), 0, 1])
        else:
            mc2 = mc.parent().one()
        gccc = gcc.conjugate(mc2)
        if return_conjugation:
            if M in FiniteFields():
                return gccc, m*mc*mc2
            else:
                return gccc, m*mc*mc2, psi
        return gccc

    def indeterminacy_locus(self):
        r"""
        Return the indeterminacy locus of this map.

        Only for rational maps on projective space defined over a field. 
        The indeterminacy locus is the set of points in projective space at which all of the defining polynomials of the rational map simultaneously vanish.

        OUTPUT:

        - subscheme of the domain of the map.  The empty subscheme is returned as the vanishing 
          of the coordinate functions of the domain.

        EXAMPLES::

            sage: P.<x,y,z> = ProjectiveSpace(QQ,2)
            sage: H = End(P)
            sage: f = H([x*z-y*z, x^2-y^2, z^2])
            sage: f.indeterminacy_locus()
            Closed subscheme of Projective Space of dimension 2 over Rational Field defined by:
                x*z - y*z,
                x^2 - y^2,
                z^2

        ::

            sage: P.<x,y,z> = ProjectiveSpace(QQ,2)
            sage: H = End(P)
            sage: f = H([x^2, y^2, z^2])
            sage: f.indeterminacy_locus()
            Closed subscheme of Projective Space of dimension 2 over Rational Field
            defined by:
                x,
                y,
                z

        ::

            sage: P1.<x,y,z> = ProjectiveSpace(RR,2)
            sage: P2.<t,u,v,w> = ProjectiveSpace(RR,3)
            sage: H = Hom(P1,P2)
            sage: h = H([y^3*z^3, x^3*z^3, y^3*z^3, x^2*y^2*z^2])
            sage: h.indeterminacy_locus()
            Closed subscheme of Projective Space of dimension 2 over Real Field with
            53 bits of precision defined by:
              y^3*z^3,
              x^3*z^3,
              y^3*z^3,
              x^2*y^2*z^2

        If defining polynomials are not normalized, output scheme will not be normalized::

            sage: P.<x,y,z>=ProjectiveSpace(QQ,2)
            sage: H=End(P)
            sage: f=H([x*x^2,x*y^2,x*z^2])
            sage: f.indeterminacy_locus()
            Closed subscheme of Projective Space of dimension 2 over Rational Field
            defined by:
              x^3,
              x*y^2,
              x*z^2
        """
        from sage.schemes.projective.projective_space import is_ProjectiveSpace
        dom = self.domain()
        if not is_ProjectiveSpace(dom):
            raise NotImplementedError("not implemented for subschemes")
        defPolys = self.defining_polynomials()
        locus = dom.subscheme(defPolys)
        if locus.dimension() < 0:
            locus = dom.subscheme(dom.gens())
        return locus 

    def indeterminacy_points(self, F=None):
        r"""
        Return the indeterminacy locus of this map defined over ``F``.

        Only for rational maps on projective space. Returns the set of points in projective space at which all of the defining polynomials of the rational map simultaneously vanish.

        INPUT:

        - ``F`` - a field (optional).

        OUTPUT:

        - indeterminacy points of the map defined over ``F``, provided the indeterminacy scheme is 0-dimensional.

        EXAMPLES::

            sage: P.<x,y,z> = ProjectiveSpace(QQ,2)
            sage: H = End(P)
            sage: f = H([x*z-y*z, x^2-y^2, z^2])
            sage: f.indeterminacy_points()
            [(-1 : 1 : 0), (1 : 1 : 0)]

        ::

            sage: P1.<x,y,z> = ProjectiveSpace(RR,2)
            sage: P2.<t,u,v,w> = ProjectiveSpace(RR,3)
            sage: H = Hom(P1,P2)
            sage: h = H([x+y, y, z+y, y])
            sage: h.indeterminacy_points()
            []
            sage: g = H([y^3*z^3, x^3*z^3, y^3*z^3, x^2*y^2*z^2])
            sage: g.indeterminacy_points()
            Traceback (most recent call last):
            ...
            ValueError: indeterminacy scheme is not dimension 0

        ::

            sage: P.<x,y,z> = ProjectiveSpace(QQ,2)
            sage: H = End(P)
            sage: f = H([x^2+y^2, x*z, x^2+y^2])
            sage: f.indeterminacy_points()
            [(0 : 0 : 1)]
            sage: R.<t> = QQ[]
            sage: K.<a> = NumberField(t^2+1)
            sage: f.indeterminacy_points(F=K)
            [(-a : 1 : 0), (0 : 0 : 1), (a : 1 : 0)]
            sage: set_verbose(None)
            sage: f.indeterminacy_points(F=QQbar)
            [(-1*I : 1 : 0), (0 : 0 : 1), (1*I : 1 : 0)]

        ::

            sage: set_verbose(None)
            sage: K.<t>=FunctionField(QQ)
            sage: P.<x,y,z>=ProjectiveSpace(K,2)
            sage: H=End(P)
            sage: f=H([x^2-t^2*y^2,y^2-z^2,x^2-t^2*z^2])
            sage: f.indeterminacy_points()
            [(-t : -1 : 1), (-t : 1 : 1), (t : -1 : 1), (t : 1 : 1)]

        ::

            sage: set_verbose(None)
            sage: P.<x,y,z>=ProjectiveSpace(Qp(3),2)
            sage: H=End(P)
            sage: f=H([x^2-7*y^2,y^2-z^2,x^2-7*z^2])
            sage: f.indeterminacy_points()
            [(2 + 3 + 3^2 + 2*3^3 + 2*3^5 + 2*3^6 + 3^8 + 3^9 + 2*3^11 + 3^15 +
            2*3^16 + 3^18 + O(3^20) : 1 + O(3^20) : 1 + O(3^20)),
            (2 + 3 + 3^2 + 2*3^3 + 2*3^5 + 2*3^6 + 3^8 + 3^9 + 2*3^11 + 3^15 +
            2*3^16 + 3^18 + O(3^20) : 2 + 2*3 + 2*3^2 + 2*3^3 + 2*3^4 + 2*3^5 +
            2*3^6 + 2*3^7 + 2*3^8 + 2*3^9 + 2*3^10 + 2*3^11 + 2*3^12 + 2*3^13 +
            2*3^14 + 2*3^15 + 2*3^16 + 2*3^17 + 2*3^18 + 2*3^19 + O(3^20) : 1 +
            O(3^20)),
             (1 + 3 + 3^2 + 2*3^4 + 2*3^7 + 3^8 + 3^9 + 2*3^10 + 2*3^12 + 2*3^13 +
            2*3^14 + 3^15 + 2*3^17 + 3^18 + 2*3^19 + O(3^20) : 1 + O(3^20) : 1 +
            O(3^20)),
             (1 + 3 + 3^2 + 2*3^4 + 2*3^7 + 3^8 + 3^9 + 2*3^10 + 2*3^12 + 2*3^13 +
            2*3^14 + 3^15 + 2*3^17 + 3^18 + 2*3^19 + O(3^20) : 2 + 2*3 + 2*3^2 +
            2*3^3 + 2*3^4 + 2*3^5 + 2*3^6 + 2*3^7 + 2*3^8 + 2*3^9 + 2*3^10 + 2*3^11
            + 2*3^12 + 2*3^13 + 2*3^14 + 2*3^15 + 2*3^16 + 2*3^17 + 2*3^18 + 2*3^19
            + O(3^20) : 1 + O(3^20))]
        """
        if F is None:
            fcn = self
        else:
            if not F.is_field():
                raise NotImplementedError("indeterminacy points only implemented for fields")
            fcn = self.change_ring(F)
        indScheme = fcn.indeterminacy_locus()
        if indScheme.dimension() > 0:
            raise ValueError("indeterminacy scheme is not dimension 0")
        # Other error checking is in indeterminacy_locus
        indPoints = indScheme.rational_points()
        return indPoints

class SchemeMorphism_polynomial_projective_space_finite_field(SchemeMorphism_polynomial_projective_space_field):

    def _fast_eval(self, x):
        """
        Evaluate projective morphism at point described by x.

        EXAMPLES::

            sage: P.<x,y,z> = ProjectiveSpace(GF(7),2)
            sage: H = Hom(P,P)
            sage: f = H([x^2+y^2, y^2, z^2 + y*z])
            sage: f._fast_eval([1,1,1])
            [2, 1, 2]
        """
        if self._is_prime_finite_field:
            p = self.base_ring().characteristic()
            P = [Integer(f(*x)) % p for f in self._fastpolys]
        else:
            P = [f(*x) for f in self._fastpolys]
        return P

    def orbit_structure(self, P):
        r"""
        Return the pair `[m,n]` where `m` is the preperiod and `n`
        is the period of the point ``P`` by this map.

        Every point is preperiodic over a finite field so every point
        will be preperiodic.

        INPUT:

        - ``P`` -- a point in the domain of this map.

        OUTPUT:

        - a list `[m,n]` of integers.

        EXAMPLES::

            sage: P.<x,y,z> = ProjectiveSpace(GF(5),2)
            sage: H = Hom(P,P)
            sage: f = H([x^2 + y^2,y^2, z^2 + y*z])
            sage: f.orbit_structure(P(2,1,2))
            [0, 6]

        ::

            sage: P.<x,y,z> = ProjectiveSpace(GF(7),2)
            sage: X = P.subscheme(x^2-y^2)
            sage: H = Hom(X,X)
            sage: f = H([x^2, y^2, z^2])
            sage: f.orbit_structure(X(1,1,2))
            [0, 2]

        ::

            sage: P.<x,y> = ProjectiveSpace(GF(13),1)
            sage: H = Hom(P,P)
            sage: f = H([x^2 - y^2, y^2])
            sage: f.orbit_structure(P(3,4))
            [2, 3]
        """
        return(P.orbit_structure(self))

    def cyclegraph(self):
        r"""
        Return the digraph of all orbits of this map.

        Over a finite field this is a finite graph. For subscheme domains, only points
        on the subscheme whose image are also on the subscheme are in the digraph.

        OUTPUT:

        - a digraph

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(GF(13),1)
            sage: H = Hom(P,P)
            sage: f = H([x^2-y^2, y^2])
            sage: f.cyclegraph()
            Looped digraph on 14 vertices

        ::

            sage: P.<x,y,z> = ProjectiveSpace(GF(5^2,'t'),2)
            sage: H = Hom(P,P)
            sage: f = H([x^2+y^2, y^2, z^2+y*z])
            sage: f.cyclegraph()
            Looped digraph on 651 vertices

        ::

            sage: P.<x,y,z> = ProjectiveSpace(GF(7),2)
            sage: X = P.subscheme(x^2-y^2)
            sage: H = Hom(X,X)
            sage: f = H([x^2, y^2, z^2])
            sage: f.cyclegraph()
            Looped digraph on 15 vertices
        """
        if self.domain() != self.codomain():
            raise NotImplementedError("domain and codomain must be equal")
        V = []
        E = []
        from sage.schemes.projective.projective_space import is_ProjectiveSpace
        if is_ProjectiveSpace(self.domain()) is True:
            for P in self.domain():
                V.append(P)
                Q = self(P)
                Q.normalize_coordinates()
                E.append([Q])
        else:
            X = self.domain()
            for P in X.ambient_space():
                try:
                    XP = X.point(P)
                    V.append(XP)
                    Q = self(XP)
                    Q.normalize_coordinates()
                    E.append([Q])
                except TypeError:  # not a point on the scheme
                    pass
        from sage.graphs.digraph import DiGraph
        g = DiGraph(dict(zip(V, E)), loops=True)
        return g

    def possible_periods(self, return_points=False):
        r"""
        Returns the list of possible minimal periods of a periodic point
        over `\QQ` and (optionally) a point in each cycle.

        REFERENCES:

        .. [Hutz-gr] \B. Hutz. Good reduction of periodic points, Illinois Journal of
           Mathematics 53 (Winter 2009), no. 4, 1109-1126.

        ALGORITHM:

        See [Hutz-gr]_.

        INPUT:

        - ``return_points`` - Boolean (optional) - a value of True returns the points as well as the possible periods.

        OUTPUT:

        - a list of positive integers, or a list of pairs of projective points and periods if ``flag`` is 1.

        Examples::

            sage: P.<x,y> = ProjectiveSpace(GF(23),1)
            sage: H = End(P)
            sage: f = H([x^2-2*y^2, y^2])
            sage: f.possible_periods()
            [1, 5, 11, 22, 110]

        ::

            sage: P.<x,y> = ProjectiveSpace(GF(13),1)
            sage: H = End(P)
            sage: f = H([x^2-y^2, y^2])
            sage: sorted(f.possible_periods(True))
            [[(0 : 1), 2], [(1 : 0), 1], [(3 : 1), 3], [(3 : 1), 36]]

        ::

            sage: PS.<x,y,z> = ProjectiveSpace(2,GF(7))
            sage: H = End(PS)
            sage: f = H([-360*x^3 + 760*x*z^2, y^3 - 604*y*z^2 + 240*z^3, 240*z^3])
            sage: f.possible_periods()
            [1, 2, 4, 6, 12, 14, 28, 42, 84]

        .. TODO::

            - do not return duplicate points

            - improve hash to reduce memory of pointtable
        """
        return _fast_possible_periods(self,return_points)

    def automorphism_group(self, **kwds):
        r"""
        Return the subgroup of `PGL2` that is the automorphism group of this map.

        Only for dimension 1. The automorphism group is the set of `PGL2` elements that
        fixed the map under conjugation. See [FMV]_ for the algorithm.

        INPUT:

        keywords:

        - ``absolute``-- Boolean - True returns the absolute automorphism group and a field of definition. default: False (optional)

        - ``iso_type`` -- Boolean - True returns the isomorphism type of the automorphism group. default: False (optional)

        - ``return_functions``-- Boolean - True returns elements as linear fractional transformations.
            False returns elements as `PGL2` matrices. default: False. (optional)

        OUTPUT:

        - list - elements of automorphism group.

        AUTHORS:

        - Original algorithm written by Xander Faber, Michelle Manes, Bianca Viray

        - Modified by Joao Alberto de Faria, Ben Hutz, Bianca Thompson

        EXAMPLES::

            sage: R.<x,y> = ProjectiveSpace(GF(7^3,'t'),1)
            sage: H = End(R)
            sage: f = H([x^2-y^2, x*y])
            sage: f.automorphism_group()
            [
            [1 0]  [6 0]
            [0 1], [0 1]
            ]

        ::

            sage: R.<x,y> = ProjectiveSpace(GF(3^2,'t'),1)
            sage: H = End(R)
            sage: f = H([x^3,y^3])
            sage: f.automorphism_group(return_functions=True, iso_type=True) # long time
            ([x, x/(x + 1), x/(2*x + 1), 2/(x + 2), (2*x + 1)/(2*x), (2*x + 2)/x,
            1/(2*x + 2), x + 1, x + 2, x/(x + 2), 2*x/(x + 1), 2*x, 1/x, 2*x + 1,
            2*x + 2, ((t + 2)*x + t + 2)/((2*t + 1)*x + t + 2), (t*x + 2*t)/(t*x +
            t), 2/x, (x + 1)/(x + 2), (2*t*x + t)/(t*x), (2*t + 1)/((2*t + 1)*x +
            2*t + 1), ((2*t + 1)*x + 2*t + 1)/((2*t + 1)*x), t/(t*x + 2*t), (2*x +
            1)/(x + 1)], 'PGL(2,3)')

        ::

            sage: R.<x,y> = ProjectiveSpace(GF(2^5,'t'),1)
            sage: H = End(R)
            sage: f=H([x^5,y^5])
            sage: f.automorphism_group(return_functions=True, iso_type=True)
            ([x, 1/x], 'Cyclic of order 2')

            ::

            sage: R.<x,y> = ProjectiveSpace(GF(3^4,'t'),1)
            sage: H = End(R)
            sage: f=H([x^2+25*x*y+y^2, x*y+3*y^2])
            sage: f.automorphism_group(absolute=True)
            [Univariate Polynomial Ring in w over Finite Field in b of size 3^4,
             [
            [1 0]
            [0 1]
            ]]
        """
        absolute = kwds.get('absolute', False)
        iso_type = kwds.get('iso_type', False)
        return_functions=kwds.get('return_functions', False)

        if self.domain().dimension_relative()!=1:
            raise NotImplementedError("must be dimension 1")
        else:
            f = self.dehomogenize(1)
            z = f[0].parent().gen()
        if f[0].denominator()!=1:
            F = (f[0].numerator().polynomial(z))/f[0].denominator().polynomial(z)
        else:
            F = f[0].numerator().polynomial(z)
        from .endPN_automorphism_group import automorphism_group_FF
        return(automorphism_group_FF(F, absolute, iso_type, return_functions))


