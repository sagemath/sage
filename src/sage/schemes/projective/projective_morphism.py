r"""
Morphisms on projective schemes

This module defines morphisms from projective schemes. A morphism from a
projective scheme to a projective scheme is defined by homogeneous polynomials
of the same degree that define what the morphism does on points in the ambient
projective space. A morphism from a projective scheme to an affine scheme is
determined by rational function, that is, quotients of homogeneous polynomials
of the same degree.

EXAMPLES::

    sage: P2.<x0,x1,x2> = ProjectiveSpace(QQ, 2)
    sage: A2.<x,y> = AffineSpace(QQ, 2)
    sage: P2.hom([x0, x1, x1 + x2], P2)
    Scheme endomorphism of Projective Space of dimension 2 over Rational Field
      Defn: Defined on coordinates by sending (x0 : x1 : x2) to (x0 : x1 : x1 + x2)
    sage: P2.hom([x1/x0, (x1 + x2)/x0], A2)
    Scheme morphism:
      From: Projective Space of dimension 2 over Rational Field
      To:   Affine Space of dimension 2 over Rational Field
      Defn: Defined on coordinates by sending (x0 : x1 : x2) to (x1/x0, (x1 + x2)/x0)

AUTHORS:

- David Kohel, William Stein: initial version

- William Stein (2006-02-11): fixed bug where P(0,0,0) was allowed as
  a projective point

- Volker Braun (2011-08-08): renamed classes, more documentation, misc
  cleanups

- Ben Hutz (2013-03): iteration functionality and new directory structure
  for affine/projective, height functionality

- Brian Stout, Ben Hutz (2013-11): added minimal model functionality

- Dillon Rose (2014-01): speed enhancements

- Ben Hutz (2015-11): iteration of subschemes

- Kwankyu Lee (2020-02): added indeterminacy_locus() and image()

- Kwankyu Lee (2022-05): added graph(), projective_degrees(), and degree()
"""

# ****************************************************************************
#       Copyright (C) 2011 Volker Braun <vbraun.name@gmail.com>
#       Copyright (C) 2006 David Kohel <kohel@maths.usyd.edu.au>
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
# ****************************************************************************

import sys

import sage.rings.abc

from sage.arith.functions import lcm
from sage.arith.misc import GCD as gcd
from sage.categories.fields import Fields
from sage.categories.finite_fields import FiniteFields
from sage.categories.homset import Hom, End
from sage.categories.number_fields import NumberFields
from sage.misc.cachefunc import cached_method
from sage.misc.lazy_attribute import lazy_attribute
from sage.misc.lazy_import import lazy_import
from sage.misc.misc_c import prod
from sage.rings.finite_rings.finite_field_base import FiniteField
from sage.rings.fraction_field import FractionField
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.quotient_ring import QuotientRing_generic
from sage.rings.rational_field import QQ
from sage.schemes.generic.morphism import SchemeMorphism_polynomial

lazy_import('sage.dynamics.arithmetic_dynamics.generic_ds', 'DynamicalSystem')
lazy_import('sage.dynamics.arithmetic_dynamics.projective_ds',
            ['DynamicalSystem_projective', 'DynamicalSystem_projective_field',
             'DynamicalSystem_projective_finite_field'])
lazy_import('sage.rings.algebraic_closure_finite_field', 'AlgebraicClosureFiniteField_generic')
lazy_import('sage.rings.number_field.number_field_ideal', 'NumberFieldFractionalIdeal')
lazy_import('sage.rings.padics.padic_base_generic', 'pAdicGeneric')
lazy_import('sage.rings.padics.padic_valuation', 'pAdicValuation_base')


_NumberFields = NumberFields()
_FiniteFields = FiniteFields()
_Fields = Fields()


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

    An example of a morphism between projective plane curves (see :issue:`10297`)::

        sage: # needs sage.schemes
        sage: P2.<x,y,z> = ProjectiveSpace(QQ, 2)
        sage: f = x^3 + y^3 + 60*z^3
        sage: g = y^2*z - (x^3 - 6400*z^3/3)
        sage: C = Curve(f)
        sage: E = Curve(g)
        sage: xbar,ybar,zbar = C.coordinate_ring().gens()
        sage: H = C.Hom(E)
        sage: H([zbar, xbar - ybar, -(xbar+ybar)/80])
        Scheme morphism:
          From: Projective Plane Curve over Rational Field defined by x^3 + y^3 + 60*z^3
          To:   Projective Plane Curve over Rational Field defined by -x^3 + y^2*z + 6400/3*z^3
          Defn: Defined on coordinates by sending (x : y : z) to
                (z : x - y : -1/80*x - 1/80*y)

    A more complicated example::

        sage: P2.<x,y,z> = ProjectiveSpace(2, QQ)
        sage: P1 = P2.subscheme(x - y)
        sage: H12 = P1.Hom(P2)
        sage: H12([x^2, x*z, z^2])
        Scheme morphism:
          From: Closed subscheme of Projective Space of dimension 2 over Rational Field
                defined by: x - y
          To:   Projective Space of dimension 2 over Rational Field
          Defn: Defined on coordinates by sending (x : y : z) to (x^2 : x*z : z^2)

    We illustrate some error checking::

        sage: R.<x,y> = QQ[]
        sage: P1 = ProjectiveSpace(R)
        sage: H = P1.Hom(P1)
        sage: f = H([x - y, x*y])
        Traceback (most recent call last):
        ...
        ValueError: polys (=[x - y, x*y]) must be of the same degree

        sage: H([x - 1, x*y + x])
        Traceback (most recent call last):
        ...
        ValueError: polys (=[x - 1, x*y + x]) must be homogeneous

        sage: H([exp(x), exp(y)])                                                       # needs sage.symbolic
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
        sage: f(f(f(X)))                                                                # needs sage.libs.singular
        Closed subscheme of Projective Space of dimension 2 over Rational Field
         defined by:
          y - z

    ::

        sage: P.<x,y,z,w> = ProjectiveSpace(QQ, 3)
        sage: H = End(P)
        sage: f = H([(x-2*y)^2, (x-2*z)^2, (x-2*w)^2, x^2])
        sage: f(P.subscheme([x,y,z]))                                                   # needs sage.libs.singular
        Closed subscheme of Projective Space of dimension 3 over Rational Field
         defined by:
          w,
          y,
          x
    """
    def __init__(self, parent, polys, check=True):
        """
        Initialize.

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
            Scheme endomorphism of Closed subscheme of Projective Space of dimension 2
             over Univariate Polynomial Ring in t over Rational Field defined by: x
              Defn: Defined on coordinates by sending (x : y : z) to
                    (x^2 : t*y^2 : x*z)

        When elements of the quotient ring is used, they are reduced::

            sage: # needs sage.rings.real_mpfr
            sage: P.<x,y,z> = ProjectiveSpace(CC, 2)
            sage: X = P.subscheme([x - y])
            sage: u,v,w = X.coordinate_ring().gens()                                    # needs sage.libs.singular
            sage: H = End(X)
            sage: H([u^2, v^2, w*u])                                                    # needs sage.libs.singular
            Scheme endomorphism of Closed subscheme of Projective Space of dimension 2
             over Complex Field with 53 bits of precision defined by: x - y
              Defn: Defined on coordinates by sending (x : y : z) to
                    (y^2 : y^2 : y*z)
        """
        if check:
            try:  # polys might be in a quotient ring
                polys = [f.lift() for f in polys]
            except (TypeError, AttributeError):
                pass

            try:
                source_ring = parent.domain().ambient_space().coordinate_ring()
                K = FractionField(source_ring)

                try:
                    polys = [K(f) for f in polys]
                except TypeError:
                    raise TypeError("polys (=%s) must be elements of %s" % (polys, source_ring))

                if parent.codomain().is_projective():
                    degs = []
                    l = 1
                    for f in polys:
                        num = f.numerator()
                        den = f.denominator()
                        if not num.is_homogeneous() or not den.is_homogeneous():
                            raise ValueError("polys (={}) must be homogeneous".format(polys))

                        if not num.is_zero():
                            l *= den
                            degs.append(num.degree() - den.degree())

                    d = degs[0]
                    if not all(d == deg for deg in degs[1:]):
                        raise ValueError("polys (={}) must be of the same degree".format(polys))

                    polys = [(l * f).numerator() for f in polys]
                elif parent.codomain().is_affine():
                    for f in polys:
                        num = f.numerator()
                        den = f.denominator()
                        if not (num.is_homogeneous() and
                                den.is_homogeneous() and
                                num.degree() == den.degree()):
                            raise ValueError("polys (={}) must be quotients of "
                                             "homogeneous polynomials of the same degree".format(polys))
                check = False
            except (NotImplementedError, TypeError, AttributeError):
                pass

        SchemeMorphism_polynomial.__init__(self, parent, polys, check)

        R = polys[0].base_ring()
        self._is_prime_finite_field = isinstance(R, FiniteField) and R.is_prime_field()

    def __call__(self, x, check=True):
        r"""
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

        - ``x`` -- a point or subscheme in domain of this map

        - ``check`` -- boolean; if ``False`` assume that ``x`` is a point

        EXAMPLES::

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: H = End(P)
            sage: f = H([x^2 + y^2, y^2, z^2 + y*z])
            sage: f(P([1,1,1]))
            (1 : 1/2 : 1)

        ::

            sage: PS.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: P1.<u,v> = ProjectiveSpace(QQ, 1)
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
            sage: X = PS.subscheme([z^2 + y*w])
            sage: f(X)                                                                  # needs sage.libs.singular
            Closed subscheme of Projective Space of dimension 3 over Rational Field
            defined by:
              x*z - w^2

        ::

            sage: PS.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: P1.<u,v> = ProjectiveSpace(ZZ, 1)
            sage: H = End(PS)
            sage: f = H([x^2, y^2, z^2])
            sage: X = P1.subscheme([u - v])
            sage: f(X)                                                                  # needs sage.libs.singular
            Traceback (most recent call last):
            ...
            TypeError: subscheme must be in ambient space of domain of map

        ::

            sage: PS.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: P1.<u,v> = ProjectiveSpace(ZZ, 1)
            sage: H = End(P1)
            sage: f = H([u^2, v^2])
            sage: f([u - v])                                                            # needs sage.libs.singular
            Closed subscheme of Projective Space of dimension 1 over Integer Ring defined by:
              u - v
            sage: X = PS.subscheme([x - z])
            sage: f([x - z])
            Traceback (most recent call last):
            ...
            TypeError: [x - z] fails to convert into the map's domain Projective Space of
            dimension 1 over Integer Ring, but a `pushforward` method is not properly implemented

        TESTS:

        Check that :issue:`32209` is fixed::

            sage: S.<x,y> = ProjectiveSpace(ZZ, 1)
            sage: T.<u,v> = ProjectiveSpace(ZZ, 1)
            sage: h = T.hom([u^2 + v^2, u*v], S); h
            Scheme morphism:
              From: Projective Space of dimension 1 over Integer Ring
              To:   Projective Space of dimension 1 over Integer Ring
              Defn: Defined on coordinates by sending (u : v) to
                    (u^2 + v^2 : u*v)

            sage: # needs sage.rings.finite_rings
            sage: F.<a> = GF(4)
            sage: P = T(F)(1, a)
            sage: h(P)                                                                  # needs sage.libs.singular
            (1 : 1)
            sage: h(P).domain()
            Spectrum of Finite Field in a of size 2^2
            sage: h.change_ring(F)(P)
            (1 : 1)
        """
        from sage.schemes.projective.projective_point import SchemeMorphism_point_projective_ring
        if check:
            from sage.schemes.projective.projective_subscheme import AlgebraicScheme_subscheme_projective
            if isinstance(x, SchemeMorphism_point_projective_ring):
                if self.domain() != x.codomain():
                    try:
                        x = self.domain()(x)
                    except (TypeError, NotImplementedError):
                        raise TypeError("%s fails to convert into the map's domain %s, but a `pushforward` method is not properly implemented" % (x, self.domain()))
                # else pass it onto the eval below
            elif isinstance(x, AlgebraicScheme_subscheme_projective):
                return x._forward_image(self)  # call subscheme eval
            else:  # not a projective point or subscheme
                try:
                    x = self.domain()(x)
                except (TypeError, NotImplementedError):
                    try:
                        x = self.domain().subscheme(x)
                        return x._forward_image(self)  # call subscheme eval
                    except (TypeError, NotImplementedError):
                        raise TypeError("%s fails to convert into the map's domain %s, but a `pushforward` method is not properly implemented" % (x, self.domain()))

        R = x.domain().coordinate_ring()
        if R is self.base_ring():
            P = self._fast_eval(x._coords)
        else:
            P = [f(x._coords) for f in self._polys]
        return self.codomain().point_homset(R)(P, check=check)

    @lazy_attribute
    def _fastpolys(self):
        """
        Lazy attribute for fast_callable polynomials for this map.

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: H = Hom(P,P)
            sage: f = H([x^2 + y^2, y^2])
            sage: [g.op_list() for g in f._fastpolys]
            [[('load_const', 0), ('load_const', 1), ('load_arg', ...), ('ipow', 2), 'mul', 'add', ('load_const', 1), ('load_arg', ...), ('ipow', 2), 'mul', 'add', 'return'],
             [('load_const', 0), ('load_const', 1), ('load_arg', 1), ('ipow', 2), 'mul', 'add', 'return']]
        """
        from sage.ext.fast_callable import fast_callable

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

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: H = Hom(P, P)
            sage: f = H([x^2 + y^2, y^2, z^2 + y*z])
            sage: f._fast_eval([1,1,1])
            [2, 1, 2]

            ::

            sage: T.<z> = LaurentSeriesRing(ZZ)
            sage: P.<x,y> = ProjectiveSpace(T, 1)
            sage: H = End(P)
            sage: f = H([x^2 + x*y, y^2])
            sage: Q = P(z,1)
            sage: f._fast_eval(list(Q))
            [z + z^2, 1]

            ::

            sage: # needs sage.libs.pari sage.rings.real_mpfr
            sage: T.<z> = PolynomialRing(CC)
            sage: I = T.ideal(z^3)
            sage: P.<x,y> = ProjectiveSpace(T.quotient_ring(I), 1)
            sage: H = End(P)
            sage: f = H([x^2 + x*y, y^2])
            sage: Q = P(z^2, 1)
            sage: f._fast_eval(list(Q))
            [zbar^2, 1.00000000000000]

            ::

            sage: # needs sage.rings.real_mpfr
            sage: T.<z> = LaurentSeriesRing(CC)
            sage: R.<t> = PolynomialRing(T)
            sage: P.<x,y> = ProjectiveSpace(R,1)
            sage: H = End(P)
            sage: f = H([x^2 + x*y, y^2])
            sage: Q = P(t^2, z)
            sage: f._fast_eval(list(Q))
            [t^4 + z*t^2, z^2]
        """
        P = [f(*x) for f in self._fastpolys]
        return P

    def __eq__(self, right):
        """
        Test the equality of two projective morphisms.

        INPUT:

        - ``right`` -- a map on projective space

        OUTPUT:

        ``True`` if ``self`` and ``right`` define the same projective map.
        ``False`` otherwise.

        EXAMPLES::

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: H = Hom(P, P)
            sage: f = H([x^2 - 2*x*y + z*x, z^2 - y^2, 5*z*y])
            sage: g = H([x^2, y^2, z^2])
            sage: f == g
            False

        ::

            sage: # needs sage.rings.real_mpfr
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
        return all(self._polys[i] * right._polys[j] == self._polys[j] * right._polys[i]
                   for i in range(n) for j in range(i + 1, n))

    def __ne__(self, right):
        """
        Test the inequality of two projective morphisms.

        INPUT:

        - ``right`` -- a map on projective space

        OUTPUT:

        ``True`` if ``self`` and ``right`` define different projective maps.
        ``False`` otherwise.

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: H = Hom(P, P)
            sage: f = H([x^3 - 2*x^2*y, 5*x*y^2])
            sage: g = f.change_ring(GF(7))
            sage: f != g
            True

        ::

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: H = Hom(P, P)
            sage: f = H([x^2 - 2*x*y + z*x, z^2 - y^2, 5*z*y])
            sage: f != f
            False
        """
        if not isinstance(right, SchemeMorphism_polynomial):
            return True
        if self.parent() != right.parent():
            return True
        n = len(self._polys)
        return any(self._polys[i] * right._polys[j] != self._polys[j] * right._polys[i]
                   for i in range(n) for j in range(i + 1, n))

    def _matrix_times_polymap_(self, mat, h):
        """
        Multiply the morphism on the left by a matrix ``mat``.

        INPUT:

        - ``mat`` -- a matrix

        OUTPUT: a scheme morphism given by ``self*mat``

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(ZZ, 1)
            sage: H = Hom(P,P)
            sage: f = H([x^2 + y^2, y^2])
            sage: matrix([[1,2], [0,1]]) * f                                            # needs sage.modules
            Scheme endomorphism of Projective Space of dimension 1 over Integer Ring
              Defn: Defined on coordinates by sending (x : y) to
                    (x^2 + 3*y^2 : y^2)

        ::

            sage: # needs sage.rings.number_field
            sage: R.<x> = PolynomialRing(QQ)
            sage: K.<i> = NumberField(x^2 + 1)
            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: H = Hom(P, P)
            sage: f = H([1/3*x^2 + 1/2*y^2, y^2])
            sage: matrix([[i,0], [0,i]]) * f                                            # needs sage.modules
            Scheme endomorphism of Projective Space of dimension 1 over
             Number Field in i with defining polynomial x^2 + 1
              Defn: Defined on coordinates by sending (x : y) to
                    ((1/3*i)*x^2 + (1/2*i)*y^2 : i*y^2)
        """
        from sage.modules.free_module_element import vector

        if not mat.is_square():
            raise ValueError("matrix must be square")
        if mat.ncols() != self.codomain().ngens():
            raise ValueError("matrix size is incompatible")
        F = mat * vector(list(self))
        if isinstance(self, DynamicalSystem):
            return h(list(F)).as_dynamical_system()
        return h(list(F))

    def _polymap_times_matrix_(self, mat, h):
        """
        Multiply the morphism on the right by a matrix ``mat``.

        INPUT:

        - ``mat`` -- a matrix

        OUTPUT: a scheme morphism given by ``mat*self``

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(ZZ, 1)
            sage: H = Hom(P, P)
            sage: f = H([x^2 + y^2, y^2])
            sage: f * matrix([[1,2], [0,1]])                                            # needs sage.modules
            Scheme endomorphism of Projective Space of dimension 1 over Integer Ring
              Defn: Defined on coordinates by sending (x : y) to
                    (x^2 + 4*x*y + 5*y^2 : y^2)

        ::

            sage: # needs sage.rings.number_field
            sage: R.<x> = PolynomialRing(QQ)
            sage: K.<i> = NumberField(x^2 + 1)
            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: H = Hom(P, P)
            sage: f = H([1/3*x^2 + 1/2*y^2, y^2])
            sage: f * matrix([[i,0], [0,i]])                                            # needs sage.modules
            Scheme endomorphism of Projective Space of dimension 1 over
             Number Field in i with defining polynomial x^2 + 1
              Defn: Defined on coordinates by sending (x : y) to
                    (-1/3*x^2 - 1/2*y^2 : -y^2)
        """
        from sage.modules.free_module_element import vector
        if not mat.is_square():
            raise ValueError("matrix must be square")
        if mat.nrows() != self.domain().ngens():
            raise ValueError("matrix size is incompatible")
        X = mat * vector(self[0].parent().gens())
        F = vector(self._polys)
        F = F(list(X))
        if isinstance(self, DynamicalSystem):
            return h(list(F)).as_dynamical_system()
        return h(list(F))

    def as_dynamical_system(self):
        """
        Return this endomorphism as a :class:`DynamicalSystem_projective`.

        OUTPUT: :class:`DynamicalSystem_projective`

        EXAMPLES::

            sage: P.<x,y,z> = ProjectiveSpace(ZZ, 2)
            sage: H = End(P)
            sage: f = H([x^2, y^2, z^2])
            sage: type(f.as_dynamical_system())                                         # needs sage.schemes
            <class 'sage.dynamics.arithmetic_dynamics.projective_ds.DynamicalSystem_projective'>

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: H = End(P)
            sage: f = H([x^2 - y^2, y^2])
            sage: type(f.as_dynamical_system())                                         # needs sage.schemes
            <class 'sage.dynamics.arithmetic_dynamics.projective_ds.DynamicalSystem_projective_field'>

        ::

            sage: P.<x,y> = ProjectiveSpace(GF(5), 1)
            sage: H = End(P)
            sage: f = H([x^2, y^2])
            sage: type(f.as_dynamical_system())                                         # needs sage.schemes
            <class 'sage.dynamics.arithmetic_dynamics.projective_ds.DynamicalSystem_projective_finite_field'>

        ::

            sage: P.<x,y> = ProjectiveSpace(RR, 1)
            sage: f = DynamicalSystem([x^2 + y^2, y^2], P)                              # needs sage.schemes
            sage: g = f.as_dynamical_system()                                           # needs sage.schemes
            sage: g is f                                                                # needs sage.schemes
            True
        """
        if isinstance(self, DynamicalSystem):
            return self
        if not self.is_endomorphism():
            raise TypeError("must be an endomorphism")
        R = self.base_ring()
        if R not in _Fields:
            return DynamicalSystem_projective(list(self), self.domain())
        if isinstance(R, FiniteField):
            return DynamicalSystem_projective_finite_field(list(self), self.domain())
        return DynamicalSystem_projective_field(list(self), self.domain())

    def scale_by(self, t):
        """
        Scale each coordinate by a factor of ``t``.

        A :exc:`TypeError` occurs if the point is not in the coordinate ring
        of the parent after scaling.

        INPUT:

        - ``t`` -- a ring element

        OUTPUT: none

        EXAMPLES::

            sage: A.<x,y> = ProjectiveSpace(QQ, 1)
            sage: H = Hom(A, A)
            sage: f = H([x^3 - 2*x*y^2, x^2*y])
            sage: f.scale_by(1/x)
            sage: f
            Scheme endomorphism of Projective Space of dimension 1 over Rational Field
              Defn: Defined on coordinates by sending (x : y) to (x^2 - 2*y^2 : x*y)

        ::

            sage: R.<t> = PolynomialRing(QQ)
            sage: P.<x,y> = ProjectiveSpace(R, 1)
            sage: H = Hom(P,P)
            sage: f = H([3/5*x^2, 6*y^2])
            sage: f.scale_by(5/3*t); f
            Scheme endomorphism of Projective Space of dimension 1 over
             Univariate Polynomial Ring in t over Rational Field
              Defn: Defined on coordinates by sending (x : y) to (t*x^2 : 10*t*y^2)

        ::

            sage: P.<x,y,z> = ProjectiveSpace(GF(7), 2)
            sage: X = P.subscheme(x^2 - y^2)
            sage: H = Hom(X, X)
            sage: f = H([x^2, y^2, z^2])
            sage: f.scale_by(x - y); f                                                  # needs sage.libs.singular
            Scheme endomorphism of Closed subscheme of Projective Space of dimension 2
             over Finite Field of size 7 defined by: x^2 - y^2
              Defn: Defined on coordinates by sending (x : y : z) to
                    (x*y^2 - y^3 : x*y^2 - y^3 : x*z^2 - y*z^2)
        """
        if t == 0:
            raise ValueError("Cannot scale by 0")
        R = self.domain().coordinate_ring()
        if isinstance(R, QuotientRing_generic):
            phi = R._internal_coerce_map_from(self.domain().ambient_space().coordinate_ring())
            for i in range(self.codomain().ambient_space().dimension_relative() + 1):
                new_polys = [phi(u * t).lift() for u in self]
        else:
            for i in range(self.codomain().ambient_space().dimension_relative() + 1):
                new_polys = [R(u * t) for u in self]
        self._polys = tuple(new_polys)

    def normalize_coordinates(self, **kwds):
        """
        Ensures that this morphism has integral coefficients.
        If the coordinate ring has a GCD, then it ensures that the
        coefficients have no common factor.

        It also makes the leading coefficients of the first polynomial
        positive (if positive has meaning in the coordinate ring).
        This is done in place.

        When ``ideal`` or ``valuation`` is specified, normalization occurs
        with respect to the absolute value defined by the ``ideal`` or
        ``valuation``. That is, the coefficients are scaled such that
        one coefficient has absolute value 1 while the others have
        absolute value less than or equal to 1.
        Only supported when the base ring is a number field.

        INPUT: keyword arguments:

        - ``ideal`` -- (optional) a prime ideal of the base ring of this
          morphism

        - ``valuation`` -- (optional) a valuation of the base ring of this
          morphism

        OUTPUT: none

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: H = Hom(P, P)
            sage: f = H([5/4*x^3, 5*x*y^2])
            sage: f.normalize_coordinates(); f
            Scheme endomorphism of Projective Space of dimension 1 over Rational Field
              Defn: Defined on coordinates by sending (x : y) to (x^2 : 4*y^2)

        ::

            sage: P.<x,y,z> = ProjectiveSpace(GF(7), 2)
            sage: X = P.subscheme(x^2 - y^2)
            sage: H = Hom(X, X)
            sage: f = H([x^3 + x*y^2, x*y^2, x*z^2])
            sage: f.normalize_coordinates(); f                                          # needs sage.libs.singular
            Scheme endomorphism of Closed subscheme of Projective Space of dimension 2
             over Finite Field of size 7 defined by: x^2 - y^2
              Defn: Defined on coordinates by sending (x : y : z) to (2*y^2 : y^2 : z^2)

        ::

            sage: R.<a,b> = QQ[]
            sage: P.<x,y,z> = ProjectiveSpace(R, 2)
            sage: H = End(P)
            sage: f = H([a*(x*z + y^2)*x^2, a*b*(x*z + y^2)*y^2, a*(x*z + y^2)*z^2])
            sage: f.normalize_coordinates(); f
            Scheme endomorphism of Projective Space of dimension 2 over
             Multivariate Polynomial Ring in a, b over Rational Field
              Defn: Defined on coordinates by sending (x : y : z) to (x^2 : b*y^2 : z^2)

        ::

            sage: # needs sage.rings.number_field
            sage: K.<w> = QuadraticField(5)
            sage: P.<x,y> = ProjectiveSpace(K, 1)
            sage: f = DynamicalSystem([w*x^2 + (1/5*w)*y^2, w*y^2])
            sage: f.normalize_coordinates(); f
            Dynamical System of Projective Space of dimension 1 over Number Field in w
             with defining polynomial x^2 - 5 with w = 2.236067977499790?
              Defn: Defined on coordinates by sending (x : y) to (5*x^2 + y^2 : 5*y^2)

        ::

            sage: # needs sage.rings.number_field
            sage: R.<t> = PolynomialRing(ZZ)
            sage: K.<b> = NumberField(t^3 - 11)
            sage: a = 7/(b - 1)
            sage: P.<x,y> = ProjectiveSpace(K, 1)
            sage: f = DynamicalSystem_projective([a*y^2 - (a*y - x)^2, y^2])
            sage: f.normalize_coordinates(); f
            Dynamical System of Projective Space of dimension 1 over
             Number Field in b with defining polynomial t^3 - 11
              Defn: Defined on coordinates by sending (x : y) to
                    (-100*x^2 + (140*b^2 + 140*b + 140)*x*y + (-77*b^2 - 567*b - 1057)*y^2
                     : 100*y^2)

        We can used ``ideal`` to scale with respect to a norm defined by an ideal::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSystem_projective([2*x^3, 2*x^2*y + 4*x*y^2])
            sage: f.normalize_coordinates(ideal=2); f
            Dynamical System of Projective Space of dimension 1 over Rational Field
              Defn: Defined on coordinates by sending (x : y) to (x^3 : x^2*y + 2*x*y^2)

        ::

            sage: # needs sage.rings.number_field
            sage: R.<w> = QQ[]
            sage: A.<a> = NumberField(w^2 + 1)
            sage: P.<x,y,z> = ProjectiveSpace(A, 2)
            sage: X = P.subscheme(x^2 - y^2)
            sage: H = Hom(X, X)
            sage: f = H([(a+1)*x^3 + 2*x*y^2, 4*x*y^2, 8*x*z^2])
            sage: f.normalize_coordinates(ideal=A.prime_above(2)); f
            Scheme endomorphism of Closed subscheme of Projective Space of dimension 2 over
             Number Field in a with defining polynomial w^2 + 1 defined by: x^2 - y^2
              Defn: Defined on coordinates by sending (x : y : z) to
                    ((-a + 2)*x*y^2 : (-2*a + 2)*x*y^2 : (-4*a + 4)*x*z^2)

        We can pass in a valuation to ``valuation``::

            sage: g = H([(a+1)*x^3 + 2*x*y^2, 4*x*y^2, 8*x*z^2])                        # needs sage.rings.number_field
            sage: g.normalize_coordinates(valuation=A.valuation(A.prime_above(2)))      # needs sage.rings.number_field
            sage: g == f                                                                # needs sage.rings.number_field
            True

        ::

            sage: P.<x,y> = ProjectiveSpace(Qp(3), 1)                                   # needs sage.rings.padics
            sage: f = DynamicalSystem_projective([3*x^2 + 6*y^2, 9*x*y])                # needs sage.rings.padics
            sage: f.normalize_coordinates(); f                                          # needs sage.rings.padics
            Dynamical System of Projective Space of dimension 1 over
             3-adic Field with capped relative precision 20
              Defn: Defined on coordinates by sending (x : y) to
                    (x^2 + (2 + O(3^20))*y^2 : (3 + O(3^21))*x*y)

        Check that #35797 is fixed::

            sage: # needs sage.rings.number_field
            sage: R.<x> = QQ[]
            sage: K.<a> = NumberField(3*x^2 + 1)
            sage: P.<z,w> = ProjectiveSpace(K, 1)
            sage: f = DynamicalSystem_projective([a*(z^2 + w^2), z*w])
            sage: f.normalize_coordinates(); f
            Dynamical System of Projective Space of dimension 1 over
            Number Field in a with defining polynomial 3*x^2 + 1
            Defn: Defined on coordinates by sending (z : w) to
                ((3/2*a + 1/2)*z^2 + (3/2*a + 1/2)*w^2 : (-3/2*a + 3/2)*z*w)

        ::

            sage: R.<a,b> = QQ[]
            sage: P.<x,y,z> = ProjectiveSpace(FractionField(R), 2)
            sage: H = End(P)
            sage: f = H([a/b*(x*z + y^2)*x^2, a*b*(x*z + y^2)*y^2, a*(x*z + y^2)*z^2])
            sage: f.normalize_coordinates(); f
            Scheme endomorphism of Projective Space of dimension 2 over Fraction
            Field of Multivariate Polynomial Ring in a, b over Rational Field
            Defn: Defined on coordinates by sending (x : y : z) to
                (x^2 : (b^2)*y^2 : b*z^2)
        """
        # If ideal or valuation is specified, we scale according the norm
        # defined by the ideal/valuation
        ideal = kwds.pop('ideal', None)
        if ideal is not None:
            if not (ideal in ZZ or isinstance(ideal, NumberFieldFractionalIdeal)):
                raise TypeError('ideal must be an ideal of a number field, not %s' % ideal)
            if isinstance(ideal, NumberFieldFractionalIdeal):
                if ideal.number_field() != self.base_ring():
                    raise ValueError('ideal must be an ideal of the base ring of this morphism ' +
                                     ', not an ideal of %s' % ideal.number_field())
                if not ideal.is_prime():
                    raise ValueError('ideal was %s, not a prime ideal' % ideal)
                for generator in ideal.gens():
                    if generator.valuation(ideal) == 1:
                        uniformizer = generator
                        break
            else:
                ideal = ZZ(ideal)
                if self.base_ring() != QQ:
                    raise ValueError('ideal was an integer, but the base ring of this ' +
                        'morphism is %s' % self.base_ring())
                if not ideal.is_prime():
                    raise ValueError('ideal must be a prime, not %s' % ideal)
                uniformizer = ideal
            valuations = []
            for poly in self:
                for coefficient, monomial in poly:
                    if coefficient != 0:
                        valuations.append(coefficient.valuation(ideal))
            min_val = min(valuations)
            self.scale_by(uniformizer**(-1 * min_val))
            return

        valuation = kwds.pop('valuation', None)
        if valuation is not None:
            if not isinstance(valuation, pAdicValuation_base):
                raise TypeError('valuation must be a valuation on a number field, not %s' % valuation)
            if valuation.domain() != self.base_ring():
                raise ValueError('the domain of valuation must be the base ring of this morphism ' +
                    'not %s' % valuation.domain())
            uniformizer = valuation.uniformizer()
            ramification_index = 1 / valuation(uniformizer)
            valuations = []
            for poly in self:
                for coefficient, monomial in poly:
                    if coefficient != 0:
                        valuations.append(valuation(coefficient) * ramification_index)
            min_val = min(valuations)
            self.scale_by(uniformizer**(-1 * min_val))
            return

        N = self.codomain().ambient_space().dimension_relative() + 1

        R = self.domain().base_ring()

        # Clear any denominators from the coefficients
        if R in NumberFields():
            if R.maximal_order() is ZZ:
                denom = lcm([self[i].denominator() for i in range(N)])
            else:
                denom = R.ideal([c for poly in self for c in poly.coefficients()]).norm().denominator()

            self.scale_by(denom)
        else:
            self.scale_by(lcm([self[i].denominator() for i in range(N)]))

        # There are cases, such as the example above over GF(7),
        # where we want to compute GCDs, but NOT in the case
        # where R is a NumberField of class number > 1.
        if R in NumberFields():
            if R.class_number() > 1:
                return

        # R is a Number Field with class number 1 (i.e., a UFD) then
        # we can compute GCDs, so we attempt to remove any common factors.

        GCD = gcd(self[0], self[1])
        index = 2

        while GCD != 1 and index < N:
            GCD = gcd(GCD, self[index])
            index += +1
        if GCD != 1:
            self.scale_by(R(1) / GCD)

        # Scale by 1/GCD of the coefficients.
        if R in _NumberFields:
            O = R.maximal_order()
        elif isinstance(R, FiniteField):
            O = R
        elif isinstance(R, QuotientRing_generic):
            O = R.ring()
        elif isinstance(R, sage.rings.abc.pAdicField):
            O = R.integer_ring()
        else:
            O = R
        GCD = gcd([O(c) for poly in self for c in poly.coefficients()])

        if GCD != 1:
            self.scale_by(1 / GCD)

        # If R is not p-adic, we make the first coordinate positive
        if not isinstance(R, pAdicGeneric):
            if self[0].lc() < 0:
                self.scale_by(-1)

    def degree(self):
        r"""
        Return the degree of this map.

        The degree is defined as the degree of the homogeneous
        polynomials that are the coordinates of this map.

        OUTPUT: positive integer

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: H = Hom(P, P)
            sage: f = H([x^2 + y^2, y^2])
            sage: f.degree()
            2

        ::

            sage: # needs sage.rings.real_mpfr
            sage: P.<x,y,z> = ProjectiveSpace(CC, 2)
            sage: H = Hom(P, P)
            sage: f = H([x^3 + y^3, y^2*z, z*x*y])
            sage: f.degree()
            3

        ::

            sage: R.<t> = PolynomialRing(QQ)
            sage: P.<x,y,z> = ProjectiveSpace(R, 2)
            sage: H = Hom(P, P)
            sage: f = H([x^2 + t*y^2, (2-t)*y^2, z^2])
            sage: f.degree()
            2

        ::

            sage: P.<x,y,z> = ProjectiveSpace(ZZ, 2)
            sage: X = P.subscheme(x^2 - y^2)
            sage: H = Hom(X, X)
            sage: f = H([x^2, y^2, z^2])
            sage: f.degree()
            2
        """
        return self._polys[0].degree()

    def dehomogenize(self, n):
        r"""
        Return the standard dehomogenization at the ``n[0]`` coordinate for the domain
        and the ``n[1]`` coordinate for the codomain.

        Note that the new function is defined over the fraction field
        of the base ring of this map.

        INPUT:

        - ``n`` -- tuple of nonnegative integers; if ``n`` is an integer, then
          the two values of the tuple are assumed to be the same

        OUTPUT: :class:`SchemeMorphism_polynomial_affine_space`

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(ZZ, 1)
            sage: H = Hom(P, P)
            sage: f = H([x^2 + y^2, y^2])
            sage: f.dehomogenize(0)
            Scheme endomorphism of Affine Space of dimension 1 over Integer Ring
              Defn: Defined on coordinates by sending (y) to (y^2/(y^2 + 1))

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: H = Hom(P, P)
            sage: f = H([x^2 - y^2, y^2])
            sage: f.dehomogenize((0,1))
            Scheme morphism:
              From: Affine Space of dimension 1 over Rational Field
              To:   Affine Space of dimension 1 over Rational Field
              Defn: Defined on coordinates by sending (y) to ((-y^2 + 1)/y^2)

        ::

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: H = Hom(P, P)
            sage: f = H([x^2 + y^2, y^2 - z^2, 2*z^2])
            sage: f.dehomogenize(2)
            Scheme endomorphism of Affine Space of dimension 2 over Rational Field
              Defn: Defined on coordinates by sending (x, y) to
                    (1/2*x^2 + 1/2*y^2, 1/2*y^2 - 1/2)

        ::

            sage: R.<t> = PolynomialRing(QQ)
            sage: P.<x,y,z> = ProjectiveSpace(FractionField(R),2)
            sage: H = Hom(P,P)
            sage: f = H([x^2 + t*y^2, t*y^2 - z^2, t*z^2])
            sage: f.dehomogenize(2)
            Scheme endomorphism of Affine Space of dimension 2 over Fraction Field
             of Univariate Polynomial Ring in t over Rational Field
              Defn: Defined on coordinates by sending (x, y) to
                    (1/t*x^2 + y^2, y^2 - 1/t)

        ::

            sage: P.<x,y,z> = ProjectiveSpace(ZZ, 2)
            sage: X = P.subscheme(x^2 - y^2)
            sage: H = Hom(X, X)
            sage: f = H([x^2, y^2, x*z])
            sage: f.dehomogenize(2)                                                     # needs sage.libs.singular
            Scheme endomorphism of Closed subscheme of Affine Space of dimension 2
             over Integer Ring defined by: x^2 - y^2
              Defn: Defined on coordinates by sending (x, y) to (x, y^2/x)

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: H = End(P)
            sage: f = H([x^2 - 2*x*y, y^2])
            sage: f.dehomogenize(0).homogenize(0) == f
            True

        ::

            sage: # needs sage.rings.number_field
            sage: K.<w> = QuadraticField(3)
            sage: O = K.ring_of_integers()
            sage: P.<x,y> = ProjectiveSpace(O, 1)
            sage: H = End(P)
            sage: f = H([x^2 - O(w)*y^2, y^2])
            sage: f.dehomogenize(1)
            Scheme endomorphism of Affine Space of dimension 1 over
             Maximal Order generated by w in Number Field in w with defining polynomial x^2 - 3
              with w = 1.732050807568878?
              Defn: Defined on coordinates by sending (x) to (x^2 - w)

        ::

            sage: P1.<x,y> = ProjectiveSpace(QQ, 1)
            sage: P2.<u,v,w> = ProjectiveSpace(QQ, 2)
            sage: H = Hom(P2, P1)
            sage: f = H([u*w, v^2 + w^2])
            sage: f.dehomogenize((2,1))
            Scheme morphism:
              From: Affine Space of dimension 2 over Rational Field
              To:   Affine Space of dimension 1 over Rational Field
              Defn: Defined on coordinates by sending (u, v) to (u/(v^2 + 1))
        """
        # the dehomogenizations are stored for future use
        try:
            return self.__dehomogenization[n]
        except AttributeError:
            self.__dehomogenization = {}
        except KeyError:
            pass
        # it is possible to dehomogenize the domain and codomain at different coordinates
        if isinstance(n, (tuple, list)):
            ind = tuple(n)
        else:
            ind = (n, n)
        PS_domain = self.domain()
        A_domain = PS_domain.ambient_space()
        if self._polys[ind[1]].substitute({A_domain.gen(ind[0]): 1}) == 0:
            raise ValueError("can't dehomogenize at 0 coordinate")
        else:
            Aff_domain = PS_domain.affine_patch(ind[0])
            S = Aff_domain.ambient_space().coordinate_ring()
            FS = FractionField(S)
            N = A_domain.dimension_relative()
            R = A_domain.coordinate_ring()
            phi = R.hom([S.gen(j) for j in range(0, ind[0])] + [1] + [S.gen(j) for j in range(ind[0], N)], FS)
            F = []
            G = phi(self._polys[ind[1]])
            # ind[1] is relative to codomain
            M = self.codomain().ambient_space().dimension_relative()
            F.extend(phi(self._polys[i]) / G
                     for i in range(M + 1) if i != ind[1])
            H = Hom(Aff_domain, self.codomain().affine_patch(ind[1]))
            # since often you dehomogenize at the same coordinate in domain
            # and codomain it should be stored appropriately.
            if ind == (n, n):
                self.__dehomogenization[ind] = H(F)
                return self.__dehomogenization[ind]
            else:
                self.__dehomogenization[n] = H(F)
                return self.__dehomogenization[n]

    @cached_method
    def is_morphism(self):
        r"""
        Return ``True`` if this map is a morphism.

        The map is a morphism if and only if the ideal generated by
        the defining polynomials is the unit ideal
        (no common zeros of the defining polynomials).

        OUTPUT: boolean

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: H = Hom(P, P)
            sage: f = H([x^2 + y^2, y^2])
            sage: f.is_morphism()                                                       # needs sage.libs.singular
            True

        ::

            sage: P.<x,y,z> = ProjectiveSpace(RR, 2)
            sage: H = Hom(P, P)
            sage: f = H([x*z - y*z, x^2 - y^2, z^2])
            sage: f.is_morphism()                                                       # needs sage.libs.singular
            False

        ::

            sage: R.<t> = PolynomialRing(GF(5))
            sage: P.<x,y,z> = ProjectiveSpace(R, 2)
            sage: H = Hom(P, P)
            sage: f = H([x*z - t*y^2, x^2 - y^2, t*z^2])
            sage: f.is_morphism()                                                       # needs sage.libs.singular
            True

        Map that is not morphism on projective space, but is over a subscheme::

            sage: P.<x,y,z> = ProjectiveSpace(RR, 2)
            sage: X = P.subscheme([x*y + y*z])
            sage: H = Hom(X, X)
            sage: f = H([x*z - y*z, x^2 - y^2, z^2])
            sage: f.is_morphism()                                                       # needs sage.libs.singular
            True
        """

        R = self.coordinate_ring()
        F = list(self._polys)
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

    def global_height(self, prec=None):
        r"""
        Return the global height of the coefficients as a projective point.

        INPUT:

        - ``prec`` -- desired floating point precision (default:
          default RealField precision)

        OUTPUT: a real number

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: H = Hom(P, P)
            sage: f = H([1/1331*x^2 + 1/4000*y^2, 210*x*y]);
            sage: f.global_height()                                                     # needs sage.symbolic
            20.8348429892146

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: H = Hom(P, P)
            sage: f = H([1/1331*x^2 + 1/4000*y^2, 210*x*y]);
            sage: f.global_height(prec=11)                                              # needs sage.symbolic
            20.8

        ::

            sage: P.<x,y,z> = ProjectiveSpace(ZZ, 2)
            sage: H = Hom(P, P)
            sage: f = H([4*x^2 + 100*y^2, 210*x*y, 10000*z^2]);
            sage: f.global_height()                                                     # needs sage.symbolic
            8.51719319141624

        ::

            sage: # needs sage.rings.number_field
            sage: R.<z> = PolynomialRing(QQ)
            sage: K.<w> = NumberField(z^2 - 2)
            sage: O = K.maximal_order()
            sage: P.<x,y> = ProjectiveSpace(O, 1)
            sage: H = Hom(P, P)
            sage: f = H([2*x^2 + 3*O(w)*y^2, O(w)*y^2])
            sage: f.global_height()
            1.09861228866811

        ::

            sage: # needs sage.rings.number_field sage.symbolic
            sage: P.<x,y> = ProjectiveSpace(QQbar, 1)
            sage: P2.<u,v,w> = ProjectiveSpace(QQbar, 2)
            sage: H = Hom(P, P2)
            sage: f = H([x^2 + QQbar(I)*x*y + 3*y^2, y^2, QQbar(sqrt(5))*x*y])
            sage: f.global_height()
            1.09861228866811

        ::

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: A.<z,w> = ProjectiveSpace(QQ, 1)
            sage: H = Hom(P, A)
            sage: f = H([1/1331*x^2 + 4000*y*z, y^2])
            sage: f.global_height()                                                     # needs sage.symbolic
            15.4877354584971

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSystem([1/25*x^2 + 25/3*x*y + y^2, 1*y^2])
            sage: exp(f.global_height())                                                # needs sage.symbolic
            625.000000000000

        Scaling should not change the result::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSystem([1/25*x^2 + 25/3*x*y + y^2, 1*y^2])
            sage: f.global_height()                                                     # needs sage.symbolic
            6.43775164973640
            sage: c = 10000
            sage: f.scale_by(c)
            sage: f.global_height()                                                     # needs sage.symbolic
            6.43775164973640
        """
        K = self.domain().base_ring()
        if K in _NumberFields or K == ZZ or isinstance(K, sage.rings.abc.Order):
            f = self
        elif isinstance(K, sage.rings.abc.AlgebraicField):
            f = self._number_field_from_algebraics()
        else:
            raise TypeError("Must be over a Numberfield or a Numberfield Order or QQbar")

        # Get the coefficients from all of the polynomials in the dynamical system
        coeffs = [x for k in f for x in k.coefficients()]

        from sage.schemes.projective.projective_space import ProjectiveSpace

        P = ProjectiveSpace(K, len(coeffs) - 1)
        return P.point(coeffs).global_height(prec=prec)

    def local_height(self, v, prec=None):
        r"""
        Return the maximum of the local height of the coefficients in any
        of the coordinate functions of this map.

        INPUT:

        - ``v`` -- a prime or prime ideal of the base ring

        - ``prec`` -- desired floating point precision (default:
          default RealField precision)

        OUTPUT: a real number

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: H = Hom(P, P)
            sage: f = H([1/1331*x^2 + 1/4000*y^2, 210*x*y])
            sage: f.local_height(1331)                                                  # needs sage.rings.real_mpfr
            7.19368581839511

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: H = Hom(P, P)
            sage: f = H([1/1331*x^2 + 1/4000*y^2, 210*x*y])
            sage: f.local_height(1331, prec=2)                                          # needs sage.rings.real_mpfr
            8.0

        This function does not automatically normalize::

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: H = Hom(P, P)
            sage: f = H([4*x^2 + 3/100*y^2, 8/210*x*y, 1/10000*z^2])
            sage: f.local_height(2)                                                     # needs sage.rings.real_mpfr
            2.77258872223978
            sage: f.normalize_coordinates()                                             # needs sage.libs.singular
            sage: f.local_height(2)                                                     # needs sage.libs.singular
            0.000000000000000

        ::

            sage: # needs sage.rings.number_field
            sage: R.<z> = PolynomialRing(QQ)
            sage: K.<w> = NumberField(z^2 - 2)
            sage: P.<x,y> = ProjectiveSpace(K, 1)
            sage: H = Hom(P, P)
            sage: f = H([2*x^2 + w/3*y^2, 1/w*y^2])
            sage: f.local_height(K.ideal(3))
            1.09861228866811
        """
        K = FractionField(self.domain().base_ring())
        if K not in _NumberFields:
            raise TypeError("must be over a number field or a number field order")
        return max([K(c).local_height(v, prec=prec) for f in self for c in f.coefficients()])

    def local_height_arch(self, i, prec=None):
        r"""
        Return the maximum of the local height at the ``i``-th infinite place of the coefficients in any
        of the coordinate functions of this map.

        INPUT:

        - ``i`` -- integer

        - ``prec`` -- desired floating point precision (default:
          default RealField precision)

        OUTPUT: a real number

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: H = Hom(P, P)
            sage: f = H([1/1331*x^2 + 1/4000*y^2, 210*x*y])
            sage: f.local_height_arch(0)                                                # needs sage.rings.real_mpfr
            5.34710753071747

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: H = Hom(P, P)
            sage: f = H([1/1331*x^2 + 1/4000*y^2, 210*x*y])
            sage: f.local_height_arch(0, prec=5)                                        # needs sage.rings.real_mpfr
            5.2

        ::

            sage: # needs sage.rings.number_field
            sage: R.<z> = PolynomialRing(QQ)
            sage: K.<w> = NumberField(z^2 - 2)
            sage: P.<x,y> = ProjectiveSpace(K, 1)
            sage: H = Hom(P, P)
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

    def wronskian_ideal(self):
        r"""
        Return the ideal generated by the critical point locus.

        This is the vanishing of the maximal minors of the Jacobian matrix.
        Not implemented for subvarieties.

        OUTPUT: an ideal in the coordinate ring of the domain of this map

        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: R.<x> = PolynomialRing(QQ)
            sage: K.<w> = NumberField(x^2 + 11)
            sage: P.<x,y> = ProjectiveSpace(K, 1)
            sage: H = End(P)
            sage: f = H([x^2 - w*y^2, w*y^2])
            sage: f.wronskian_ideal()
            Ideal ((4*w)*x*y) of Multivariate Polynomial Ring in x, y
             over Number Field in w with defining polynomial x^2 + 11

        ::

            sage: # needs sage.rings.number_field
            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: P2.<u,v,t> = ProjectiveSpace(K, 2)
            sage: H = Hom(P, P2)
            sage: f = H([x^2 - 2*y^2, y^2, x*y])
            sage: f.wronskian_ideal()
            Ideal (4*x*y, 2*x^2 + 4*y^2, -2*y^2) of
             Multivariate Polynomial Ring in x, y over Rational Field
        """
        from sage.calculus.functions import jacobian

        dom = self.domain()
        from sage.schemes.projective.projective_space import ProjectiveSpace_ring
        if not (isinstance(dom, ProjectiveSpace_ring) and isinstance(self.codomain(), ProjectiveSpace_ring)):
            raise NotImplementedError("not implemented for subschemes")
        N = dom.dimension_relative() + 1
        R = dom.coordinate_ring()
        J = jacobian(self.defining_polynomials(), dom.gens())
        return R.ideal(J.minors(N))


class SchemeMorphism_polynomial_projective_space_field(SchemeMorphism_polynomial_projective_space):

    def rational_preimages(self, Q, k=1):
        r"""
        Determine all of the rational `k`-th preimages of ``Q`` by this map.

        Given a rational point ``Q`` in the domain of this map, return all the rational points ``P``
        in the domain with `f^k(P)==Q`. In other words, the set of `k`-th preimages of ``Q``.
        The map must be defined over a number field and be an endomorphism for `k > 1`.

        If ``Q`` is a subscheme, then return the subscheme that maps to ``Q`` by this map.
        In particular, `f^{-k}(V(h_1,\ldots,h_t)) = V(h_1 \circ f^k, \ldots, h_t \circ f^k)`.

        INPUT:

        - ``Q`` -- a rational point or subscheme in the domain of this map

        - ``k`` -- positive integer

        OUTPUT: a list of rational points or a subscheme in the domain of this map

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: H = End(P)
            sage: f = H([16*x^2 - 29*y^2, 16*y^2])
            sage: f.rational_preimages(P(-1, 4))                                        # needs sage.libs.singular
            [(-5/4 : 1), (5/4 : 1)]

        ::

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: H = End(P)
            sage: f = H([76*x^2 - 180*x*y + 45*y^2 + 14*x*z + 45*y*z - 90*z^2,
            ....:        67*x^2 - 180*x*y - 157*x*z + 90*y*z,
            ....:        -90*z^2])
            sage: f.rational_preimages(P(-9, -4, 1))                                    # needs sage.libs.singular
            [(0 : 4 : 1)]

        A non-periodic example ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: H = End(P)
            sage: f = H([x^2 + y^2, 2*x*y])
            sage: f.rational_preimages(P(17, 15))                                       # needs sage.libs.singular
            [(3/5 : 1), (5/3 : 1)]

        ::

            sage: P.<x,y,z,w> = ProjectiveSpace(QQ, 3)
            sage: H = End(P)
            sage: f = H([x^2 - 2*y*w - 3*w^2, -2*x^2 + y^2 - 2*x*z + 4*y*w + 3*w^2,
            ....:        x^2 - y^2 + 2*x*z + z^2 - 2*y*w - w^2,
            ....:        w^2])
            sage: f.rational_preimages(P(0, -1, 0, 1))                                  # needs sage.libs.singular
            []

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: H = End(P)
            sage: f = H([x^2 + y^2, 2*x*y])
            sage: f.rational_preimages([CC.0, 1])                                       # needs sage.libs.singular
            Traceback (most recent call last):
            ...
            TypeError: point must be in codomain of self

        A number field example ::

            sage: # needs sage.rings.number_field
            sage: z = QQ['z'].0
            sage: K.<a> = NumberField(z^2 - 2)
            sage: P.<x,y> = ProjectiveSpace(K, 1)
            sage: H = End(P)
            sage: f = H([x^2 + y^2, y^2])
            sage: f.rational_preimages(P(3, 1))                                         # needs sage.libs.singular
            [(-a : 1), (a : 1)]

        ::

            sage: # needs sage.rings.number_field
            sage: z = QQ['z'].0
            sage: K.<a> = NumberField(z^2 - 2)
            sage: P.<x,y,z> = ProjectiveSpace(K, 2)
            sage: X = P.subscheme([x^2 - z^2])
            sage: H = End(X)
            sage: f= H([x^2 - z^2, a*y^2, z^2 - x^2])
            sage: f.rational_preimages(X([1, 2, -1]))                                   # needs sage.libs.singular
            []

        ::

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: X = P.subscheme([x^2 - z^2])
            sage: H = End(X)
            sage: f = H([x^2-z^2, y^2, z^2-x^2])
            sage: f.rational_preimages(X([0, 1, 0]))                                    # needs sage.libs.singular
            Closed subscheme of Projective Space of dimension 2 over Rational Field defined by:
              x^2 - z^2,
              -x^2 + z^2,
              0,
              -x^2 + z^2

        ::

            sage: P.<x, y> = ProjectiveSpace(QQ, 1)
            sage: H = End(P)
            sage: f = H([x^2 - y^2, y^2])
            sage: f.rational_preimages(P.subscheme([x]))                                # needs sage.libs.singular
            Closed subscheme of Projective Space of dimension 1 over Rational Field
             defined by: x^2 - y^2

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: H = End(P)
            sage: f = H([x^2 - 29/16*y^2, y^2])
            sage: f.rational_preimages(P(5/4, 1), k=4)                                  # needs sage.libs.singular
            [(-3/4 : 1), (3/4 : 1), (-7/4 : 1), (7/4 : 1)]

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: P2.<u,v,w> = ProjectiveSpace(QQ, 2)
            sage: H = Hom(P, P2)
            sage: f = H([x^2, y^2, x^2-y^2])
            sage: f.rational_preimages(P2(1, 1, 0))                                     # needs sage.libs.singular
            [(-1 : 1), (1 : 1)]
        """
        k = ZZ(k)
        if k <= 0:
            raise ValueError("k (=%s) must be a positive integer" % k)
        # first check if subscheme
        from sage.schemes.projective.projective_subscheme import AlgebraicScheme_subscheme_projective
        if isinstance(Q, AlgebraicScheme_subscheme_projective):
            return Q.preimage(self, k)

        # else assume a point
        BR = self.base_ring()
        if k > 1 and not self.is_endomorphism():
            raise TypeError("must be an endomorphism of projective space")
        if Q not in self.codomain():
            raise TypeError("point must be in codomain of self")
        if isinstance(BR.base_ring(), (sage.rings.abc.ComplexField, sage.rings.abc.RealField,
                                       sage.rings.abc.RealIntervalField, sage.rings.abc.ComplexIntervalField)):
            raise NotImplementedError("not implemented over precision fields")
        PS = self.domain().ambient_space()
        N = PS.dimension_relative()
        L = [Q]
        for n in range(k):
            L2 = []
            for P in L:
                I = list(self.domain().defining_polynomials())
                I.extend(P[i] * self[j] - P[j] * self[i]
                         for i in range(N + 1) for j in range(i + 1, N + 1))
                X = PS.subscheme(I)
                if X.dimension() > 0:
                    return X
                L2.extend(PS(T) for T in X.rational_points()
                          if not all(g(tuple(T)) == 0 for g in self))
            L = L2
        return L

    def _number_field_from_algebraics(self):
        r"""
        Given a projective map defined over `\QQbar`, return the same map, but defined
        over a number field.

        This is only implemented for maps of projective space.

        OUTPUT: scheme morphism

        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: R.<x> = PolynomialRing(QQ)
            sage: P.<x,y> = ProjectiveSpace(QQbar, 1)
            sage: H = End(P)
            sage: f = H([QQbar(3^(1/3))*x^2 + QQbar(sqrt(-2))*y^2, y^2])                # needs sage.symbolic
            sage: f._number_field_from_algebraics()                                     # needs sage.symbolic
            Scheme endomorphism of Projective Space of dimension 1 over Number
             Field in a with defining polynomial y^6 + 6*y^4 - 6*y^3 + 12*y^2 + 36*y + 17
             with a = 1.442249570307409? - 1.414213562373095?*I
              Defn: Defined on coordinates by sending (x : y) to
                    ((-48/269*a^5 + 27/269*a^4 - 320/269*a^3 + 468/269*a^2 - 772/269*a
                    - 1092/269)*x^2 + (-48/269*a^5 + 27/269*a^4 - 320/269*a^3 + 468/269*a^2
                    - 1041/269*a - 1092/269)*y^2 : y^2)

        ::

            sage: # needs sage.rings.number_field
            sage: P.<x,y> = ProjectiveSpace(QQbar, 1)
            sage: P2.<u,v,w> = ProjectiveSpace(QQbar, 2)
            sage: H = Hom(P, P2)
            sage: f = H([x^2 + QQbar(I)*x*y + 3*y^2, y^2, QQbar(sqrt(5))*x*y])          # needs sage.symbolic
            sage: f._number_field_from_algebraics()                                     # needs sage.symbolic
            Scheme morphism:
              From: Projective Space of dimension 1 over Number Field in a
                    with defining polynomial y^4 + 3*y^2 + 1
                    with a = 0.?e-166 + 1.618033988749895?*I
              To:   Projective Space of dimension 2 over Number Field in a
                    with defining polynomial y^4 + 3*y^2 + 1
                    with a = 0.?e-166 + 1.618033988749895?*I
              Defn: Defined on coordinates by sending (x : y) to
                    (x^2 + (-a^3 - 2*a)*x*y + 3*y^2 : y^2 : (-2*a^2 - 3)*x*y)

        The following was fixed in :issue:`23808`::

            sage: # needs sage.rings.number_field
            sage: R.<t> = PolynomialRing(QQ)
            sage: s = (t^3 + t + 1).roots(QQbar)[0][0]
            sage: P.<x,y> = ProjectiveSpace(QQbar, 1)
            sage: H = Hom(P, P)
            sage: f = H([s*x^3 - 13*y^3, y^3 - 15*y^3])
            sage: f
            Scheme endomorphism of Projective Space of dimension 1 over Algebraic Field
              Defn: Defined on coordinates by sending (x : y) to
                    ((-0.6823278038280193?)*x^3 + (-13)*y^3 : (-14)*y^3)
            sage: f_alg = f._number_field_from_algebraics()
            sage: f_alg.change_ring(QQbar)  # Used to fail
            Scheme endomorphism of Projective Space of dimension 1 over Algebraic Field
              Defn: Defined on coordinates by sending (x : y) to
                    ((-0.6823278038280193?)*x^3 + (-13)*y^3 : (-14)*y^3)
        """
        from sage.rings.qqbar import number_field_elements_from_algebraics
        from sage.schemes.projective.projective_space import ProjectiveSpace_ring

        if not (isinstance(self.domain(), ProjectiveSpace_ring) and isinstance(self.domain(), ProjectiveSpace_ring)):
            raise NotImplementedError("not implemented for subschemes")

        K_pre, C, phi = number_field_elements_from_algebraics([c for f in self
            for c in f.coefficients()], minimal=True)
        # check if the same field
        if K_pre is QQ:
            if K_pre is self.base_ring():
                return self
        elif not isinstance(self.base_ring(), sage.rings.abc.AlgebraicField) and K_pre.is_isomorphic(self.base_ring()):
            return self
        # Issue 23808: The field K_pre returned above does not have its embedding set to be phi
        # and phi is forgotten, so we redefine K_pre to be a field K with phi as the specified
        # embedding:
        if K_pre is QQ:
            K = QQ
        else:
            from sage.rings.number_field.number_field import NumberField
            K = NumberField(K_pre.polynomial(), embedding=phi(K_pre.gen()), name='a')
            psi = K_pre.hom([K.gen()], K)  # Identification of K_pre with K
            C = [psi(c) for c in C]  # The elements of C were in K_pre, move them to K
        from sage.schemes.projective.projective_space import ProjectiveSpace
        N = self.domain().dimension_relative()
        PS = ProjectiveSpace(K, N, self.domain().variable_names())
        if self.is_endomorphism():
            H = End(PS)
        else:
            PS2 = ProjectiveSpace(K, self.codomain().dimension_relative(),
                                  self.codomain().variable_names())
            H = Hom(PS, PS2)
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
        return H(F)

    def base_indeterminacy_locus(self):
        r"""
        Return the base indeterminacy locus of this map.

        The base indeterminacy locus is the set of points in projective space
        at which all of the defining polynomials of the rational map
        simultaneously vanish.

        OUTPUT: a subscheme of the domain of the map

        EXAMPLES::

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: H = End(P)
            sage: f = H([x*z - y*z, x^2 - y^2, z^2])
            sage: f.base_indeterminacy_locus()
            Closed subscheme of Projective Space of dimension 2 over Rational Field defined by:
              x*z - y*z,
              x^2 - y^2,
              z^2

        ::

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: H = End(P)
            sage: f = H([x^2, y^2, z^2])
            sage: f.base_indeterminacy_locus()
            Closed subscheme of Projective Space of dimension 2 over Rational Field
             defined by:
              x^2,
              y^2,
              z^2

        ::

            sage: P1.<x,y,z> = ProjectiveSpace(RR, 2)
            sage: P2.<t,u,v,w> = ProjectiveSpace(RR, 3)
            sage: H = Hom(P1, P2)
            sage: h = H([y^3*z^3, x^3*z^3, y^3*z^3, x^2*y^2*z^2])
            sage: h.base_indeterminacy_locus()                                          # needs sage.rings.real_mpfr
            Closed subscheme of Projective Space of dimension 2 over
             Real Field with 53 bits of precision defined by:
              y^3*z^3,
              x^3*z^3,
              y^3*z^3,
              x^2*y^2*z^2

        If defining polynomials are not normalized, output scheme will not be normalized::

            sage: P.<x,y,z> = ProjectiveSpace(QQ,2)
            sage: H = End(P)
            sage: f = H([x*x^2,x*y^2,x*z^2])
            sage: f.base_indeterminacy_locus()
            Closed subscheme of Projective Space of dimension 2 over Rational Field
             defined by:
              x^3,
              x*y^2,
              x*z^2
        """
        dom = self.domain()
        AS = dom.ambient_space()
        return AS.subscheme(list(dom.defining_polynomials()) + list(self.defining_polynomials()))

    def indeterminacy_locus(self):
        r"""
        Return the indeterminacy locus of this map as a rational map on the domain.

        The indeterminacy locus is the intersection of all the base indeterminacy
        locuses of maps that define the same rational map as by this map.

        OUTPUT: a subscheme of the domain of the map

        EXAMPLES::

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: H = End(P)
            sage: f = H([x^2, y^2, z^2])
            sage: f.indeterminacy_locus()                                               # needs sage.libs.singular
            ... DeprecationWarning: The meaning of indeterminacy_locus() has changed.
            Read the docstring. See https://github.com/sagemath/sage/issues/29145 for details.
            Closed subscheme of Projective Space of dimension 2 over Rational Field defined by:
              z,
              y,
              x

        ::

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: H = End(P)
            sage: f = H([x*z - y*z, x^2 - y^2, z^2])
            sage: f.indeterminacy_locus()                                               # needs sage.libs.singular
            Closed subscheme of Projective Space of dimension 2 over Rational Field defined by:
              z,
              x^2 - y^2

        There is related :meth:`base_indeterminacy_locus()` method. This
        computes the indeterminacy locus only from the defining polynomials of
        the map::

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: H = End(P)
            sage: f = H([x*z - y*z, x^2 - y^2, z^2])
            sage: f.base_indeterminacy_locus()
            Closed subscheme of Projective Space of dimension 2 over Rational Field defined by:
              x*z - y*z,
              x^2 - y^2,
              z^2
        """
        from sage.misc.superseded import deprecation
        deprecation(29145, "The meaning of indeterminacy_locus() has changed. Read the docstring.")
        P = self.domain()
        X = P.subscheme(0)  # projective space as a subscheme
        return (self*X.hom(P.gens(), P)).indeterminacy_locus()

    def indeterminacy_points(self, F=None, base=False):
        r"""
        Return the points in the indeterminacy locus of this map.

        If the dimension of the indeterminacy locus is not zero, an error is raised.

        INPUT:

        - ``F`` -- a field; if not given, the base ring of the domain is assumed

        - ``base`` -- if ``True``, the base indeterminacy locus is used

        OUTPUT: indeterminacy points of the map defined over ``F``

        EXAMPLES::

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: H = End(P)
            sage: f = H([x*z - y*z, x^2 - y^2, z^2])
            sage: f.indeterminacy_points()                                              # needs sage.libs.singular
            ... DeprecationWarning: The meaning of indeterminacy_locus() has changed.
            Read the docstring. See https://github.com/sagemath/sage/issues/29145 for details.
            [(-1 : 1 : 0), (1 : 1 : 0)]

        ::

            sage: P1.<x,y,z> = ProjectiveSpace(RR, 2)
            sage: P2.<t,u,v,w> = ProjectiveSpace(RR, 3)
            sage: H = Hom(P1, P2)
            sage: h = H([x + y, y, z + y, y])
            sage: set_verbose(None)
            sage: h.indeterminacy_points(base=True)                                     # needs sage.libs.singular
            []
            sage: g = H([y^3*z^3, x^3*z^3, y^3*z^3, x^2*y^2*z^2])
            sage: g.indeterminacy_points(base=True)                                     # needs sage.libs.singular
            Traceback (most recent call last):
            ...
            ValueError: indeterminacy scheme is not dimension 0

        ::

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: H = End(P)
            sage: f = H([x^2 + y^2, x*z, x^2 + y^2])
            sage: f.indeterminacy_points()                                              # needs sage.libs.singular
            [(0 : 0 : 1)]

            sage: R.<t> = QQ[]
            sage: K.<a> = NumberField(t^2 + 1)                                          # needs sage.rings.number_field
            sage: f.indeterminacy_points(F=K)                                           # needs sage.libs.singular sage.rings.number_field
            [(-a : 1 : 0), (0 : 0 : 1), (a : 1 : 0)]
            sage: set_verbose(None)
            sage: f.indeterminacy_points(F=QQbar, base=True)                            # needs sage.libs.singular sage.rings.number_field
            [(-1*I : 1 : 0), (0 : 0 : 1), (1*I : 1 : 0)]

        ::

            sage: set_verbose(None)
            sage: K.<t> = FunctionField(QQ)
            sage: P.<x,y,z> = ProjectiveSpace(K, 2)
            sage: H = End(P)
            sage: f = H([x^2 - t^2*y^2, y^2 - z^2, x^2 - t^2*z^2])
            sage: f.indeterminacy_points(base=True)                                     # needs sage.libs.singular
            [(-t : -1 : 1), (-t : 1 : 1), (t : -1 : 1), (t : 1 : 1)]

        ::

            sage: # needs sage.rings.padics
            sage: set_verbose(None)
            sage: P.<x,y,z> = ProjectiveSpace(Qp(3), 2)
            sage: H = End(P)
            sage: f = H([x^2 - 7*y^2, y^2 - z^2, x^2 - 7*z^2])
            sage: f.indeterminacy_points(base=True)                                     # needs sage.libs.singular
            [(2 + 3 + 3^2 + 2*3^3 + 2*3^5 + 2*3^6 + 3^8
                + 3^9 + 2*3^11 + 3^15 + 2*3^16 + 3^18 + O(3^20) : 1 + O(3^20) : 1 + O(3^20)),
             (2 + 3 + 3^2 + 2*3^3 + 2*3^5 + 2*3^6 + 3^8 + 3^9 + 2*3^11 + 3^15
                + 2*3^16 + 3^18 + O(3^20) : 2 + 2*3 + 2*3^2 + 2*3^3 + 2*3^4 + 2*3^5
                + 2*3^6 + 2*3^7 + 2*3^8 + 2*3^9 + 2*3^10 + 2*3^11 + 2*3^12 + 2*3^13
                + 2*3^14 + 2*3^15 + 2*3^16 + 2*3^17 + 2*3^18 + 2*3^19 + O(3^20) : 1 + O(3^20)),
             (1 + 3 + 3^2 + 2*3^4 + 2*3^7 + 3^8 + 3^9 + 2*3^10 + 2*3^12 + 2*3^13
                + 2*3^14 + 3^15 + 2*3^17 + 3^18 + 2*3^19 + O(3^20) : 1 + O(3^20) : 1 + O(3^20)),
             (1 + 3 + 3^2 + 2*3^4 + 2*3^7 + 3^8 + 3^9 + 2*3^10 + 2*3^12 + 2*3^13
                + 2*3^14 + 3^15 + 2*3^17 + 3^18 + 2*3^19 + O(3^20) : 2 + 2*3 + 2*3^2
                + 2*3^3 + 2*3^4 + 2*3^5 + 2*3^6 + 2*3^7 + 2*3^8 + 2*3^9 + 2*3^10 + 2*3^11
                + 2*3^12 + 2*3^13 + 2*3^14 + 2*3^15 + 2*3^16 + 2*3^17 + 2*3^18 + 2*3^19
                + O(3^20) : 1 + O(3^20))]
        """
        if F is None:
            fcn = self
        else:
            if not F.is_field():
                raise NotImplementedError("indeterminacy points only implemented for fields")
            fcn = self.change_ring(F)
        if base:
            indScheme = fcn.base_indeterminacy_locus()
        else:
            indScheme = fcn.indeterminacy_locus()
        if indScheme.dimension() > 0:
            raise ValueError("indeterminacy scheme is not dimension 0")
        # Other error checking is in indeterminacy_locus
        indPoints = indScheme.rational_points()
        return indPoints

    def reduce_base_field(self):
        """
        Return this map defined over the field of definition of the coefficients.

        The base field of the map could be strictly larger than
        the field where all of the coefficients are defined. This function
        reduces the base field to the minimal possible. This can be done when
        the base ring is a number field, QQbar, a finite field, or algebraic
        closure of a finite field.

        OUTPUT: a scheme morphism

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<t> = GF(3^4)
            sage: P.<x,y> = ProjectiveSpace(K, 1)
            sage: P2.<a,b,c> = ProjectiveSpace(K, 2)
            sage: H = End(P)
            sage: H2 = Hom(P, P2)
            sage: H3 = Hom(P2, P)
            sage: f = H([x^2 + (2*t^3 + 2*t^2 + 1)*y^2, y^2])
            sage: f.reduce_base_field()                                                 # needs sage.libs.singular sage.modules
            Scheme endomorphism of Projective Space of dimension 1
             over Finite Field in t2 of size 3^2
              Defn: Defined on coordinates by sending (x : y) to (x^2 + t2*y^2 : y^2)
            sage: f2 = H2([x^2 + 5*y^2, y^2, 2*x*y])
            sage: f2.reduce_base_field()                                                # needs sage.libs.singular sage.modules
            Scheme morphism:
              From: Projective Space of dimension 1 over Finite Field of size 3
              To:   Projective Space of dimension 2 over Finite Field of size 3
              Defn: Defined on coordinates by sending (x : y) to (x^2 - y^2 : y^2 : -x*y)
            sage: f3 = H3([a^2 + t*b^2, c^2])
            sage: f3.reduce_base_field()                                                # needs sage.libs.singular sage.modules
            Scheme morphism:
              From: Projective Space of dimension 2 over Finite Field in t of size 3^4
              To:   Projective Space of dimension 1 over Finite Field in t of size 3^4
              Defn: Defined on coordinates by sending (a : b : c) to (a^2 + t*b^2 : c^2)

        ::

            sage: # needs sage.rings.number_field
            sage: K.<v> = CyclotomicField(4)
            sage: P.<x,y> = ProjectiveSpace(K, 1)
            sage: H = End(P)
            sage: f = H([x^2 + 2*y^2, y^2])
            sage: f.reduce_base_field()                                                 # needs sage.libs.singular
            Scheme endomorphism of Projective Space of dimension 1 over Rational Field
              Defn: Defined on coordinates by sending (x : y) to (x^2 + 2*y^2 : y^2)

        ::

            sage: # needs sage.rings.finite_rings
            sage: K.<v> = GF(5)
            sage: L = K.algebraic_closure()
            sage: P.<x,y> = ProjectiveSpace(L, 1)
            sage: H = End(P)
            sage: f = H([(L.gen(2))*x^2 + L.gen(4)*y^2, x*y])
            sage: f.reduce_base_field()                                                 # needs sage.libs.singular
            Scheme endomorphism of Projective Space of dimension 1
             over Finite Field in z4 of size 5^4
              Defn: Defined on coordinates by sending (x : y) to
                    ((z4^3 + z4^2 + z4 - 2)*x^2 + z4*y^2 : x*y)
            sage: f = DynamicalSystem_projective([L.gen(3)*x^2 + L.gen(2)*y^2, x*y])    # needs sage.schemes
            sage: f.reduce_base_field()                                                 # needs sage.libs.singular sage.schemes
            Dynamical System of Projective Space of dimension 1
             over Finite Field in z6 of size 5^6
              Defn: Defined on coordinates by sending (x : y) to
                    ((-z6^5 + z6^4 - z6^3 - z6^2 - 2*z6 - 2)*x^2
                     + (z6^5 - 2*z6^4 + z6^2 - z6 + 1)*y^2 : x*y)

        TESTS::

            sage: # needs sage.rings.finite_rings
            sage: F = GF(3).algebraic_closure()
            sage: P.<x,y> = ProjectiveSpace(F, 1)
            sage: H = Hom(P, P)
            sage: f = H([x^2 + y^2, y^2])
            sage: f.reduce_base_field()                                                 # needs sage.libs.singular
            Scheme endomorphism of Projective Space of dimension 1 over Finite Field of size 3
              Defn: Defined on coordinates by sending (x : y) to
                    (x^2 + y^2 : y^2)
        """
        K = self.base_ring()
        if K in NumberFields() or isinstance(K, sage.rings.abc.AlgebraicField):
            return self._number_field_from_algebraics()
        if K in FiniteFields():
            #find the degree of the extension containing the coefficients
            c = [v for g in self for v in g.coefficients()]
            d = lcm([a.minpoly().degree() for a in c])
            if d == 1:
                from sage.rings.finite_rings.finite_field_constructor import GF

                return self.change_ring(GF(K.characteristic()))
            if d == K.degree():
                return self
            # otherwise we are not in the prime subfield so coercion
            # to it does not work
            for L, phi in K.subfields():
                # find the right subfield and its embedding
                if L.degree() == d:
                    break
            # we need to rewrite each of the coefficients in terms of the generator
            # of L. To do this, we'll set-up an ideal and use elimination
            R = PolynomialRing(K.prime_subfield(), 2, 'a')
            a, b = R.gens()
            from sage.schemes.projective.projective_space import ProjectiveSpace
            new_domain = ProjectiveSpace(L, self.domain().dimension_relative(),
                                         self.domain().variable_names())
            new_R = new_domain.coordinate_ring()
            u = phi(L.gen())  # gen of L in terms of gen of K
            g = R(str(u).replace(K.variable_name(), R.variable_names()[0])) #converted to R
            new_f = []
            for fi in self:
                mon = fi.monomials()
                mon_deg = [m.degrees() for m in mon]
                coef = fi.coefficients()
                new_c = []
                for c in coef:
                    # for each coefficient do the elimination
                    w = R(str(c).replace(K.variable_name(), R.variable_names()[0]))
                    I = R.ideal([b-g, w])
                    v = I.elimination_ideal([a]).gen(0)
                    # elimination can change scale the result, so correct the leading coefficient
                    # and convert back to L
                    if v.subs({b:g}).lc() == w.lc():
                        new_c.append(L(str(v).replace(R.variable_names()[1], L.variable_name())))
                    else:
                        new_c.append(L(str(w.lc()*v).replace(R.variable_names()[1], L.variable_name())))
                # reconstruct as a poly in the new domain
                new_f.append(sum(new_c[i]*prod(new_R.gen(j)**mon_deg[i][j]
                                               for j in range(new_R.ngens()))
                                 for i in range(len(mon))))
            # return the correct type of map
            if self.is_endomorphism():
                H = Hom(new_domain, new_domain)
            else:
                new_codomain = ProjectiveSpace(L, self.codomain().dimension_relative(), self.codomain().variable_names())
                H = Hom(new_domain, new_codomain)
            return H(new_f)
        elif isinstance(K, AlgebraicClosureFiniteField_generic):
            self.domain().coordinate_ring()
            # find the degree of the extension containing the coefficients
            c = [v for g in self for v in g.coefficients()]
            d = lcm([a.minpoly().degree() for a in c])
            # get the appropriate subfield
            L, L_to_K = K.subfield(d)
            from sage.schemes.projective.projective_space import ProjectiveSpace
            new_domain = ProjectiveSpace(L, self.domain().dimension_relative(),
                                         self.domain().variable_names())
            new_R = new_domain.coordinate_ring()
            # we need to rewrite each of the coefficients in terms of the generator
            # of L. To do this, we'll set-up an ideal and use elimination
            new_f = []
            for fi in self:
                mon = fi.monomials()
                mon_deg = [m.degrees() for m in mon]
                coef = fi.coefficients()
                new_c = []
                for c in coef:
                    # for each coefficient move to the correct base field
                    da = c.minpoly().degree()
                    for M, M_to_L in L.subfields():
                        # find the right subfield and it's embedding
                        if M.degree() == da:
                            break
                    c = M(str(c).replace(c.as_finite_field_element()[0].variable_name(),
                                          M.variable_name()))
                    new_c.append(M_to_L(c))
                # reconstruct as a poly in the new domain
                new_f.append(sum([new_c[i] * prod(new_R.gen(j)**mon_deg[i][j]
                                                  for j in range(new_R.ngens()))
                                  for i in range(len(mon))]))
            # return the correct type of map
            if self.is_endomorphism():
                H = Hom(new_domain, new_domain)
            else:
                new_codomain = ProjectiveSpace(L, self.codomain().dimension_relative(), self.codomain().variable_names())
                H = Hom(new_domain, new_codomain)
            return H(new_f)
        raise NotImplementedError("only implemented for number fields and finite fields")

    def image(self):
        """
        Return the scheme-theoretic image of the morphism.

        OUTPUT: a subscheme of the ambient space of the codomain

        EXAMPLES::

            sage: P2.<x0,x1,x2> = ProjectiveSpace(QQ, 2)
            sage: f = P2.hom([x0^3, x0^2*x1, x0*x1^2], P2)
            sage: f.image()                                                             # needs sage.libs.singular
            Closed subscheme of Projective Space of dimension 2 over Rational Field defined by:
              x1^2 - x0*x2
            sage: f = P2.hom([x0 - x1, x0 - x2, x1 - x2], P2)
            sage: f.image()                                                             # needs sage.libs.singular
            Closed subscheme of Projective Space of dimension 2 over Rational Field defined by:
              x0 - x1 + x2

        ::

            sage: P2.<x0,x1,x2> = ProjectiveSpace(QQ, 2)
            sage: A2.<x,y> = AffineSpace(QQ, 2)
            sage: f = P2.hom([1, x0/x1], A2)
            sage: f.image()                                                             # needs sage.libs.singular
            Closed subscheme of Affine Space of dimension 2 over Rational Field defined by:
              -x + 1
        """
        X = self.domain().subscheme(0)
        e = X.embedding_morphism()
        return (self*e).image()


class SchemeMorphism_polynomial_projective_space_finite_field(SchemeMorphism_polynomial_projective_space_field):

    def _fast_eval(self, x):
        """
        Evaluate projective morphism at point described by x.

        EXAMPLES::

            sage: P.<x,y,z> = ProjectiveSpace(GF(7), 2)
            sage: H = Hom(P, P)
            sage: f = H([x^2 + y^2, y^2, z^2 + y*z])
            sage: f._fast_eval([1,1,1])
            [2, 1, 2]
        """
        if self._is_prime_finite_field:
            p = self.base_ring().characteristic()
            P = [Integer(f(*x)) % p for f in self._fastpolys]
        else:
            P = [f(*x) for f in self._fastpolys]
        return P


class SchemeMorphism_polynomial_projective_subscheme_field(SchemeMorphism_polynomial_projective_space_field):
    """
    Morphisms from subschemes of projective spaces defined over fields.
    """
    def __call__(self, x):
        """
        Apply this morphism to the point ``x``.

        INPUT:

        - ``x`` -- a point in the domain of definition

        OUTPUT: the image of the point ``x`` under the morphism

        TESTS::

            sage: # needs sage.libs.pari sage.schemes
            sage: R.<x,y,z> = QQ[]
            sage: C = Curve(7*x^2 + 2*y*z + z^2)
            sage: f, g = C.parametrization()
            sage: g([0, -1, 2])
            (1 : 0)
            sage: f([1, 0])
            (0 : -1/2 : 1)
            sage: _ == C([0, -1, 2])
            True
        """
        try:
            reprs = self.representatives()
        except NotImplementedError:  # Singular does not support the base field
            try:
                return super().__call__(x)
            except ValueError:
                raise ValueError('cannot apply the morphism to this point')

        for m in reprs:
            try:
                return super(SchemeMorphism_polynomial_projective_subscheme_field, m).__call__(x)
            except ValueError:
                pass
        raise ValueError('the morphism is not defined at this point')

    def __eq__(self, other):
        """
        EXAMPLES::

            sage: R.<x,y,z> = QQ[]

            sage: # needs sage.libs.pari sage.schemes
            sage: C = Curve(7*x^2 + 2*y*z + z^2)  # conic
            sage: f, g = C.parametrization()
            sage: f*g == C.identity_morphism()
            True

            sage: # needs sage.schemes
            sage: C = Curve(x^2 + y^2 - z^2)
            sage: P.<u, v> = ProjectiveSpace(QQ, 1)
            sage: f = C.hom([x + z, y], P)
            sage: g = C.hom([y, z - x], P)
            sage: f == g
            True
            sage: h = C.hom([z, x - y], P)
            sage: f == h
            False
        """
        Y = self.codomain()

        if not isinstance(other, SchemeMorphism_polynomial):
            return False
        if self.domain() != other.domain() or Y != other.codomain():
            return False

        if not Y.is_projective():  # codomain is affine
            e = Y.projective_embedding(0)
            return (e * self) == (e * other)

        from sage.matrix.constructor import matrix

        R = self.domain().coordinate_ring()
        mat = matrix([self.defining_polynomials(), other.defining_polynomials()])
        return all(R(minor).is_zero() for minor in mat.minors(2))

    @cached_method
    def representatives(self):
        """
        Return all maps representing the same rational map as by this map.

        EXAMPLES::

            sage: P2.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: X = P2.subscheme(0)
            sage: f = X.hom([x^2*y, x^2*z, x*y*z], P2)
            sage: f.representatives()                                                   # needs sage.libs.singular
            [Scheme morphism:
               From: Closed subscheme of Projective Space of dimension 2
                     over Rational Field defined by: 0
               To:   Projective Space of dimension 2 over Rational Field
               Defn: Defined on coordinates by sending (x : y : z) to (x*y : x*z : y*z)]

        ::

            sage: P2.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: P1.<a,b> = ProjectiveSpace(QQ, 1)
            sage: X = P2.subscheme([x^2 - y^2 - y*z])
            sage: f = X.hom([x, y], P1)
            sage: f.representatives()                                                   # needs sage.libs.singular
            [Scheme morphism:
               From: Closed subscheme of Projective Space of dimension 2
                     over Rational Field defined by: x^2 - y^2 - y*z
               To:   Projective Space of dimension 1 over Rational Field
               Defn: Defined on coordinates by sending (x : y : z) to (y + z : x),
             Scheme morphism:
               From: Closed subscheme of Projective Space of dimension 2
                     over Rational Field defined by: x^2 - y^2 - y*z
               To:   Projective Space of dimension 1 over Rational Field
               Defn: Defined on coordinates by sending (x : y : z) to (x : y)]
            sage: g = _[0]                                                              # needs sage.libs.singular
            sage: g.representatives()                                                   # needs sage.libs.singular
            [Scheme morphism:
               From: Closed subscheme of Projective Space of dimension 2
                     over Rational Field defined by: x^2 - y^2 - y*z
               To:   Projective Space of dimension 1 over Rational Field
               Defn: Defined on coordinates by sending (x : y : z) to (y + z : x),
             Scheme morphism:
               From: Closed subscheme of Projective Space of dimension 2
                     over Rational Field defined by: x^2 - y^2 - y*z
               To:   Projective Space of dimension 1 over Rational Field
               Defn: Defined on coordinates by sending (x : y : z) to (x : y)]

        ::

            sage: P2.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: X = P2.subscheme([x^2 - y^2 - y*z])
            sage: A1.<a> = AffineSpace(QQ, 1)
            sage: g = X.hom([y/x], A1)
            sage: g.representatives()                                                   # needs sage.libs.singular
            [Scheme morphism:
               From: Closed subscheme of Projective Space of dimension 2
                     over Rational Field defined by: x^2 - y^2 - y*z
               To:   Affine Space of dimension 1 over Rational Field
               Defn: Defined on coordinates by sending (x : y : z) to (x/(y + z)),
             Scheme morphism:
               From: Closed subscheme of Projective Space of dimension 2
                     over Rational Field defined by: x^2 - y^2 - y*z
               To:   Affine Space of dimension 1 over Rational Field
               Defn: Defined on coordinates by sending (x : y : z) to (y/x)]
            sage: g0, g1 = _                                                            # needs sage.libs.singular
            sage: emb = A1.projective_embedding(0)
            sage: emb*g0                                                                # needs sage.libs.singular
            Scheme morphism:
              From: Closed subscheme of Projective Space of dimension 2
                    over Rational Field defined by: x^2 - y^2 - y*z
              To:   Projective Space of dimension 1 over Rational Field
              Defn: Defined on coordinates by sending (x : y : z) to (y + z : x)
            sage: emb*g1                                                                # needs sage.libs.singular
            Scheme morphism:
              From: Closed subscheme of Projective Space of dimension 2
                    over Rational Field defined by: x^2 - y^2 - y*z
              To:   Projective Space of dimension 1 over Rational Field
              Defn: Defined on coordinates by sending (x : y : z) to (x : y)

        ALGORITHM:

        The algorithm is from Proposition 1.1 in [Sim2004]_.
        """
        X = self.domain()
        Y = self.codomain()

        if not Y.is_projective():  # Y is affine
            emb = Y.projective_embedding(0)
            hom = self.parent()
            reprs = []
            for r in (emb*self).representatives():
                f0 = r[0]
                reprs.append(hom([f / f0 for f in r[1:]]))
            return reprs

        if not (X.base_ring() in _NumberFields or
                X.base_ring() in _FiniteFields):
            raise NotImplementedError("base ring {} is not supported by Singular".format(X.base_ring()))

        if not X.is_irreducible():
            raise ValueError("domain is not an irreducible scheme")

        # prepare homogeneous coordinate ring of X in Singular
        from sage.interfaces.singular import singular
        from sage.rings.polynomial.term_order import TermOrder
        T = TermOrder('degrevlex')
        T._singular_ringorder_column = 1  # (c,dp) in Singular
        S = X.ambient_space().coordinate_ring().change_ring(order=T)
        R = S.quotient_ring(X.defining_ideal().change_ring(S))

        if R is S:  # true when the defining ideal is zero
            def lift(x):
                return x.numerator()
        else:  # R is an ordinary quotient ring
            def lift(x):
                return x.lift()

        F = [R(f) for f in self.defining_polynomials()]
        n = len(F)
        I = R.ideal(F)

        # find r with nonzero F[r]
        r = 0
        while not F[r]:
            r = r + 1

        # This is a minimal free presentation of the ideal I:
        #
        #          phi         F
        #    R^m ------> R^n -----> I
        #
        # where n is the number of defining polynomials and m is the number of
        # syzygy relations of I.

        # compute the kernel of the transpose of phi in Singular
        phi = I._singular_().syz()
        phi_trans = phi.matrix().transpose()
        kernel = singular.modulo(phi_trans, singular.module())

        M = kernel.sage_matrix(R)  # m * n matrix over R
        reprs = []
        for i in range(M.ncols()):
            lifts = [lift(F[j] * M[r][i] / F[r]) for j in range(n)]
            reprs.append(X.hom(lifts, Y))

        return reprs

    def indeterminacy_locus(self):
        """
        Return the indeterminacy locus of this map.

        The map defines a rational map on the domain. The output is the
        subscheme of the domain on which the rational map is not defined by any
        representative of the rational map. See :meth:`representatives()`.

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: P2.<x0,x1,x2> = ProjectiveSpace(QQ, 2)
            sage: X = P2.subscheme(0)
            sage: f = X.hom([x1,x0], P)
            sage: L = f.indeterminacy_locus()                                           # needs sage.libs.singular
            sage: L.rational_points()                                                   # needs sage.libs.singular
            [(0 : 0 : 1)]

        ::

            sage: P2.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: P1.<a,b> = ProjectiveSpace(QQ, 1)
            sage: X = P2.subscheme([x^2 - y^2 - y*z])
            sage: f = X.hom([x,y], P1)
            sage: f.indeterminacy_locus()                                               # needs sage.libs.singular
            Closed subscheme of Projective Space of dimension 2 over Rational Field defined by:
              z,
              y,
              x

        ::

            sage: P3.<x,y,z,w> = ProjectiveSpace(QQ, 3)
            sage: P2.<a,b,c> = ProjectiveSpace(QQ, 2)
            sage: X = P3.subscheme(x^2 - w*y - x*z)
            sage: f = X.hom([x*y, y*z, z*x], P2)
            sage: L = f.indeterminacy_locus()                                           # needs sage.libs.singular
            sage: L.dimension()                                                         # needs sage.libs.singular
            0
            sage: L.degree()                                                            # needs sage.libs.singular
            2
            sage: L.rational_points()                                                   # needs sage.libs.singular
            [(0 : 0 : 0 : 1), (0 : 1 : 0 : 0)]

        ::

            sage: P3.<x,y,z,w> = ProjectiveSpace(QQ, 3)
            sage: A2.<a,b> = AffineSpace(QQ, 2)
            sage: X = P3.subscheme(x^2 - w*y - x*z)
            sage: f = X.hom([x/z, y/x], A2)
            sage: L = f.indeterminacy_locus()                                           # needs sage.libs.singular
            sage: L.rational_points()                                                   # needs sage.libs.singular
            [(0 : 0 : 0 : 1), (0 : 1 : 0 : 0)]

        ::

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: X = P.subscheme(x - y)
            sage: H = End(X)
            sage: f = H([x^2 - 4*y^2, y^2 - z^2, 4*z^2 - x^2])
            sage: Z = f.indeterminacy_locus(); Z                                        # needs sage.libs.singular
            Closed subscheme of Projective Space of dimension 2 over Rational Field defined by:
              z,
              y,
              x
        """
        X = self.domain()
        Y = self.codomain()
        Amb = X.ambient_space()

        if not X.is_irreducible():
            components = X.irreducible_components()

            def self_with_domain(C):
                return self*C.hom(Amb.gens(), X)

            locus = self_with_domain(components[0]).indeterminacy_locus()
            for C in components[1:]:
                locus = locus.union(self_with_domain(C).indeterminacy_locus())

            return locus

        if not Y.is_projective():  # Y is affine
            emb = Y.projective_embedding(0)
        else:
            emb = None

        polys = list(X.defining_polynomials())

        for r in self.representatives():
            r_proj = r if emb is None else emb * r
            polys.extend(r_proj)

        return Amb.subscheme(polys).reduce()

    def is_morphism(self):
        """
        Return ``True`` if the map is defined everywhere on the domain.

        EXAMPLES::

            sage: P2.<x,y,z> = ProjectiveSpace(QQ,2)
            sage: P1.<a,b> = ProjectiveSpace(QQ, 1)
            sage: X = P2.subscheme([x^2 - y^2 - y*z])
            sage: f = X.hom([x,y], P1)
            sage: f.is_morphism()                                                       # needs sage.libs.singular
            True
        """
        return self.indeterminacy_locus().dimension() < 0

    def image(self):
        """
        Return the scheme-theoretic image of the morphism.

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: P2.<x0,x1,x2> = ProjectiveSpace(QQ, 2)
            sage: X = P2.subscheme(0)
            sage: f = X.hom([x1,x0], P)
            sage: f.image()                                                             # needs sage.libs.singular
            Closed subscheme of Projective Space of dimension 1 over Rational Field defined by:
              (no polynomials)

        ::

            sage: P2.<x,y,z> = ProjectiveSpace(QQ,2)
            sage: X = P2.subscheme([z^3 - x*y^2 + y^3])
            sage: f = X.hom([x*z, x*y, x^2 + y*z], P2)
            sage: f.image()                                                             # needs sage.libs.singular
            Closed subscheme of Projective Space of dimension 2 over Rational Field defined by:
              x^6 + 2*x^3*y^3 + x*y^5 + y^6 - x^3*y^2*z - y^5*z
        """
        X = self.domain()
        Y = self.codomain()

        if not Y.is_projective():
            e = Y.projective_embedding(0)
            return (e * self).image().affine_patch(0, Y.ambient_space())

        k = self.base_ring()

        AX = X.ambient_space()
        AY = Y.ambient_space()

        S = AX.coordinate_ring()
        T = AY.coordinate_ring()

        n = S.ngens()
        m = T.ngens()

        dummy_names = ['d{}__'.format(i) for i in range(m)]
        D = PolynomialRing(k, names=dummy_names)

        names = list(S.variable_names()) + dummy_names  # this order of variables is important
        R = PolynomialRing(k, names=names, order='degrevlex({}),degrevlex({})'.format(m,n))

        # compute the ideal of the image by elimination
        i = R.ideal(list(X.defining_ideal().gens()) + [self._polys[i] - R.gen(n + i) for i in range(m)])
        j = [g for g in i.groebner_basis() if g in D]

        gens = [g.subs(dict(zip(R.gens()[n:],T.gens()))) for g in j]
        return AY.subscheme(gens)

    @cached_method
    def graph(self):
        """
        Return the graph of this morphism.

        The graph is a subscheme of the product of the ambient spaces of the
        domain and the codomain. If the ambient space of the codomain is an
        affine space, it is first embedded into a projective space.

        EXAMPLES:

        We get the standard quadratic curve as the graph of a quadratic function
        of an affine line. ::

            sage: A1.<x> = AffineSpace(1, QQ)
            sage: X = A1.subscheme(0)  # affine line
            sage: phi = X.hom([x^2], A1)
            sage: mor = phi.homogenize(0)                                               # needs sage.libs.singular
            sage: G = mor.graph(); G                                                    # needs sage.libs.singular
            Closed subscheme of Product of projective spaces P^1 x P^1
              over Rational Field defined by: x1^2*x2 - x0^2*x3
            sage: G.affine_patch([0, 0])                                                # needs sage.libs.singular
            Closed subscheme of Affine Space of dimension 2
              over Rational Field defined by: x0^2 - x1
        """
        X = self.domain()
        Y = self.codomain()

        if not Y.is_projective():
            e = Y.projective_embedding(0)
            return (e * self).graph()

        AX = X.ambient_space()
        AY = Y.ambient_space()

        n = AX.dimension()
        m = AY.dimension()

        if any(v in AX.variable_names() for v in AY.variable_names()):
            from sage.schemes.product_projective.space import ProductProjectiveSpaces
            AXY = ProductProjectiveSpaces([n, m], self.base_ring())
        else:
            AXY = AX * AY  # product of projective spaces

        R = AXY.coordinate_ring()
        F = [R(f) for f in self.defining_polynomials()]
        g = R.gens()

        # Suppose R = k[x_0, ..., x_n, y_0, ..., y_m]. Then the bihomogeneous
        # ideal of the graph is
        #
        #    I + (y_iF_j - y_jF_i : 0 <= i, j <= m)
        #
        # saturated with respect to (F_0, F_1, ..., F_m).
        n1 = n + 1
        m1 = m + 1
        I = X.defining_ideal().change_ring(R)
        h = [g[n1 + i] * F[j] - g[n1 + j] * F[i] for i in range(m1) for j in range(i + 1, m1)]
        J, _ = (I + R.ideal(h)).saturation(R.ideal(F))

        return AXY.subscheme(J)

    @cached_method
    def projective_degrees(self):
        """
        Return the projective degrees of this rational map.

        EXAMPLES::

            sage: # needs sage.schemes
            sage: k = GF(11)
            sage: E = EllipticCurve(k, [1,1])
            sage: Q = E(6, 5)
            sage: phi = E.scalar_multiplication(2)
            sage: mor = phi.as_morphism()
            sage: mor.projective_degrees()
            (12, 3)
        """
        X = self.domain()
        Y = self.codomain()

        if not Y.is_projective():  # Y is affine
            e = Y.projective_embedding(0)
            return (e * self).projective_degrees()

        AX = X.ambient_space()
        AY = Y.ambient_space()
        xn = AX.ngens()
        yn = AY.ngens()

        G = self.graph()
        I = G.defining_ideal()  # a bihomogeneous ideal

        from sage.modules.free_module_element import vector

        degrees = xn * [vector([1, 0])] + yn * [vector([0, 1])]
        res = I.graded_free_resolution(degrees=degrees, algorithm='shreyer')
        kpoly = res.K_polynomial()

        L = kpoly.parent()
        t1, t2 = L.gens()
        poly = kpoly.substitute({t1: 1 - t1, t2: 1 - t2})

        n = AX.dimension()
        m = AY.dimension()
        k = X.dimension()
        return tuple(poly.monomial_coefficient(L.monomial(n - i, m - k + i)) for i in range(k + 1))

    def degree(self):
        """
        Return the degree of this rational map.

        EXAMPLES::

            sage: # needs sage.schemes
            sage: k = GF(11)
            sage: E = EllipticCurve(k, [1,1])
            sage: Q = E(6, 5)
            sage: phi = E.scalar_multiplication(2)
            sage: mor = phi.as_morphism()
            sage: mor.degree()
            4
        """
        return self.projective_degrees()[0] // self.image().degree()
