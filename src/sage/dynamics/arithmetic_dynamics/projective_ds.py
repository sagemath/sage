r"""
Dynamical systems on projective schemes

A dynamical system of projective schemes determined by homogeneous
polynomials functions that define what the morphism does on points
in the ambient projective space.

The main constructor functions are given by :class:`DynamicalSystem` and
:class:`DynamicalSystem_projective`. The constructors function can take either
polynomials or a morphism from which to construct a dynamical system.
If the domain is not specified, it is constructed. However, if you plan on
working with points or subvarieties in the domain, it recommended to specify
the domain.

The initialization checks are always performed by the constructor functions.
It is possible, but not recommended, to skip these checks by calling the
class initialization directly.

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

- Ben Hutz (2017-7): relocate code and create class
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
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from copy import copy
from itertools import count, product

import sage.rings.abc

from sage.arith.functions import lcm
from sage.arith.misc import binomial, gcd, integer_ceil as ceil, is_prime, moebius, next_prime, primes
from sage.calculus.functions import jacobian
from sage.categories.fields import Fields
from sage.categories.finite_fields import FiniteFields
from sage.categories.function_fields import FunctionFields
from sage.categories.homset import End
from sage.categories.number_fields import NumberFields
from sage.dynamics.arithmetic_dynamics.endPN_automorphism_group import (
    automorphism_group_QQ_CRT,
    automorphism_group_QQ_fixedpoints,
    conjugating_set_helper,
    conjugating_set_initializer,
    is_conjugate_helper)
from sage.dynamics.arithmetic_dynamics.endPN_automorphism_group import automorphism_group_FF
from sage.dynamics.arithmetic_dynamics.generic_ds import DynamicalSystem
from sage.dynamics.arithmetic_dynamics.projective_ds_helper import (
    _fast_possible_periods,
    _all_periodic_points)
from sage.matrix.constructor import matrix, identity_matrix
from sage.misc.cachefunc import cached_method
from sage.misc.classcall_metaclass import typecall
from sage.misc.functional import sqrt
from sage.misc.lazy_import import lazy_import
from sage.misc.mrange import xmrange
from sage.modules.free_module_element import vector
from sage.parallel.ncpus import ncpus
from sage.parallel.use_fork import p_iter_fork
from sage.rings.complex_mpfr import ComplexField
from sage.rings.finite_rings.finite_field_base import FiniteField
from sage.rings.finite_rings.finite_field_constructor import GF
from sage.rings.finite_rings.integer_mod_ring import Zmod
from sage.rings.fraction_field import FractionField, FractionField_generic, FractionField_1poly_field
from sage.rings.fraction_field_element import FractionFieldElement
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.flatten import FlatteningMorphism, UnflatteningMorphism
from sage.rings.morphism import RingHomomorphism_im_gens
from sage.rings.polynomial.multi_polynomial_ring_base import MPolynomialRing_base
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.polynomial.polynomial_ring import PolynomialRing_generic
from sage.rings.quotient_ring import QuotientRing_generic
from sage.rings.rational_field import QQ
from sage.rings.real_mpfr import RealField
from sage.schemes.generic.morphism import SchemeMorphism_polynomial
from sage.schemes.product_projective.space import ProductProjectiveSpaces_ring
from sage.schemes.projective.projective_morphism import (
    SchemeMorphism_polynomial_projective_space,
    SchemeMorphism_polynomial_projective_space_field,
    SchemeMorphism_polynomial_projective_space_finite_field)
from sage.schemes.projective.projective_space import ProjectiveSpace, ProjectiveSpace_ring
from sage.schemes.projective.projective_subscheme import AlgebraicScheme_subscheme_projective
from sage.structure.element import get_coercion_model
from sage.rings.qqbar import QQbar, number_field_elements_from_algebraics
from sage.schemes.elliptic_curves.constructor import EllipticCurve

lazy_import('sage.rings.algebraic_closure_finite_field', 'AlgebraicClosureFiniteField_generic')
lazy_import('sage.rings.number_field.number_field_ideal', 'NumberFieldFractionalIdeal')
lazy_import('sage.rings.padics.factory', 'Qp')
lazy_import('sage.rings.qqbar', 'number_field_elements_from_algebraics')

from sage.libs.pari.all import PariError


class DynamicalSystem_projective(SchemeMorphism_polynomial_projective_space,
                                      DynamicalSystem):
    r"""A dynamical system of projective schemes determined by homogeneous
    polynomials that define what the morphism does on points in the
    ambient projective space.

    .. WARNING::

        You should not create objects of this class directly because
        no type or consistency checking is performed. The preferred
        method to construct such dynamical systems is to use
        :func:`~sage.dynamics.arithmetic_dynamics.generic_ds.DynamicalSystem_projective`
        function

    INPUT:

    - ``morphism_or_polys`` -- a SchemeMorphism, a polynomial, a
      rational function, or a list or tuple of homogeneous polynomials

    - ``domain`` -- (optional) projective space or projective subscheme

    - ``names`` -- tuple of strings (default: ``'X','Y'``) to be used as coordinate
      names for a projective space that is constructed

      The following combinations of ``morphism_or_polys`` and
      ``domain`` are meaningful:

      * ``morphism_or_polys`` is a SchemeMorphism; ``domain`` is
        ignored in this case.

      * ``morphism_or_polys`` is a list of homogeneous polynomials
        that define a rational endomorphism of ``domain``.

      * ``morphism_or_polys`` is a list of homogeneous polynomials and
        ``domain`` is unspecified; ``domain`` is then taken to be the
        projective space of appropriate dimension over the common base ring,
        if one exists, of the elements of ``morphism_or_polys``.

      * ``morphism_or_polys`` is a single polynomial or rational
        function; ``domain`` is ignored and taken to be a
        1-dimensional projective space over the base ring of
        ``morphism_or_polys`` with coordinate names given by ``names``.

    OUTPUT: :class:`DynamicalSystem_projective`

    EXAMPLES::

        sage: P1.<x,y> = ProjectiveSpace(QQ,1)
        sage: DynamicalSystem_projective([y, 2*x])
        Dynamical System of Projective Space of dimension 1 over Rational Field
          Defn: Defined on coordinates by sending (x : y) to
                (y : 2*x)

    We can define dynamical systems on `P^1` by giving a polynomial or
    rational function::

        sage: R.<t> = QQ[]
        sage: DynamicalSystem_projective(t^2 - 3)
        Dynamical System of Projective Space of dimension 1 over Rational Field
          Defn: Defined on coordinates by sending (X : Y) to
                (X^2 - 3*Y^2 : Y^2)
        sage: DynamicalSystem_projective(1/t^2)
        Dynamical System of Projective Space of dimension 1 over Rational Field
          Defn: Defined on coordinates by sending (X : Y) to
                (Y^2 : X^2)

    ::

        sage: R.<x> = PolynomialRing(QQ,1)
        sage: DynamicalSystem_projective(x^2, names=['a','b'])
        Dynamical System of Projective Space of dimension 1 over Rational Field
          Defn: Defined on coordinates by sending (a : b) to
                (a^2 : b^2)

    Symbolic Ring elements are not allowed::

        sage: x,y = var('x,y')                                                          # needs sage.symbolic
        sage: DynamicalSystem_projective([x^2, y^2])                                    # needs sage.symbolic
        Traceback (most recent call last):
        ...
        ValueError: [x^2, y^2] must be elements of a polynomial ring

    ::

        sage: R.<x> = PolynomialRing(QQ,1)
        sage: DynamicalSystem_projective(x^2)
        Dynamical System of Projective Space of dimension 1 over Rational Field
          Defn: Defined on coordinates by sending (X : Y) to
                (X^2 : Y^2)

    ::

        sage: R.<t> = PolynomialRing(QQ)
        sage: P.<x,y,z> = ProjectiveSpace(R, 2)
        sage: X = P.subscheme([x])
        sage: DynamicalSystem_projective([x^2, t*y^2, x*z], domain=X)
        Dynamical System of Closed subscheme of Projective Space of dimension
        2 over Univariate Polynomial Ring in t over Rational Field defined by:
          x
          Defn: Defined on coordinates by sending (x : y : z) to
                (x^2 : t*y^2 : x*z)

    When elements of the quotient ring are used, they are reduced::

        sage: P.<x,y,z> = ProjectiveSpace(CC, 2)
        sage: X = P.subscheme([x - y])
        sage: u,v,w = X.coordinate_ring().gens()                                        # needs sage.rings.function_field
        sage: DynamicalSystem_projective([u^2, v^2, w*u], domain=X)                     # needs sage.rings.function_field
        Dynamical System of Closed subscheme of Projective Space of dimension
        2 over Complex Field with 53 bits of precision defined by:
          x - y
          Defn: Defined on coordinates by sending (x : y : z) to
                (y^2 : y^2 : y*z)

    We can also compute the forward image of subschemes through
    elimination. In particular, let `X = V(h_1,\ldots, h_t)` and define the ideal
    `I = (h_1,\ldots,h_t,y_0-f_0(\bar{x}), \ldots, y_n-f_n(\bar{x}))`.
    Then the elimination ideal `I_{n+1} = I \cap K[y_0,\ldots,y_n]` is a homogeneous
    ideal and `f(X) = V(I_{n+1})`::

        sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
        sage: f = DynamicalSystem_projective([(x-2*y)^2, (x-2*z)^2, x^2])
        sage: X = P.subscheme(y - z)
        sage: f(f(f(X)))                                                                # needs sage.rings.function_field
        Closed subscheme of Projective Space of dimension 2 over Rational Field
        defined by:
          y - z

    ::

        sage: P.<x,y,z,w> = ProjectiveSpace(QQ, 3)
        sage: f = DynamicalSystem_projective([(x-2*y)^2, (x-2*z)^2, (x-2*w)^2, x^2])
        sage: f(P.subscheme([x, y, z]))                                                 # needs sage.rings.function_field
        Closed subscheme of Projective Space of dimension 3 over Rational Field
        defined by:
          w,
          y,
          x

    ::

        sage: T.<x,y,z,w,u> = ProductProjectiveSpaces([2, 1], QQ)
        sage: DynamicalSystem_projective([x^2*u, y^2*w, z^2*u, w^2, u^2], domain=T)
        Dynamical System of Product of projective spaces P^2 x P^1 over Rational Field
          Defn: Defined by sending (x : y : z , w : u) to
                (x^2*u : y^2*w : z^2*u , w^2 : u^2).

    ::

        sage: # needs sage.rings.number_field
        sage: K.<v> = QuadraticField(-7)
        sage: P.<x,y> = ProjectiveSpace(K, 1)
        sage: f = DynamicalSystem([x^3 + v*x*y^2, y^3])
        sage: fbar = f.change_ring(QQbar)
        sage: fbar.is_postcritically_finite()
        False
    """

    @staticmethod
    def __classcall_private__(cls, morphism_or_polys, domain=None, names=None):
        r"""
        Return the appropriate dynamical system on a projective scheme.

        TESTS::

            sage: R.<x,y> = QQ[]
            sage: P1 = ProjectiveSpace(R)
            sage: f = DynamicalSystem_projective([x - y, x*y])
            Traceback (most recent call last):
            ...
            ValueError: polys (=[x - y, x*y]) must be of the same degree
            sage: DynamicalSystem_projective([x - 1, x*y + x])
            Traceback (most recent call last):
            ...
            ValueError: polys (=[x - 1, x*y + x]) must be homogeneous

        ::

            sage: DynamicalSystem_projective([exp(x), exp(y)])                          # needs sage.symbolic
            Traceback (most recent call last):
            ...
            ValueError: [e^x, e^y] must be elements of a polynomial ring

        ::

            sage: T.<x,y,z,w,u> = ProductProjectiveSpaces([2, 1], QQ)
            sage: DynamicalSystem_projective([x^2*u, y^2*w, z^2*u, w^2, u*z], domain=T)
            Traceback (most recent call last):
            ...
            TypeError: polys (=[x^2*u, y^2*w, z^2*u, w^2, z*u]) must be
            multi-homogeneous of the same degrees (by component)

        ::

            sage: A.<x,y> = AffineSpace(ZZ, 2)
            sage: DynamicalSystem_projective([x^2, y^2], A)
            Traceback (most recent call last):
            ...
            ValueError: "domain" must be a projective scheme
            sage: H = End(A)
            sage: f = H([x,y])
            sage: DynamicalSystem_projective(f)
            Traceback (most recent call last):
            ...
            ValueError: "domain" must be a projective scheme

        ::

            sage: R.<x> = PolynomialRing(QQ,1)
            sage: DynamicalSystem_projective(x^2, names='t')
            Traceback (most recent call last):
            ...
            ValueError: specify 2 variable names

        ::

            sage: P1.<x,y> = ProjectiveSpace(QQ,1)
            sage: DynamicalSystem_projective([y, x, y], domain=P1)
            Traceback (most recent call last):
            ...
            ValueError: number of polys does not match dimension of Projective Space of dimension 1 over Rational Field

        ::

            sage: A.<x,y> = AffineSpace(QQ,2)
            sage: DynamicalSystem_projective([y,x], domain=A)
            Traceback (most recent call last):
            ...
            ValueError: "domain" must be a projective scheme

        ::

            sage: R.<x> = QQ[]
            sage: DynamicalSystem([x^2])
            Traceback (most recent call last):
            ...
            ValueError: list/tuple must have at least 2 polynomials

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: f = DynamicalSystem([CC.0*x^2 + 2*y^2, 1*y^2], domain=P)
            Traceback (most recent call last):
            ...
            TypeError: coefficients of polynomial not in Rational Field
        """
        from sage.dynamics.arithmetic_dynamics.product_projective_ds import DynamicalSystem_product_projective

        if isinstance(morphism_or_polys, SchemeMorphism_polynomial):
            R = morphism_or_polys.base_ring()
            domain = morphism_or_polys.domain()
            polys = list(morphism_or_polys)
            if domain != morphism_or_polys.codomain():
                raise ValueError('domain and codomain do not agree')
            if not isinstance(domain, ProjectiveSpace_ring) and not isinstance(domain, AlgebraicScheme_subscheme_projective):
                raise ValueError('"domain" must be a projective scheme')
            if R not in Fields():
                return typecall(cls, polys, domain)
            if isinstance(R, FiniteField):
                return DynamicalSystem_projective_finite_field(polys, domain)
            return DynamicalSystem_projective_field(polys, domain)

        if isinstance(morphism_or_polys, (list, tuple)):
            polys = list(morphism_or_polys)
            if len(polys) == 1:
                raise ValueError("list/tuple must have at least 2 polynomials")
            test = lambda x: isinstance(x, (PolynomialRing_generic, MPolynomialRing_base))
            if not all(test(poly.parent()) for poly in polys):
                try:
                    polys = [poly.lift() for poly in polys]
                except AttributeError:
                    raise ValueError('{} must be elements of a polynomial ring'.format(morphism_or_polys))
        else:
            # homogenize!
            f = morphism_or_polys
            aff_CR = f.parent()
            if (not isinstance(aff_CR, (PolynomialRing_generic, FractionField_generic))
                and not (isinstance(aff_CR, MPolynomialRing_base) and aff_CR.ngens() == 1)):
                msg = '{} is not a single variable polynomial or rational function'
                raise ValueError(msg.format(f))
            if isinstance(aff_CR, FractionField_generic):
                polys = [f.numerator(),f.denominator()]
            else:
                polys = [f, aff_CR(1)]
            d = max(poly.degree() for poly in polys)
            if names is None:
                names = ('X','Y')
            elif len(names) != 2:
                raise ValueError('specify 2 variable names')
            proj_CR = PolynomialRing(aff_CR.base_ring(), names=names)
            X,Y = proj_CR.gens()
            polys = [proj_CR(Y**d * poly(X/Y)) for poly in polys]

        if domain is None:
            PR = get_coercion_model().common_parent(*polys)
            polys = [PR(poly) for poly in polys]
            domain = ProjectiveSpace(PR)
        else:
            # Check if we can coerce the given polynomials over the given domain
            PR = domain.ambient_space().coordinate_ring()
            try:
                polys = [PR(poly) for poly in polys]
            except TypeError:
                raise TypeError('coefficients of polynomial not in {}'.format(domain.base_ring()))
        if len(polys) != domain.ambient_space().coordinate_ring().ngens():
            raise ValueError(f'number of polys does not match dimension of {domain}')
        R = domain.base_ring()
        if isinstance(R, sage.rings.abc.SymbolicRing):
            raise TypeError("the base ring cannot be the Symbolic Ring or a symbolic subring")

        if isinstance(domain, ProductProjectiveSpaces_ring):
            splitpolys = domain._factors(polys)
            for split_poly in splitpolys:
                split_d = domain._degree(split_poly[0])
                if not all(split_d == domain._degree(f) for f in split_poly):
                    msg = 'polys (={}) must be multi-homogeneous of the same degrees (by component)'
                    raise TypeError(msg.format(polys))
            if isinstance(R, FiniteField):
                from sage.dynamics.arithmetic_dynamics.product_projective_ds import DynamicalSystem_product_projective_finite_field
                return DynamicalSystem_product_projective_finite_field(polys, domain)
            return DynamicalSystem_product_projective(polys, domain)

        # Now polys define an endomorphism of a scheme in P^n
        if not all(poly.is_homogeneous() for poly in polys):
            msg = 'polys (={}) must be homogeneous'
            raise ValueError(msg.format(polys))
        d = polys[0].degree()
        if not all(poly.degree() == d for poly in polys):
            msg = 'polys (={}) must be of the same degree'
            raise ValueError(msg.format(polys))

        if not isinstance(domain, (ProjectiveSpace_ring,
                                   AlgebraicScheme_subscheme_projective)):
            raise ValueError('"domain" must be a projective scheme')
        if R not in Fields():
            return typecall(cls, polys, domain)
        if isinstance(R, FiniteField):
                return DynamicalSystem_projective_finite_field(polys, domain)
        return DynamicalSystem_projective_field(polys, domain)

    def __init__(self, polys, domain):
        r"""
        The Python constructor.

        See :class:`DynamicalSystem` for details.

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: DynamicalSystem_projective([3/5*x^2, y^2], domain=P)
            Dynamical System of Projective Space of dimension 1 over Rational Field
              Defn: Defined on coordinates by sending (x : y) to
                    (3/5*x^2 : y^2)
        """
        # Next attribute needed for _fast_eval and _fastpolys
        R = polys[0].base_ring()
        self._is_prime_finite_field = isinstance(R, FiniteField) and R.degree() == 1
        DynamicalSystem.__init__(self, polys, domain)

    def __copy__(self):
        r"""
        Return a copy of this dynamical system.

        OUTPUT: :class:`DynamicalSystem_projective`

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: f = DynamicalSystem_projective([3/5*x^2,6*y^2])
            sage: g = copy(f)
            sage: f == g
            True
            sage: f is g
            False
        """
        return DynamicalSystem_projective(self._polys, self.domain())

    def _number_field_from_algebraics(self):
        r"""
        Return a dynamical system defined over the number field of its coefficients.

        OUTPUT: dynamical system

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQbar,1)                                    # needs sage.rings.number_field
            sage: f = DynamicalSystem_projective([x^2 + QQbar(sqrt(2)) * y^2, y^2])     # needs sage.rings.number_field sage.symbolic
            sage: f._number_field_from_algebraics()                                     # needs sage.rings.number_field sage.symbolic
            Dynamical System of Projective Space of dimension 1 over Number Field in a
             with defining polynomial y^2 - 2 with a = 1.414213562373095?
              Defn: Defined on coordinates by sending (x : y) to
                    (x^2 + a*y^2 : y^2)
        """
        return self.as_scheme_morphism()._number_field_from_algebraics().as_dynamical_system()

    def dehomogenize(self, n):
        r"""
        Return the standard dehomogenization at the ``n[0]`` coordinate
        for the domain and the ``n[1]`` coordinate for the codomain.

        Note that the new function is defined over the fraction field
        of the base ring of this map.

        INPUT:

        - ``n`` -- tuple of nonnegative integers; if ``n`` is an integer,
          then the two values of the tuple are assumed to be the same

        OUTPUT:

        If the dehomogenizing indices are the same for the domain and
        codomain, then a :class:`DynamicalSystem_affine` given by
        dehomogenizing the source and target of ``self`` with respect to
        the given indices is returned. If the dehomogenizing indices
        for the domain and codomain are different then the resulting
        affine patches are different and a scheme morphism is returned.

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(ZZ,1)
            sage: f = DynamicalSystem_projective([x^2 + y^2, y^2])
            sage: f.dehomogenize(0)
            Dynamical System of Affine Space of dimension 1 over Integer Ring
              Defn: Defined on coordinates by sending (y) to
                    (y^2/(y^2 + 1))
            sage: f.dehomogenize((0, 1))
            Scheme morphism:
              From: Affine Space of dimension 1 over Integer Ring
              To:   Affine Space of dimension 1 over Integer Ring
              Defn: Defined on coordinates by sending (y) to
                    ((y^2 + 1)/y^2)
        """
        F = self.as_scheme_morphism().dehomogenize(n)
        if F.domain() == F.codomain():
            return F.as_dynamical_system()
        else:
            return F

    def dynatomic_polynomial(self, period):
        r"""
        For a dynamical system of `\mathbb{P}^1` compute the dynatomic
        polynomial.

        The dynatomic polynomial is the analog of the cyclotomic
        polynomial and its roots are the points of formal period `period`.
        If possible the division is done in the coordinate ring of this
        map and a polynomial is returned. In rings where that is not
        possible, a :class:`FractionField` element will be returned.
        In certain cases, when the conversion back to a polynomial fails,
        a :class:`SymbolRing` element will be returned.

        ALGORITHM:

        For a positive integer `n`, let `[F_n,G_n]` be the coordinates of the `n`-th
        iterate of `f`. Then construct

        .. MATH::

            \Phi^{\ast}_n(f)(x,y) = \sum_{d \mid n}
                (yF_d(x,y) - xG_d(x,y))^{\mu(n/d)},

        where `\mu` is the Möbius function.

        For a pair `[m,n]`, let `f^m = [F_m,G_m]`. Compute

        .. MATH::

            \Phi^{\ast}_{m,n}(f)(x,y) = \Phi^{\ast}_n(f)(F_m,G_m) /
                \Phi^{\ast}_n(f)(F_{m-1},G_{m-1})

        REFERENCES:

        - [Hutz2015]_
        - [MoPa1994]_

        INPUT:

        - ``period`` -- positive integer or a list/tuple `[m,n]` where
          `m` is the preperiod and `n` is the period

        OUTPUT:

        If possible, a two variable polynomial in the coordinate ring
        of this map. Otherwise a fraction field element of the coordinate
        ring of this map. Or, a :class:`SymbolicRing` element.

        .. TODO::

            - Do the division when the base ring is `p`-adic so that
              the output is a polynomial.

            - Convert back to a polynomial when the base ring is a
              function field (not over `\QQ` or `F_p`).

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: f = DynamicalSystem_projective([x^2 + y^2, y^2])
            sage: f.dynatomic_polynomial(2)                                             # needs sage.libs.pari
            x^2 + x*y + 2*y^2

        ::

            sage: P.<x,y> = ProjectiveSpace(ZZ,1)
            sage: f = DynamicalSystem_projective([x^2 + y^2, x*y])
            sage: f.dynatomic_polynomial(4)                                             # needs sage.libs.pari
            2*x^12 + 18*x^10*y^2 + 57*x^8*y^4 + 79*x^6*y^6 + 48*x^4*y^8 + 12*x^2*y^10 + y^12

        ::

            sage: P.<x,y> = ProjectiveSpace(CC,1)
            sage: f = DynamicalSystem_projective([x^2 + y^2, 3*x*y])
            sage: f.dynatomic_polynomial(3)                                             # needs sage.libs.pari
            13.0000000000000*x^6 + 117.000000000000*x^4*y^2 +
            78.0000000000000*x^2*y^4 + y^6

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: f = DynamicalSystem_projective([x^2 - 10/9*y^2, y^2])
            sage: f.dynatomic_polynomial([2,1])
            x^4*y^2 - 11/9*x^2*y^4 - 80/81*y^6

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: f = DynamicalSystem_projective([x^2 - 29/16*y^2, y^2])
            sage: f.dynatomic_polynomial([2,3])                                         # needs sage.libs.pari
            x^12 - 95/8*x^10*y^2 + 13799/256*x^8*y^4 - 119953/1024*x^6*y^6 +
            8198847/65536*x^4*y^8 - 31492431/524288*x^2*y^10 +
            172692729/16777216*y^12

        ::

            sage: P.<x,y> = ProjectiveSpace(ZZ,1)
            sage: f = DynamicalSystem_projective([x^2 - y^2, y^2])
            sage: f.dynatomic_polynomial([1,2])                                         # needs sage.libs.pari
            x^2 - x*y

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: f = DynamicalSystem_projective([x^3 - y^3, 3*x*y^2])
            sage: f.dynatomic_polynomial([0,4])==f.dynatomic_polynomial(4)              # needs sage.libs.pari
            True

        ::

            sage: P.<x,y,z> = ProjectiveSpace(QQ,2)
            sage: f = DynamicalSystem_projective([x^2 + y^2, x*y, z^2])
            sage: f.dynatomic_polynomial(2)
            Traceback (most recent call last):
            ...
            TypeError: does not make sense in dimension >1

        ::

            sage: P.<x,y> = ProjectiveSpace(Qp(5),1)                                    # needs sage.rings.padics
            sage: f = DynamicalSystem_projective([x^2 + y^2, y^2])                      # needs sage.rings.padics
            sage: f.dynatomic_polynomial(2)                                             # needs sage.rings.padics
            (x^4*y + (2 + O(5^20))*x^2*y^3 - x*y^4 + (2 + O(5^20))*y^5)/(x^2*y - x*y^2 + y^3)

        ::

            sage: L.<t> = PolynomialRing(QQ)
            sage: P.<x,y> = ProjectiveSpace(L,1)
            sage: f = DynamicalSystem_projective([x^2 + t*y^2, y^2])
            sage: f.dynatomic_polynomial(2)                                             # needs sage.libs.pari
            x^2 + x*y + (t + 1)*y^2

        ::

            sage: K.<c> = PolynomialRing(ZZ)
            sage: P.<x,y> = ProjectiveSpace(K,1)
            sage: f = DynamicalSystem_projective([x^2 + c*y^2, y^2])
            sage: f.dynatomic_polynomial([1, 2])                                        # needs sage.libs.pari
            x^2 - x*y + (c + 1)*y^2

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: f = DynamicalSystem_projective([x^2 + y^2, y^2])
            sage: f.dynatomic_polynomial(2)                                             # needs sage.libs.pari
            x^2 + x*y + 2*y^2
            sage: R.<X> = PolynomialRing(QQ)
            sage: K.<c> = NumberField(X^2 + X + 2)                                      # needs sage.rings.number_field
            sage: PP = P.change_ring(K)
            sage: ff = f.change_ring(K)
            sage: p = PP((c, 1))
            sage: ff(ff(p)) == p
            True

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: f = DynamicalSystem_projective([x^2 + y^2, x*y])
            sage: f.dynatomic_polynomial([2, 2])                                        # needs sage.libs.pari
            x^4 + 4*x^2*y^2 + y^4
            sage: R.<X> = PolynomialRing(QQ)
            sage: K.<c> = NumberField(X^4 + 4*X^2 + 1)                                  # needs sage.rings.number_field
            sage: PP = P.change_ring(K)
            sage: ff = f.change_ring(K)
            sage: p = PP((c, 1))
            sage: ff.nth_iterate(p, 4) == ff.nth_iterate(p, 2)                          # needs sage.rings.number_field
            True

        ::

            sage: P.<x,y> = ProjectiveSpace(CC, 1)
            sage: f = DynamicalSystem_projective([x^2 - CC.0/3*y^2, y^2])
            sage: f.dynatomic_polynomial(2)                                             # needs sage.libs.pari
            (x^4*y + (-0.666666666666667*I)*x^2*y^3 - x*y^4
             + (-0.111111111111111 - 0.333333333333333*I)*y^5)/(x^2*y - x*y^2
                                                                 + (-0.333333333333333*I)*y^3)

        ::

            sage: P.<x,y> = ProjectiveSpace(CC, 1)
            sage: f = DynamicalSystem_projective([x^2 - CC.0/5*y^2, y^2])
            sage: f.dynatomic_polynomial(2)                                             # needs sage.libs.pari
            x^2 + x*y + (1.00000000000000 - 0.200000000000000*I)*y^2

        ::

            sage: L.<t> = PolynomialRing(QuadraticField(2).maximal_order())             # needs sage.rings.number_field
            sage: P.<x, y> = ProjectiveSpace(L.fraction_field(), 1)
            sage: f = DynamicalSystem_projective([x^2 + (t^2 + 1)*y^2, y^2])
            sage: f.dynatomic_polynomial(2)                                             # needs sage.libs.pari sage.rings.number_field
            x^2 + x*y + (t^2 + 2)*y^2

        ::

            sage: P.<x,y> = ProjectiveSpace(ZZ, 1)
            sage: f = DynamicalSystem_projective([x^2 - 5*y^2, y^2])
            sage: f.dynatomic_polynomial([3,0 ])
            0

        TESTS:

        We check that the dynatomic polynomial has the right
        parent (see :issue:`18409`)::

            sage: # needs sage.rings.number_field
            sage: P.<x,y> = ProjectiveSpace(QQbar,1)
            sage: f = DynamicalSystem_projective([x^2 - 1/3*y^2, y^2])
            sage: f.dynatomic_polynomial(2).parent()                                    # needs sage.libs.pari
            Multivariate Polynomial Ring in x, y over Algebraic Field

        ::

            sage: # needs sage.rings.number_field
            sage: T.<v> = QuadraticField(33)
            sage: S.<t> = PolynomialRing(T)
            sage: P.<x,y> = ProjectiveSpace(FractionField(S),1)
            sage: f = DynamicalSystem_projective([t*x^2 - 1/t*y^2, y^2])
            sage: f.dynatomic_polynomial([1, 2]).parent()                               # needs sage.libs.pari
            Multivariate Polynomial Ring in x, y over Fraction Field of Univariate Polynomial
            Ring in t over Number Field in v with defining polynomial x^2 - 33 with v = 5.744562646538029?

        ::

            sage: P.<x, y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSystem_projective([x^3 - y^3*2, y^3])
            sage: f.dynatomic_polynomial(1).parent()
            Multivariate Polynomial Ring in x, y over Rational Field

        ::

            sage: R.<c> = QQ[]
            sage: P.<x,y> = ProjectiveSpace(R,1)
            sage: f = DynamicalSystem_projective([x^2 + c*y^2, y^2])
            sage: f.dynatomic_polynomial([1,2]).parent()                                # needs sage.libs.pari
            Multivariate Polynomial Ring in x, y over Univariate
            Polynomial Ring in c over Rational Field

        ::

            sage: R.<c> = QQ[]
            sage: P.<x,y> = ProjectiveSpace(ZZ,1)
            sage: f = DynamicalSystem_projective([x^2 + y^2, (1)*y^2 + (1)*x*y])
            sage: f.dynatomic_polynomial([1,2]).parent()                                # needs sage.libs.pari
            Multivariate Polynomial Ring in x, y over Integer Ring

        ::

            sage: P.<x, y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSystem_projective([x^2 + y^2, y^2])
            sage: f.dynatomic_polynomial(0)
            0
            sage: f.dynatomic_polynomial([0,0])
            0
            sage: f.dynatomic_polynomial(-1)
            Traceback (most recent call last):
            ...
            TypeError: period must be a positive integer

        ::

            sage: R.<c> = QQ[]
            sage: P.<x,y> = ProjectiveSpace(R,1)
            sage: f = DynamicalSystem_projective([x^2 + c*y^2, y^2])
            sage: f.dynatomic_polynomial([1,2]).parent()                                # needs sage.libs.pari
            Multivariate Polynomial Ring in x, y over Univariate Polynomial Ring in
            c over Rational Field

        Some rings still return :class:`SymoblicRing` elements::

            sage: S.<t> = FunctionField(CC)
            sage: P.<x,y> = ProjectiveSpace(S,1)
            sage: f = DynamicalSystem_projective([t*x^2-1*y^2, t*y^2])
            sage: f.dynatomic_polynomial([1, 2]).parent()                               # needs sage.libs.pari
            Symbolic Ring

        ::

            sage: R.<x,y> = PolynomialRing(QQ)
            sage: S = R.quo(R.ideal(y^2 - x + 1))
            sage: P.<u,v> = ProjectiveSpace(FractionField(S),1)
            sage: f = DynamicalSystem_projective([u^2 + S(x^2)*v^2, v^2])
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
        n = int(period[1])
        if n < 0:
            raise TypeError("period must be a positive integer")
        if n == 0:
            return self[0].parent().zero()
        if m == 0 and n == 1:
            return y*F0 - x*F1
        for d in range(1, n):
            if n % d == 0:
                PHI = PHI * ((y*F0 - x*F1)**moebius(n//d))
            F0, F1 = f0(F0, F1), f1(F0, F1)
        PHI = PHI * (y*F0 - x*F1)
        if m != 0:
            fm = self.nth_iterate_map(m)
            fm1 = self.nth_iterate_map(m - 1)
        try:
            QR = PHI.numerator().quo_rem(PHI.denominator())
            if not QR[1]:
                PHI = QR[0]
            if m != 0:
                PHI = PHI(fm._polys)/(PHI(fm1._polys))
                QR = PHI.numerator().quo_rem(PHI.denominator())
                if QR[1] == 0:
                    PHI = QR[0]
            return PHI
        except (TypeError, NotImplementedError): # something Singular can't handle
            if m != 0:
                PHI = PHI(fm._polys) / PHI(fm1._polys)
        #even when the ring can be passed to singular in quo_rem,
        #it can't always do the division, so we call Maxima
        if period != [0,1]: #period==[0,1] we don't need to do any division
            BR = self.domain().base_ring().base_ring()
            if not isinstance(BR, (sage.rings.abc.pAdicRing, sage.rings.abc.pAdicField)):
                try:
                    QR2 = PHI.numerator()._maxima_().divide(PHI.denominator())
                    if not QR2[1].sage():
                        # do it again to divide out by denominators of coefficients
                        PHI = QR2[0].sage()
                        PHI = PHI.numerator()._maxima_().divide(PHI.denominator())[0].sage()
                    if not isinstance(PHI, FractionFieldElement):
                        from sage.symbolic.expression_conversions import polynomial
                        PHI = polynomial(PHI, ring=self.coordinate_ring())
                except (TypeError, NotImplementedError): #something Maxima, or the conversion, can't handle
                    pass
        return PHI

    def nth_iterate_map(self, n, normalize=False):
        r"""
        Return the ``n``-th iterate of this dynamical system.

        ALGORITHM:

        Uses a form of successive squaring to reducing computations.

        .. TODO:: This could be improved.

        INPUT:

        - ``n`` -- positive integer

        - ``normalize`` -- boolean; remove gcd's during iteration

        OUTPUT: a projective dynamical system

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: f = DynamicalSystem_projective([x^2 + y^2, y^2])
            sage: f.nth_iterate_map(2)
            Dynamical System of Projective Space of dimension 1 over Rational Field
              Defn: Defined on coordinates by sending (x : y) to
                    (x^4 + 2*x^2*y^2 + 2*y^4 : y^4)

        ::

            sage: P.<x,y> = ProjectiveSpace(CC,1)
            sage: f = DynamicalSystem_projective([x^2 - y^2, x*y])
            sage: f.nth_iterate_map(3)
            Dynamical System of Projective Space of dimension 1
             over Complex Field with 53 bits of precision
              Defn: Defined on coordinates by sending (x : y) to
                    (x^8 + (-7.00000000000000)*x^6*y^2 + 13.0000000000000*x^4*y^4
                       + (-7.00000000000000)*x^2*y^6 + y^8
                     : x^7*y + (-4.00000000000000)*x^5*y^3 + 4.00000000000000*x^3*y^5 - x*y^7)

        ::

            sage: P.<x,y,z> = ProjectiveSpace(ZZ,2)
            sage: f = DynamicalSystem_projective([x^2 - y^2, x*y, z^2 + x^2])
            sage: f.nth_iterate_map(2)
            Dynamical System of Projective Space of dimension 2 over Integer Ring
              Defn: Defined on coordinates by sending (x : y : z) to
                    (x^4 - 3*x^2*y^2 + y^4 : x^3*y - x*y^3
                     : 2*x^4 - 2*x^2*y^2 + y^4 + 2*x^2*z^2 + z^4)

        ::

            sage: P.<x,y,z> = ProjectiveSpace(QQ,2)
            sage: X = P.subscheme(x*z-y^2)
            sage: f = DynamicalSystem_projective([x^2, x*z, z^2], domain=X)
            sage: f.nth_iterate_map(2)                                                  # needs sage.rings.function_field
            Dynamical System of Closed subscheme of Projective Space of dimension
            2 over Rational Field defined by:
              -y^2 + x*z
              Defn: Defined on coordinates by sending (x : y : z) to
                    (x^4 : x^2*z^2 : z^4)

        ::

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: f = DynamicalSystem_projective([y^2 * z^3, y^3 * z^2, x^5])
            sage: f.nth_iterate_map( 5, normalize=True)
            Dynamical System of Projective Space of dimension 2 over Rational Field
            Defn: Defined on coordinates by sending (x : y : z) to
            (y^202*z^443 : x^140*y^163*z^342 : x^645)
        """
        D = int(n)
        if D < 0:
            raise TypeError("iterate number must be a positive integer")
        if D == 1:
            return self
        H = End(self.domain())
        N = self.codomain().ambient_space().dimension_relative() + 1
        F = copy(self)
        Coord_ring = self.codomain().coordinate_ring()
        if isinstance(Coord_ring, QuotientRing_generic):
            PHI = H([Coord_ring.gen(i).lift() for i in range(N)])#makes a mapping
        else:
            PHI = H([Coord_ring.gen(i) for i in range(N)])
        while D:
            if D & 1:
                PHI = PHI*F
                if normalize:
                    PHI.normalize_coordinates()
            if D > 1: #avoid extra iterate
                F = F*F
            if normalize:
                F.normalize_coordinates()
            D >>= 1
        return PHI.as_dynamical_system()

    def nth_iterate(self, P, n, **kwds):
        r"""
        Return the ``n``-th iterate of the point ``P`` by this
        dynamical system.

        If ``normalize`` is ``True``, then the coordinates are
        automatically normalized.

        .. TODO:: Is there a more efficient way to do this?

        INPUT:

        - ``P`` -- a point in this map's domain

        - ``n`` -- positive integer

        kwds:

        - ``normalize`` -- boolean (default: ``False``)

        OUTPUT: a point in this map's codomain

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(ZZ,1)
            sage: f = DynamicalSystem_projective([x^2 + y^2, 2*y^2])
            sage: Q = P(1,1)
            sage: f.nth_iterate(Q,4)
            (32768 : 32768)

        ::

            sage: P.<x,y> = ProjectiveSpace(ZZ,1)
            sage: f = DynamicalSystem_projective([x^2 + y^2, 2*y^2])
            sage: Q = P(1,1)
            sage: f.nth_iterate(Q, 4, normalize=True)
            (1 : 1)

        ::

            sage: P.<x,y,z> = ProjectiveSpace(QQ,2)
            sage: f = DynamicalSystem_projective([x^2, 2*y^2, z^2 - x^2])
            sage: Q = P(2,7,1)
            sage: f.nth_iterate(Q,2)
            (-16/7 : -2744 : 1)

        ::

            sage: R.<t> = PolynomialRing(QQ)
            sage: P.<x,y,z> = ProjectiveSpace(R,2)
            sage: f = DynamicalSystem_projective([x^2 + t*y^2, (2-t)*y^2, z^2])
            sage: Q = P(2 + t, 7, t)
            sage: f.nth_iterate(Q,2)
            (t^4 + 2507*t^3 - 6787*t^2 + 10028*t + 16
             : -2401*t^3 + 14406*t^2 - 28812*t + 19208 : t^4)

        ::

            sage: P.<x,y,z> = ProjectiveSpace(ZZ,2)
            sage: X = P.subscheme(x^2-y^2)
            sage: f = DynamicalSystem_projective([x^2, y^2, z^2], domain=X)
            sage: f.nth_iterate(X(2,2,3), 3)
            (256 : 256 : 6561)

        ::

            sage: K.<c> = FunctionField(QQ)
            sage: P.<x,y> = ProjectiveSpace(K,1)
            sage: f = DynamicalSystem_projective([x^3 - 2*x*y^2 - c*y^3, x*y^2])
            sage: f.nth_iterate(P(c,1), 2)
            ((c^6 - 9*c^4 + 25*c^2 - c - 21)/(c^2 - 3) : 1)

            sage: P.<x,y,z> = ProjectiveSpace(QQ,2)
            sage: f = DynamicalSystem_projective([x^2 + 3*y^2, 2*y^2, z^2])
            sage: f.nth_iterate(P(2, 7, 1), -2)
            Traceback (most recent call last):
            ...
            TypeError: must be a forward orbit

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSystem_projective([x^3, x*y^2], domain=P)
            sage: f.nth_iterate(P(0, 1), 3, check=False)
            (0 : 0)
            sage: f.nth_iterate(P(0, 1), 3)
            Traceback (most recent call last):
            ...
            ValueError: [0, 0] does not define a valid projective point since all entries are zero

        ::

            sage: P.<x,y> = ProjectiveSpace(ZZ, 1)
            sage: f = DynamicalSystem_projective([x^3, x*y^2], domain=P)
            sage: f.nth_iterate(P(2,1), 3, normalize=False)
            (134217728 : 524288)
            sage: f.nth_iterate(P(2,1), 3, normalize=True)
            (256 : 1)

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: f = DynamicalSystem([x + y, y])
            sage: Q = (3,1)
            sage: f.nth_iterate(Q,0)
            (3 : 1)

        TESTS::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: f = DynamicalSystem([x^2 + y^2, y^2])
            sage: f.nth_iterate(0,0)
            (0 : 1)
        """
        n = Integer(n)
        if n < 0:
            raise TypeError("must be a forward orbit")
        return self.orbit(P, [n,n+1], **kwds)[0]

    def arakelov_zhang_pairing(self, g, **kwds):
        r"""
        Return an estimate of the Arakelov-Zhang pairing of the rational
        maps ``self`` and ``g`` on `\mathbb{P}^1` over a number field.

        The Arakelov-Zhang pairing was introduced by Petsche, Szpiro, and
        Tucker in 2012, which measures the dynamical closeness of two rational
        maps. They prove inter alia that if one takes a sequence of small points
        for one map (for example, preperiodic points for ``self``) and measure
        their dynamical height with respect to the other map (say, ``g``), then
        the values of the height will tend to the value of the Arakelov-Zhang pairing.

        The Arakelov-Zhang pairing involves mutual energy integrals between dynamical
        measures, which are in the case of polynomials, the equilibrium measures
        of the associated Julia sets at each place. As a result, these pairings
        are very difficult to compute exactly via analytic methods. We use a
        discrete approximation to these energy integrals.

        ALGORITHM:

        We select periodic points of order `n`, or ``n``-th preimages of a
        specified starting value given by ``f_starting_point`` and ``g_starting_point``.
        At the archimedean places and the places of bad reduction of the two maps,
        we compute the discrete approximations to the energy integrals involved
        using these points.

        INPUT:

        - ``g`` -- a rational map of `\mathbb{P}^1` given as a projective morphism
          ``g`` and ``self`` should have the same field of definition

        kwds:

        - ``n`` -- positive integer (default: 5); order of periodic points to
          use or preimages to take if starting points are specified

        - ``f_starting_point`` -- (default: ``None``) value in the base number field or None.
          If ``f_starting_point`` is None, we solve for points of period ``n`` for ``self``.
          Otherwise, we take ``n``-th preimages of the point given by ``f_starting_point``
          under ``f`` on the affine line.

        - ``g_starting_point`` -- (default: ``None``) value in the base number field or None.
          If ``g_starting_point`` is None, we solve for points of period ``n`` for ``g``.
          Otherwise, we take ``n``-th preimages of the point given by ``g_starting_point``
          under ``g`` on the affine line.

        - ``check_primes_of_bad_reduction`` -- boolean (default: ``False``);
          passed to the ``primes_of_bad_reduction`` function for ``self`` and ``g``

        - ``prec`` -- (default: ``RealField`` default);
          default precision for RealField values which are returned

        - ``noise_multiplier`` -- (default: 2) a real number.
          Discriminant terms involved in the computation at the archimedean places
          are often not needed, particularly if the capacity of the Julia sets is 1,
          and introduce a lot of error. By a well-known result of Mahler (see
          also M. Baker, ""A lower bound for averages of dynamical Green's
          functions") such error (for a set of `N` points) is on the order of
          `\log(N)/N` after our normalization. We check if the value of the
          archimedean discriminant terms is within ``2*noise_multiplier`` of
          `\log(N)/N`. If so, we discard it. In practice this greatly improves
          the accuracy of the estimate of the pairing. If desired,
          ``noise_multiplier`` can be set to 0, and no terms will be ignored.

        OUTPUT: a real number estimating the Arakelov-Zhang pairing of the two rational maps

        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: K.<k> = CyclotomicField(3)
            sage: P.<x,y> = ProjectiveSpace(K, 1)
            sage: f = DynamicalSystem_projective([x^2 + (2*k + 2)*y^2, y^2])
            sage: g = DynamicalSystem_projective([x^2, y^2])
            sage: pairingval = f.arakelov_zhang_pairing(g, n=5); pairingval
            0.409598197761958

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSystem_projective([x^2 + 4*y^2, y^2])
            sage: g = DynamicalSystem_projective([x^2, y^2])
            sage: pairingval = f.arakelov_zhang_pairing(g, n=6); pairingval             # needs sage.rings.function_field
            0.750178391443644
            sage: # Compare to the exact value:
            sage: dynheight = f.canonical_height(P(0, 1)); dynheight                    # needs sage.libs.pari
            0.75017839144364417318023000563
            sage: dynheight - pairingval                                                # needs sage.libs.pari sage.rings.function_field
            0.000000000000000

        Notice that if we set the noise_multiplier to 0, the accuracy is diminished::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSystem_projective([x^2 + 4*y^2, y^2])
            sage: g = DynamicalSystem_projective([x^2, y^2])
            sage: pairingval = f.arakelov_zhang_pairing(g, n=6, noise_multiplier=0)     # needs sage.rings.function_field
            sage: pairingval                                                            # needs sage.rings.number_field
            0.650660018921632
            sage: dynheight = f.canonical_height(P(0, 1)); dynheight                    # needs sage.libs.pari
            0.75017839144364417318023000563
            sage: pairingval - dynheight                                                # needs sage.libs.pari sage.rings.function_field
            -0.0995183725220122

        We compute the example of Prop. 18(d) from Petsche, Szpiro and Tucker::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSystem_projective([y^2 - (y - x)^2, y^2])
            sage: g = DynamicalSystem_projective([x^2, y^2])
            sage: f.arakelov_zhang_pairing(g)                                           # needs sage.rings.function_field
            0.326954667248466
            sage: # Correct value should be = 0.323067...
            sage: f.arakelov_zhang_pairing(g, n=9)      # long time                     # needs sage.rings.function_field
            0.323091061918965
            sage: _ - 0.323067                          # long time                     # needs sage.rings.function_field
            0.0000240619189654789

        Also from Prop. 18 of Petsche, Szpiro and Tucker, includes places of bad reduction::

            sage: # needs sage.rings.number_field
            sage: R.<z> = PolynomialRing(ZZ)
            sage: K.<b> = NumberField(z^3 - 11)
            sage: P.<x,y> = ProjectiveSpace(K,1)
            sage: a = 7/(b - 1)
            sage: f = DynamicalSystem_projective([a*y^2 - (a*y - x)^2, y^2])
            sage: g = DynamicalSystem_projective([x^2, y^2])

        If all archimedean absolute values of a have modulus > 2,
        then the pairing should be h(a).::

            sage: f.arakelov_zhang_pairing(g, n=6)      # long time                     # needs sage.rings.number_field
            1.93846423207664
            sage: _ - a.global_height()                 # long time                     # needs sage.rings.number_field
            -0.00744591697867292
        """
        n = kwds.pop('n', 5)
        f_starting_point = kwds.pop('f_starting_point', None)
        g_starting_point = kwds.pop('g_starting_point', None)
        check_primes_of_bad_reduction = kwds.pop('check_primes_of_bad_reduction', False)
        prec = kwds.pop('prec', None)
        noise_multiplier = kwds.pop('noise_multiplier', 2)

        f_domain = self.domain()
        R = f_domain.base_ring()
        g_domain = g.domain()

        if f_domain != g_domain:
            raise TypeError("Implemented only for rational maps of the same projective line.")

        if n <= 0:
            raise ValueError("Period must be a positive integer.")

        if not (isinstance(f_domain, ProjectiveSpace_ring) and isinstance(g_domain, ProjectiveSpace_ring)):
            raise NotImplementedError("Not implemented for subschemes.")

        if f_domain.dimension_relative() > 1:
            raise NotImplementedError("Only implemented for dimension 1.")

        if not self.is_endomorphism():
            raise TypeError("Self must be an endomorphism.")

        if R not in NumberFields() and not isinstance(R, sage.rings.abc.AlgebraicField):
            raise NotImplementedError("Only implemented for number fields.")

        f_iterate_map = self.nth_iterate_map(n)
        f_iter_map_poly = f_iterate_map.defining_polynomials()
        if f_starting_point is None:
            f_poly_hom = f_iter_map_poly[0] * f_domain.gens()[1] - f_iter_map_poly[1] * f_domain.gens()[0]
        else:
            f_poly_hom = f_iter_map_poly[0] - f_starting_point * f_iter_map_poly[1]

        g_iterate_map = g.nth_iterate_map(n)
        g_iter_map_poly = g_iterate_map.defining_polynomials()
        if g_starting_point is None:
            g_poly_hom = g_iter_map_poly[0] * g_domain.gens()[1] - g_iter_map_poly[1] * g_domain.gens()[0]
        else:
            g_poly_hom = g_iter_map_poly[0] - g_starting_point * g_iter_map_poly[1]

        f_poly = f_poly_hom([(f_domain.gens()[0]), 1]).univariate_polynomial().monic()
        g_poly = g_poly_hom([(g_domain.gens()[0]), 1]).univariate_polynomial().monic()

        # If f_poly and g_poly are not square-free, make them square-free.
        if not f_poly.is_squarefree():
            f_poly = f_poly.quo_rem(gcd(f_poly, f_poly.derivative()))[0]
        if not g_poly.is_squarefree():
            g_poly = g_poly.quo_rem(gcd(g_poly, g_poly.derivative()))[0]

        if f_poly.degree() <= 2 or g_poly.degree() <= 2:
            # f_point or g_point is exceptional
            raise ValueError("One of the starting points is exceptional. \
                              Please specify a non-exceptional initial point.")

        if gcd(f_poly, g_poly).degree() > 0:
            if f_poly.degree() > g_poly.degree():
                f_poly = f_poly.quo_rem(gcd(f_poly, g_poly))[0]
            else:
                g_poly = g_poly.quo_rem(gcd(f_poly, g_poly))[0]

            if f_poly.degree() <= 2 or g_poly.degree() <= 2:
                raise ValueError("After removing common factors, the n-th \
                                  iterates of 'self' and 'g' have too many \
                                  roots in common. Try another 'n' or starting \
                                  values.")

        # We want higher precision here temporarily, since resultants are
        # usually very large. This is not to say that the computation is
        # very accurate, merely that we want to keep track of potentially
        # very large height integers/rationals.
        old_prec = prec
        if prec is None:
            Real = RealField(512)
        elif prec < 512:
            prec = 512
            Real = RealField(prec)

        bad_primes = list(set(self.primes_of_bad_reduction(check=check_primes_of_bad_reduction))
                    .union(g.primes_of_bad_reduction(check=check_primes_of_bad_reduction)))

        f_deg = f_poly.degree()
        g_deg = g_poly.degree()

        f_disc = f_poly.discriminant()
        g_disc = g_poly.discriminant()

        res = f_poly.resultant(g_poly)

        # The code below actually computes -( mu_f - mu_g, mu_f - mu_g ),
        # so flip the sign at the end.
        AZ_pairing = Real(0)
        if R is QQ:
            for p in bad_primes:
                temp = (ZZ(1)/2) * (-f_disc.ord(p)) * Real(p).log() / (f_deg**2)
                if abs(temp) > noise_multiplier * Real(f_deg).log() / Real(f_deg):
                    AZ_pairing += temp

                temp = (ZZ(1)/2) * (-g_disc.ord(p)) * Real(p).log() / (g_deg**2)
                if abs(temp) > noise_multiplier * Real(g_deg).log() / Real(g_deg):
                    AZ_pairing += temp

                AZ_pairing -= (-res.ord(p)) * Real(p).log() / (f_deg * g_deg)

            temp = (ZZ(1)/2) * (Real(f_disc).abs().log()) / (f_deg**2)
            if abs(temp) > noise_multiplier * Real(f_deg).log() / Real(f_deg):
                AZ_pairing += temp

            temp = (ZZ(1)/2) * (Real(g_disc).abs().log()) / (g_deg**2)
            if abs(temp) > noise_multiplier * Real(g_deg).log() / Real(g_deg):
                AZ_pairing += temp

            AZ_pairing -= Real(res).abs().log() / (f_deg * g_deg)

        # For number fields
        else:
            K = self.base_ring()
            d = K.absolute_degree()

            for v in bad_primes:
                Nv = v.absolute_ramification_index() * v.residue_class_degree() / d

                temp = Nv * ((ZZ(1)/2) * K(f_disc).abs_non_arch(v, prec=prec).log() / (f_deg**2))
                if abs(temp) > noise_multiplier * Real(f_deg).log() / Real(f_deg):
                    AZ_pairing += temp

                temp = Nv * ((ZZ(1)/2) * K(g_disc).abs_non_arch(v, prec=prec).log() / (g_deg**2))
                if abs(temp) > noise_multiplier * Real(g_deg).log() / Real(g_deg):
                    AZ_pairing += temp

                AZ_pairing -= Nv * (K(res).abs_non_arch(v, prec=prec).log() / (f_deg * g_deg))

            if f_disc.is_rational():
                f_disc = QQ(f_disc)
                temp = (ZZ(1)/2) * (Real(f_disc).abs().log()) / (f_deg**2)
                if abs(temp) > noise_multiplier * Real(f_deg).log() / Real(f_deg):
                    AZ_pairing += temp
            else:
                temp = (ZZ(1)/d) * (ZZ(1)/2) * (Real(K(f_disc).norm()).abs().log()) / (f_deg**2)
                if abs(temp) > noise_multiplier * Real(f_deg).log() / Real(f_deg):
                    AZ_pairing += temp

            if g_disc.is_rational():
                g_disc = QQ(g_disc)
                temp = (ZZ(1)/2) * (Real(g_disc).abs().log()) / (g_deg**2)
                if abs(temp) > noise_multiplier * Real(g_deg).log() / Real(g_deg):
                    AZ_pairing += temp
            else:
                temp = (ZZ(1)/d) * (ZZ(1)/2) * (Real(K(g_disc).norm()).abs().log()) / (g_deg**2)
                if abs(temp) > noise_multiplier * Real(g_deg).log() / Real(g_deg):
                    AZ_pairing += temp

            if res.is_rational():
                AZ_pairing -= (Real(res).abs().log()) / (f_deg * g_deg)
            else:
                AZ_pairing -= (ZZ(1)/d) * (Real(K(res).norm()).abs().log()) / (f_deg * g_deg)

        if old_prec is None:
            Real = RealField()
        else:
            Real = RealField(old_prec)

        return Real(-AZ_pairing)

    def degree_sequence(self, iterates=2):
        r"""
        Return sequence of degrees of normalized iterates starting with
        the degree of this dynamical system.

        INPUT:

        - ``iterates`` -- (default: 2) positive integer

        OUTPUT: list of integers

        EXAMPLES::

            sage: P2.<X,Y,Z> = ProjectiveSpace(QQ, 2)
            sage: f = DynamicalSystem_projective([Z^2, X*Y, Y^2])
            sage: f.degree_sequence(15)                                                 # needs sage.rings.function_field
            [2, 3, 5, 8, 11, 17, 24, 31, 45, 56, 68, 91, 93, 184, 275]

        ::

            sage: F.<t> = PolynomialRing(QQ)
            sage: P2.<X,Y,Z> = ProjectiveSpace(F, 2)
            sage: f = DynamicalSystem_projective([Y*Z, X*Y, Y^2 + t*X*Z])
            sage: f.degree_sequence(5)                                                  # needs sage.rings.function_field
            [2, 3, 5, 8, 13]

        ::

            sage: P2.<X,Y,Z> = ProjectiveSpace(QQ, 2)
            sage: f = DynamicalSystem_projective([X^2, Y^2, Z^2])
            sage: f.degree_sequence(10)                                                 # needs sage.rings.function_field
            [2, 4, 8, 16, 32, 64, 128, 256, 512, 1024]

        ::

            sage: P2.<X,Y,Z> = ProjectiveSpace(ZZ, 2)
            sage: f = DynamicalSystem_projective([X*Y, Y*Z+Z^2, Z^2])
            sage: f.degree_sequence(10)                                                 # needs sage.rings.function_field
            [2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
        """
        if int(iterates) < 1:
            raise TypeError("number of iterates must be a positive integer")

        if self.is_morphism():
            d = self.degree()
            D = [d**t for t in range(1, iterates+1)]
        else:
            F = self
            F.normalize_coordinates()
            D = [F.degree()]
            for n in range(2, iterates+1):
                F = F*self
                F.normalize_coordinates()
                D.append(F.degree())
        return D

    def dynamical_degree(self, N=3, prec=53):
        r"""
        Return an approximation to the dynamical degree of this dynamical
        system. The dynamical degree is defined as
        `\lim_{n \to \infty} \sqrt[n]{\deg(f^n)}`.

        INPUT:

        - ``N`` -- (default: 3) positive integer, iterate to use
          for approximation

        - ``prec`` -- (default: 53) positive integer, real precision
          to use when computing root

        OUTPUT: real number

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSystem_projective([x^2 + x*y, y^2])
            sage: f.dynamical_degree()                                                  # needs sage.rings.function_field
            2.00000000000000

        ::

            sage: P2.<X,Y,Z> = ProjectiveSpace(ZZ, 2)
            sage: f = DynamicalSystem_projective([X*Y, Y*Z + Z^2, Z^2])
            sage: f.dynamical_degree(N=5, prec=100)                                     # needs sage.rings.function_field
            1.4309690811052555010452244131
        """
        if int(N) < 1:
            raise TypeError("number of iterates must be a positive integer")

        R = RealField(prec=prec)
        if self.is_morphism():
            return R(self.degree())
        else:
            D = self.nth_iterate_map(N, normalize=True).degree()
            return R(D).nth_root(N)

    def orbit(self, P, N, **kwds):
        r"""
        Return the orbit of the point ``P`` by this dynamical system.

        Let `F` be this dynamical system. If ``N`` is an integer return
        `[P,F(P),\ldots,F^N(P)]`. If ``N`` is a list or tuple `N=[m,k]`
        return `[F^m(P),\ldots,F^k(P)]`.
        Automatically normalize the points if ``normalize=True``. Perform
        the checks on point initialization if ``check=True``.

        INPUT:

        - ``P`` -- a point in this dynamical system's domain

        - ``n`` -- nonnegative integer or list or tuple of two
          nonnegative integers

        kwds:

        - ``check`` -- boolean (default: ``True``)

        - ``normalize`` -- boolean (default: ``False``)

        OUTPUT: list of points in this dynamical system's codomain

        EXAMPLES::

            sage: P.<x,y,z> = ProjectiveSpace(ZZ,2)
            sage: f = DynamicalSystem_projective([x^2 + y^2, y^2 - z^2, 2*z^2])
            sage: f.orbit(P(1,2,1), 3)
            [(1 : 2 : 1), (5 : 3 : 2), (34 : 5 : 8), (1181 : -39 : 128)]

        ::

            sage: P.<x,y,z> = ProjectiveSpace(ZZ,2)
            sage: f = DynamicalSystem_projective([x^2 + y^2, y^2 - z^2, 2*z^2])
            sage: f.orbit(P(1,2,1), [2,4])
            [(34 : 5 : 8), (1181 : -39 : 128), (1396282 : -14863 : 32768)]

        ::

            sage: P.<x,y,z> = ProjectiveSpace(ZZ,2)
            sage: X = P.subscheme(x^2 - y^2)
            sage: f = DynamicalSystem_projective([x^2, y^2, x*z], domain=X)
            sage: f.orbit(X(2,2,3), 3, normalize=True)
            [(2 : 2 : 3), (2 : 2 : 3), (2 : 2 : 3), (2 : 2 : 3)]

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: f = DynamicalSystem_projective([x^2 + y^2, y^2])
            sage: f.orbit(P.point([1,2], False), 4, check=False)
            [(1 : 2), (5 : 4), (41 : 16), (1937 : 256), (3817505 : 65536)]

        ::

            sage: K.<c> = FunctionField(QQ)
            sage: P.<x,y> = ProjectiveSpace(K,1)
            sage: f = DynamicalSystem_projective([x^2 + c*y^2, y^2])
            sage: f.orbit(P(0,1), 3)
            [(0 : 1), (c : 1), (c^2 + c : 1), (c^4 + 2*c^3 + c^2 + c : 1)]

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: f = DynamicalSystem_projective([x^2 + y^2, y^2], domain=P)
            sage: f.orbit(P.point([1, 2], False), 4, check=False)
            [(1 : 2), (5 : 4), (41 : 16), (1937 : 256), (3817505 : 65536)]

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: f = DynamicalSystem_projective([x^2, 2*y^2], domain=P)
            sage: f.orbit(P(2, 1),[-1, 4])
            Traceback (most recent call last):
            ...
            TypeError: orbit bounds must be nonnegative
            sage: f.orbit(P(2, 1), 0.1)
            Traceback (most recent call last):
            ...
            TypeError: Attempt to coerce non-integral RealNumber to Integer

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: f = DynamicalSystem_projective([x^3, x*y^2], domain=P)
            sage: f.orbit(P(0, 1), 3)
            Traceback (most recent call last):
            ...
            ValueError: [0, 0] does not define a valid projective point since all entries are zero
            sage: f.orbit(P(0, 1), 3, check=False)
            [(0 : 1), (0 : 0), (0 : 0), (0 : 0)]

        ::

            sage: P.<x,y> = ProjectiveSpace(ZZ, 1)
            sage: f = DynamicalSystem_projective([x^3, x*y^2], domain=P)
            sage: f.orbit(P(2,1), 3, normalize=False)
            [(2 : 1), (8 : 2), (512 : 32), (134217728 : 524288)]
            sage: f.orbit(P(2, 1), 3, normalize=True)
            [(2 : 1), (4 : 1), (16 : 1), (256 : 1)]

        ::

            sage: P.<x,y,z> = ProjectiveSpace(QQ,2)
            sage: f = DynamicalSystem_projective([x^2, y^2, x*z])
            sage: f.orbit((2/3, 1/3), 3)
            [(2/3 : 1/3 : 1), (2/3 : 1/6 : 1), (2/3 : 1/24 : 1), (2/3 : 1/384 : 1)]

        TESTS::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: f = DynamicalSystem([x^2 + y^2, y^2])
            sage: f.orbit(0, 0)
            [(0 : 1)]

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: f = DynamicalSystem([x^2 - y^2, y^2])
            sage: f.orbit(0,2)
            [(0 : 1), (-1 : 1), (0 : 1)]
        """
        if not isinstance(N,(list,tuple)):
            N = [0,N]
        N[0] = Integer(N[0])
        N[1] = Integer(N[1])
        if N[0] < 0 or N[1] < 0:
            raise TypeError("orbit bounds must be nonnegative")
        if N[0] > N[1]:
            return []

        R = self.domain()(P)
        if R in self.domain():  # Check whether R is a zero-dimensional point
            Q = R
        else:
            Q = P
        check = kwds.pop("check",True)
        normalize = kwds.pop("normalize",False)
        if normalize:
            Q.normalize_coordinates()
        for i in range(1, N[0]+1):
            Q = self(Q, check)
            if normalize:
                Q.normalize_coordinates()
        orb = [Q]
        for i in range(N[0]+1, N[1]+1):
            Q = self(Q, check)
            if normalize:
                Q.normalize_coordinates()
            orb.append(Q)
        return orb

    def resultant(self, normalize=False):
        r"""
        Compute the resultant of the defining polynomials of
        this dynamical system.

        If ``normalize`` is ``True``, then first normalize the coordinate
        functions with :meth:`normalize_coordinates`.

        INPUT:

        - ``normalize`` -- boolean (default: ``False``)

        OUTPUT: an element of the base ring of this map

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: f = DynamicalSystem_projective([x^2 + y^2, 6*y^2])
            sage: f.resultant()                                                         # needs sage.libs.pari
            36

        ::

            sage: R.<t> = PolynomialRing(GF(17))
            sage: P.<x,y> = ProjectiveSpace(R,1)
            sage: f = DynamicalSystem_projective([t*x^2 + t*y^2, 6*y^2])
            sage: f.resultant()                                                         # needs sage.libs.pari
            2*t^2

        ::

            sage: R.<t> = PolynomialRing(GF(17))
            sage: P.<x,y,z> = ProjectiveSpace(R,2)
            sage: f = DynamicalSystem_projective([t*x^2 + t*y^2, 6*y^2, 2*t*z^2])
            sage: f.resultant()
            13*t^8

        ::

            sage: P.<x,y,z> = ProjectiveSpace(QQ,2)
            sage: F = DynamicalSystem_projective([x^2 + y^2, 6*y^2, 10*x*z + z^2 + y^2])
            sage: F.resultant()
            1296

        ::

            sage: # needs sage.rings.number_field
            sage: R.<t> = PolynomialRing(QQ)
            sage: s = (t^3 + t + 1).roots(QQbar)[0][0]
            sage: P.<x,y> = ProjectiveSpace(QQbar, 1)
            sage: f = DynamicalSystem_projective([s*x^3 - 13*y^3, y^3 - 15*y^3])
            sage: f.resultant()
            871.6925062959149?
        """
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
                res = (f.lc() ** (d - g.degree()) * g.lc() ** (d - f.degree())
                       * f.__pari__().polresultant(g, x))
                return self.domain().base_ring()(res)
            except (TypeError, PariError):
                pass
        #Otherwise, use Macaulay
        R = F[0].parent()
        res = R.macaulay_resultant(list(F._polys))
        return res #Coercion here is not necessary as it is already done in Macaulay Resultant

    @cached_method
    def primes_of_bad_reduction(self, check=True):
        r"""
        Determine the primes of bad reduction for this dynamical system.

        Must be defined over a number field.

        If ``check`` is ``True``, each prime is verified to be of
        bad reduction.

        ALGORITHM:

        `p` is a prime of bad reduction if and only if the defining
        polynomials of ``self`` have a common zero. Or stated another way,
        `p` is a prime of bad reduction if and only if the radical of
        the ideal defined by the defining polynomials of ``self`` is not
        `(x_0,x_1,\ldots,x_N)`.  This happens if and only if some
        power of each `x_i` is not in the ideal defined by the
        defining polynomials of ``self``. This last condition is what is
        checked. The lcm of the coefficients of the monomials `x_i` in
        a Groebner basis is computed. This may return extra primes.

        INPUT:

        - ``check`` -- boolean (default: ``True``)

        OUTPUT: list of primes

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: f = DynamicalSystem_projective([1/3*x^2 + 1/2*y^2, y^2])
            sage: f.primes_of_bad_reduction()                                           # needs sage.rings.function_field
            [2, 3]

        ::

            sage: P.<x,y,z,w> = ProjectiveSpace(QQ,3)
            sage: f = DynamicalSystem_projective([12*x*z - 7*y^2, 31*x^2 - y^2,
            ....:                                 26*z^2, 3*w^2 - z*w])
            sage: f.primes_of_bad_reduction()                                           # needs sage.rings.function_field
            [2, 3, 7, 13, 31]

        A number field example::

            sage: # needs sage.rings.number_field
            sage: R.<z> = QQ[]
            sage: K.<a> = NumberField(z^2 - 2)
            sage: P.<x,y> = ProjectiveSpace(K,1)
            sage: f = DynamicalSystem_projective([1/3*x^2+1/a*y^2, y^2])
            sage: f.primes_of_bad_reduction()                                           # needs sage.rings.function_field
            [Fractional ideal (a), Fractional ideal (3)]

        This is an example where ``check=False`` returns extra primes::

            sage: P.<x,y,z> = ProjectiveSpace(ZZ,2)
            sage: f = DynamicalSystem_projective([3*x*y^2 + 7*y^3 - 4*y^2*z + 5*z^3,
            ....:                                 -5*x^3 + x^2*y + y^3 + 2*x^2*z,
            ....:                                 -2*x^2*y + x*y^2 + y^3 - 4*y^2*z + x*z^2])
            sage: f.primes_of_bad_reduction(False)                                      # needs sage.rings.function_field
            [2, 5, 37, 2239, 304432717]
            sage: f.primes_of_bad_reduction()                                           # needs sage.rings.function_field
            [5, 37, 2239, 304432717]
        """
        if (not isinstance(self.domain(), ProjectiveSpace_ring)) or (not isinstance(self.codomain(), ProjectiveSpace_ring)):
            raise NotImplementedError("not implemented for subschemes")
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
                return badprimes
        else:
            raise TypeError("base ring must be number field or number field ring")

    def conjugate(self, M, adjugate=False, normalize=False):
        r"""
        Conjugate this dynamical system by ``M``, i.e. `M^{-1} \circ f \circ M`.

        If possible the new map will be defined over the same space.
        Otherwise, will try to coerce to the base ring of ``M``.

        INPUT:

        - ``M`` -- a square invertible matrix

        - ``adjugate`` -- boolean (default: ``False``); also classically called
          adjoint, takes a square matrix ``M`` and finds the transpose of its
          cofactor matrix. Used for conjugation in place of inverse when
          specified ``True``. Functionality is the same in projective space.

        - ``normalize`` -- boolean (default: ``False``); if ``normalize`` is
          ``True``, then the method ``normalize_coordinates`` is called

        OUTPUT: a dynamical system

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(ZZ,1)
            sage: f = DynamicalSystem_projective([x^2 + y^2, y^2])
            sage: f.conjugate(matrix([[1,2], [0,1]]))
            Dynamical System of Projective Space of dimension 1 over Integer Ring
              Defn: Defined on coordinates by sending (x : y) to
                    (x^2 + 4*x*y + 3*y^2 : y^2)

        ::

            sage: R.<x> = PolynomialRing(QQ)
            sage: K.<i> = NumberField(x^2 + 1)                                          # needs sage.rings.number_field
            sage: P.<x,y> = ProjectiveSpace(ZZ,1)
            sage: f = DynamicalSystem_projective([x^3 + y^3, y^3])
            sage: f.conjugate(matrix([[i,0], [0,-i]]))                                  # needs sage.rings.number_field
            Dynamical System of Projective Space of dimension 1 over Integer Ring
              Defn: Defined on coordinates by sending (x : y) to
                    (-x^3 + y^3 : -y^3)

        ::

            sage: P.<x,y,z> = ProjectiveSpace(ZZ,2)
            sage: f = DynamicalSystem_projective([x^2 + y^2, y^2, y*z])
            sage: f.conjugate(matrix([[1,2,3], [0,1,2], [0,0,1]]))
            Dynamical System of Projective Space of dimension 2 over Integer Ring
              Defn: Defined on coordinates by sending (x : y : z) to
                    (x^2 + 4*x*y + 3*y^2 + 6*x*z + 9*y*z + 7*z^2 : y^2 + 2*y*z : y*z + 2*z^2)

        ::

            sage: P.<x,y> = ProjectiveSpace(ZZ,1)
            sage: f = DynamicalSystem_projective([x^2+y^2, y^2])
            sage: f.conjugate(matrix([[2,0], [0,1/2]]))
            Dynamical System of Projective Space of dimension 1 over Rational Field
              Defn: Defined on coordinates by sending (x : y) to
                    (2*x^2 + 1/8*y^2 : 1/2*y^2)

        ::

            sage: R.<x> = PolynomialRing(QQ)
            sage: K.<i> = NumberField(x^2 + 1)                                          # needs sage.rings.number_field
            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: f = DynamicalSystem_projective([1/3*x^2 + 1/2*y^2, y^2])
            sage: f.conjugate(matrix([[i,0], [0,-i]]))                                  # needs sage.rings.number_field
            Dynamical System of Projective Space of dimension 1
             over Number Field in i with defining polynomial x^2 + 1
              Defn: Defined on coordinates by sending (x : y) to
                    ((1/3*i)*x^2 + (1/2*i)*y^2 : (-i)*y^2)

        TESTS::

            sage: R = ZZ
            sage: P.<x,y> = ProjectiveSpace(R,1)
            sage: f = DynamicalSystem_projective([x^2 + y^2, y^2])
            sage: m = matrix(R, 2, [4, 3, 2, 1])
            sage: f.conjugate(m, normalize=False)
            Dynamical System of Projective Space of dimension 1 over Rational Field
              Defn: Defined on coordinates by sending (x : y) to
                    (-4*x^2 - 8*x*y - 7/2*y^2 : 12*x^2 + 20*x*y + 8*y^2)
            sage: f.conjugate(m, adjugate=True)
            Dynamical System of Projective Space of dimension 1 over Integer Ring
              Defn: Defined on coordinates by sending (x : y) to
                    (8*x^2 + 16*x*y + 7*y^2 : -24*x^2 - 40*x*y - 16*y^2)

        .. TODO::

            Use the left and right action functionality to replace the code below with
            #return DynamicalSystem_projective(M.inverse()*self*M, domain=self.codomain())
            once there is a function to pass to the smallest field of definition.
        """
        if not (M.is_square() == 1 and M.determinant() != 0
            and M.ncols() == self.domain().ambient_space().dimension_relative() + 1):
            raise TypeError("matrix must be invertible and size dimension + 1")
        X = M * vector(self[0].parent().gens())
        F = vector(self._polys)
        F = F(list(X))
        if adjugate:
            N = M.adjugate()
        else:
            N = M.inverse()
        F = N * F
        R = self.codomain().ambient_space().coordinate_ring()
        try:
            F = [R(f) for f in F]
            PS = self.codomain()
        except TypeError: #no longer defined over same ring
            R = R.change_ring(N.base_ring())
            F = [R(f) for f in F]
            PS = self.codomain().change_ring(N.base_ring())
        G = DynamicalSystem_projective(F, domain=PS)
        if normalize:
            G.normalize_coordinates()
        return G

    def green_function(self, P, v, **kwds):
        r"""
        Evaluate the local Green's function at the place ``v`` for ``P``
        with ``N`` terms of the series or to within a given error bound.

        Must be over a number field or order of a number field. Note that
        this is the absolute local Green's function so is scaled by the
        degree of the base field.

        Use ``v=0`` for the archimedean place over `\QQ` or field embedding.
        Non-archimedean places are prime ideals for number fields or primes
        over `\QQ`.

        ALGORITHM:

        See Exercise 5.29 and Figure 5.6 of [Sil2007]_.

        INPUT:

        - ``P`` -- a projective point

        - ``v`` -- nonnegative integer; a place, use ``0`` for the
          archimedean place

        kwds:

        - ``N`` -- (default: 10) positive integer; number of
          terms of the series to use

        - ``prec`` -- (default: 100) positive integer, float point or
          `p`-adic precision

        - ``error_bound`` -- (optional) a positive real number

        OUTPUT: a real number

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: f = DynamicalSystem_projective([x^2 + y^2, x*y]);
            sage: Q = P(5, 1)
            sage: f.green_function(Q, 0, N=30)
            1.6460930159932946233759277576

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: f = DynamicalSystem_projective([x^2 + y^2, x*y]);
            sage: Q = P(5, 1)
            sage: f.green_function(Q, 0, N=200, prec=200)
            1.6460930160038721802875250367738355497198064992657997569827

        ::

            sage: # needs sage.rings.number_field
            sage: K.<w> = QuadraticField(3)
            sage: P.<x,y> = ProjectiveSpace(K,1)
            sage: f = DynamicalSystem_projective([17*x^2 + 1/7*y^2, 17*w*x*y])
            sage: f.green_function(P.point([w, 2], False), K.places()[1])
            1.7236334013785676107373093775
            sage: f.green_function(P([2, 1]), K.ideal(7), N=7)
            0.48647753726382832627633818586
            sage: f.green_function(P([w, 1]), K.ideal(17), error_bound=0.001)
            -0.70821687320448199545278619351

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: f = DynamicalSystem_projective([x^2 + y^2, x*y])
            sage: f.green_function(P.point([5,2], False), 0, N=30)
            1.7315451844777407992085512000
            sage: f.green_function(P.point([2,1], False), 0, N=30)
            0.86577259223181088325226209926
            sage: f.green_function(P.point([1,1], False), 0, N=30)
            0.43288629610862338612700146098
        """
        N = kwds.get('N', 10)                       #Get number of iterates (if entered)
        err = kwds.get('error_bound', None)         #Get error bound (if entered)
        prec = kwds.get('prec', 100)                #Get precision (if entered)
        R = RealField(prec)
        localht = R(0)
        BR = FractionField(P.codomain().base_ring())
        GBR = self.change_ring(BR) #so the heights work

        if BR not in NumberFields():
            raise NotImplementedError("must be over a number field or a number field order")
        if not BR.is_absolute():
            raise TypeError("must be an absolute field")

        #For QQ the 'flip-trick' works better over RR or Qp
        if isinstance(v, (NumberFieldFractionalIdeal, RingHomomorphism_im_gens)):
            K = BR
        elif is_prime(v):
            K = Qp(v, prec)
        elif v == 0 and BR == QQ:
            K = R
            v = BR.places(prec=prec)[0]
        else:
            raise ValueError("invalid valuation (=%s) entered" % v)

        #Coerce all polynomials in F into polynomials with coefficients in K
        F = self.change_ring(K, check=False)
        d = F.degree()
        dim = F.codomain().ambient_space().dimension_relative()
        Q = P.change_ring(K, check=False)

        if err is not None:
            err = R(err)
            if not err > 0:
                raise ValueError("error bound (=%s) must be positive" % err)

            #if doing error estimates, compute needed number of iterates
            D = (dim + 1) * (d - 1) + 1
            # compute upper bound
            if isinstance(v, RingHomomorphism_im_gens): #archimedean
                vindex = BR.places(prec=prec).index(v)
                emb = BR.places(prec=prec)[vindex]
                U = GBR.local_height_arch(vindex, prec=prec) + R(binomial(dim + d, d)).log()
            else: #non-archimedean
                U = GBR.local_height(v, prec=prec)

            #compute lower bound - from explicit polynomials of Nullstellensatz
            CR = GBR.codomain().ambient_space().coordinate_ring() #.lift() only works over fields
            I = CR.ideal(GBR.defining_polynomials())
            maxh = 0
            for k in range(dim + 1):
                CoeffPolys = (CR.gen(k) ** D).lift(I)
                h = 1
                for poly in CoeffPolys:
                    if poly != 0:
                        if isinstance(v, RingHomomorphism_im_gens): #archimedean
                            if BR == QQ:
                                h = max([R(K(c).abs()) for c in poly.coefficients()])
                            else:
                                h = max([R(emb(c).abs()) for c in poly.coefficients()])
                        else: #non-archimedean
                            if BR == QQ:
                                h = max([R(v)**(-R(c.valuation(v)))  for c in poly.coefficients()])
                            else:
                                h = max([R(c.abs_non_arch(v, prec=prec)) for c in poly.coefficients()])
                        maxh = max(h, maxh)
            if maxh == 0:
                maxh = 1  #avoid division by 0
            if isinstance(v, RingHomomorphism_im_gens): #archimedean
                L = R(1 / ((dim + 1) * binomial(dim + D - d, D - d) * maxh)).log().abs()
            else: #non-archimedean
                if BR == QQ:
                    L = ((-self.resultant().valuation(v))*R(v).log()).abs()
                else:
                    L = (self.resultant().abs_non_arch(v, prec=prec)).log().abs()
            C = max([U, L])
            if C != 0:
                N = R(C / (err*(d-1))).log(d).abs().ceil()
            else: #we just need log||P||_v
                N = 1

        #START GREEN FUNCTION CALCULATION
        if isinstance(v, RingHomomorphism_im_gens):  #embedding for archimedean local height
            for i in range(N+1):
                Qv = [ (v(t).abs()) for t in Q ]
                m = -1
                #compute the maximum absolute value of entries of a, and where it occurs
                for n in range(dim + 1):
                    if Qv[n] > m:
                        j = n
                        m = Qv[n]
                # add to sum for the Green's function
                localht += ((1/R(d))**R(i)) * (R(m).log())
                #get the next iterate
                if i < N:
                    Q.scale_by(1/Q[j])
                    Q = F(Q, False)
            return (1/BR.absolute_degree()) * localht

        #else - prime or prime ideal for non-archimedean
        for i in range(N + 1):
            if BR == QQ:
                Qv = [ R(K(t).abs()) for t in Q ]
            else:
                Qv = [ R(t.abs_non_arch(v)) for t in Q ]
            m = -1
            #compute the maximum absolute value of entries of a, and where it occurs
            for n in range(dim + 1):
                if Qv[n] > m:
                    j = n
                    m = Qv[n]
            # add to sum for the Green's function
            localht += (1/R(d))**R(i) * (R(m).log())
            #get the next iterate
            if i < N:
                Q.scale_by(1 / Q[j])
                Q = F(Q, False)
        return (1 / BR.absolute_degree()) * localht

    def canonical_height(self, P, **kwds):
        r"""
        Evaluate the (absolute) canonical height of ``P`` with respect to
        this dynamical system.

        Must be over number field or order of a number field. Specify
        either the number of terms of the series to evaluate or the
        error bound required.

        ALGORITHM:

        The sum of the Green's function at the archimedean places and
        the places of bad reduction.

        If function is defined over `\QQ` uses Wells' Algorithm, which
        allows us to not have to factor the resultant.

        INPUT:

        - ``P`` -- a projective point

        kwds:

        - ``badprimes`` -- (optional) a list of primes of bad reduction

        - ``N`` -- (default: 10) positive integer; number of
          terms of the series to use in the local green functions

        - ``prec`` -- (default: 100) positive integer, float point or
          `p`-adic precision

        - ``error_bound`` -- (optional) a positive real number

        OUTPUT: a real number

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(ZZ,1)
            sage: f = DynamicalSystem_projective([x^2 + y^2, 2*x*y]);
            sage: f.canonical_height(P.point([5,4]), error_bound=0.001)                 # needs sage.libs.pari
            2.1970553519503404898926835324
            sage: f.canonical_height(P.point([2,1]), error_bound=0.001)                 # needs sage.libs.pari
            1.0984430632822307984974382955

        Notice that preperiodic points may not return exactly 0::

            sage: # needs sage.rings.number_field
            sage: R.<X> = PolynomialRing(QQ)
            sage: K.<a> = NumberField(X^2 + X - 1)
            sage: P.<x,y> = ProjectiveSpace(K,1)
            sage: f = DynamicalSystem_projective([x^2 - 2*y^2, y^2])
            sage: Q = P.point([a,1])
            sage: f.canonical_height(Q, error_bound=0.000001)  # Answer only within error_bound of 0
            5.7364919788790160119266380480e-8
            sage: f.nth_iterate(Q, 2) == Q  # but it is indeed preperiodic
            True

        ::

            sage: P.<x,y,z> = ProjectiveSpace(QQ,2)
            sage: X = P.subscheme(x^2 - y^2);
            sage: f = DynamicalSystem_projective([x^2, y^2, 4*z^2], domain=X);
            sage: Q = X([4,4,1])
            sage: f.canonical_height(Q, badprimes=[2])                                  # needs sage.rings.function_field
            0.0013538030870311431824555314882

        ::

            sage: P.<x,y,z> = ProjectiveSpace(QQ,2)
            sage: X = P.subscheme(x^2 - y^2);
            sage: f = DynamicalSystem_projective([x^2, y^2, 30*z^2], domain=X)
            sage: Q = X([4, 4, 1])
            sage: f.canonical_height(Q, badprimes=[2,3,5], prec=200)                    # needs sage.rings.function_field
            2.7054056208276961889784303469356774912979228770208655455481

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSystem_projective([1000*x^2 - 29*y^2, 1000*y^2])
            sage: Q = P(-1/4, 1)
            sage: f.canonical_height(Q, error_bound=0.01)                               # needs sage.libs.pari
            3.7979215342343045582800170705

        ::

            sage: RSA768 = 123018668453011775513049495838496272077285356959533479219732245215\
            ....: 1726400507263657518745202199786469389956474942774063845925192557326303453731548\
            ....: 2685079170261221429134616704292143116022212404792747377940806653514195974598569\
            ....: 02143413
            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: f = DynamicalSystem_projective([RSA768*x^2 + y^2, x*y])
            sage: Q = P(RSA768,1)
            sage: f.canonical_height(Q, error_bound=0.00000000000000001)                # needs sage.libs.pari
            931.18256422718241278672729195

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSystem([2*(-2*x^3 + 3*(x^2*y)) + 3*y^3, 3*y^3])
            sage: f.canonical_height(P(1,0))                                            # needs sage.libs.pari
            0.00000000000000000000000000000
        """
        bad_primes = kwds.get("badprimes", None)
        prec = kwds.get("prec", 100)
        error_bound = kwds.get("error_bound", None)
        K = FractionField(self.codomain().base_ring())

        if K not in NumberFields():
            if not isinstance(K, sage.rings.abc.AlgebraicField):
                raise NotImplementedError("must be over a number field or a number field order or QQbar")
            else:
                #since this an absolute height, we can compute the height of a QQbar point
                #by choosing any number field it is defined over.
                Q = P._number_field_from_algebraics()
                K = Q.codomain().base_ring()
                f = self._number_field_from_algebraics().as_dynamical_system()
                if K == QQ:
                    K = f.base_ring()
                    Q = Q.change_ring(K)
                elif f.base_ring() == QQ:
                    f = f.change_ring(K)
                else:
                    K, phi, psi, b = K.composite_fields(f.base_ring(), both_maps=True)[0]
                    Q = Q.change_ring(phi)
                    f = f.change_ring(psi)
        else:
            if not K.is_absolute():
                raise TypeError("must be an absolute field")
            Q = P
            f = self

        # After moving from QQbar to K being something like QQ, we need
        # to renormalize f, especially to match the normalized resultant.
        f.normalize_coordinates()

        # If our map and point are defined on P^1(QQ), use Wells' Algorithm
        # instead of the usual algorithm using local Green's functions:
        if K is QQ and self.codomain().ambient_space().dimension_relative() == 1:
            # write our point with coordinates whose gcd is 1
            Q.normalize_coordinates()
            if Q.parent().value_ring() is QQ:
                Q.clear_denominators()
            # assures integer coefficients
            coeffs = f[0].coefficients() + f[1].coefficients()
            t = 1
            for c in coeffs:
                t = lcm(t, c.denominator())
            A = t * f[0]
            B = t * f[1]
            Res = f.resultant(normalize=True).abs()
            H = 0
            x_i = Q[0]
            y_i = Q[1]
            d = self.degree()
            R = RealField(prec)
            N = kwds.get('N', 10)
            err = kwds.get('error_bound', None)
            # computes the error bound as defined in Algorithm 3.1 of [WELLS]
            if Res > 1:
                if err is not None:
                    err = err / 2
                    N = ceil((R(Res).log().log() - R(d-1).log() - R(err).log())/(R(d).log()))
                    N = max(N, 1)
                    kwds.update({'error_bound': err})
                    kwds.update({'N': N})
                for n in range(N):
                    x = A(x_i,y_i) % Res**(N-n)
                    y = B(x_i,y_i) % Res**(N-n)
                    g = gcd([x, y, Res])
                    H = H + R(g).abs().log() / (d**(n+1))
                    x_i = x / g
                    y_i = y / g
            # this looks different than Wells' Algorithm because of the difference
            # between what Wells' calls H_infty,
            # and what Green's Function returns for the infinite place
            h = f.green_function(Q, 0, **kwds) - H + R(t).log()
            # The value returned by Well's algorithm may be negative. As the canonical height
            # is always nonnegative, so if this value is within -err of 0, return 0.
            if h < 0:
                assert h > -err, "A negative height less than -error_bound was computed. " + \
                 "This should be impossible, please report bug on https://github.com/sagemath/sage/issues"
                    # This should be impossible. The error bound for Wells' is rigorous
                    # and the actual height is always >= 0. If we see something less than -err,
                    # something has g one very wrong.
                h = R(0)
            return h

        if bad_primes is None:
            bad_primes = []
            for b in Q:
                if K == QQ:
                    bad_primes += b.denominator().prime_factors()
                else:
                    bad_primes += b.denominator_ideal().prime_factors()
            bad_primes += K(f.resultant(normalize=True)).support()
            bad_primes = list(set(bad_primes))

        emb = K.places(prec=prec)
        num_places = len(emb) + len(bad_primes)
        if error_bound is not None:
            error_bound /= num_places
        R = RealField(prec)
        h = R.zero()

        ##update the keyword dictionary for use in green_function
        kwds.update({"badprimes": bad_primes})
        kwds.update({"error_bound": error_bound})

        # Archimedean local heights
        # :: WARNING: If places is fed the default Sage precision of 53 bits,
        # it uses Real or Complex Double Field in place of RealField(prec) or ComplexField(prec).
        # RDF is an instance of a separate class.
        for v in emb:
            if isinstance(v.codomain(), (sage.rings.abc.RealField, sage.rings.abc.RealDoubleField)):
                dv = R.one()
            else:
                dv = R(2)
            h += dv * f.green_function(Q, v, **kwds)     #arch Green function

        # Non-Archimedean local heights
        for v in bad_primes:
            if K == QQ:
                dv = R.one()
            else:
                dv = R(v.residue_class_degree() * v.absolute_ramification_index())
            h += dv * f.green_function(Q, v, **kwds)  #non-arch Green functions
        return h

    def height_difference_bound(self, prec=None):
        r"""
        Return an upper bound on the different between the canonical
        height of a point with respect to this dynamical system and the
        absolute height of the point.

        This map must be a morphism.

        ALGORITHM:

        Uses a Nullstellensatz argument to compute the constant.
        For details: see [Hutz2015]_.

        INPUT:

        - ``prec`` -- (default: :class:`RealField` default)
          positive integer, float point precision

        OUTPUT: a real number

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSystem_projective([x^2 + y^2, x*y])
            sage: f.height_difference_bound()                                           # needs sage.symbolic
            1.38629436111989

            sage: P.<x,y,z> = ProjectiveSpace(ZZ, 2)
            sage: f = DynamicalSystem_projective([4*x^2 + 100*y^2, 210*x*y, 10000*z^2])
            sage: f.height_difference_bound()                                           # needs sage.symbolic
            10.3089526606443

        A number field example::

            sage: # needs sage.rings.number_field
            sage: R.<x> = QQ[]
            sage: K.<c> = NumberField(x^3 - 2)
            sage: P.<x,y,z> = ProjectiveSpace(K, 2)
            sage: f = DynamicalSystem_projective([1/(c+1)*x^2 + c*y^2, 210*x*y, 10000*z^2])
            sage: f.height_difference_bound()                                           # needs sage.symbolic
            11.3683039374269

        ::

            sage: # needs sage.rings.number_field sage.symbolic
            sage: P.<x,y,z> = ProjectiveSpace(QQbar, 2)
            sage: f = DynamicalSystem_projective([x^2, QQbar(sqrt(-1))*y^2,
            ....:                                 QQbar(sqrt(3))*z^2])
            sage: f.height_difference_bound()
            2.89037175789616

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSystem([5*x^2 + 3*x*y , y^2 + 3*x^2])
            sage: f.height_difference_bound(prec=100)                                   # needs sage.symbolic
            5.3375380797013179737224159274
        """
        FF = FractionField(self.domain().base_ring()) #lift will only work over fields, so coercing into FF
        if FF not in NumberFields():
            if isinstance(FF, sage.rings.abc.AlgebraicField):
                #since this is absolute height, we can choose any number field over which the
                #function is defined.
                f = self._number_field_from_algebraics()
            else:
                raise NotImplementedError("fraction field of the base ring must be a number field or QQbar")
        else:
            f = self.change_ring(FF)
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
        maxh = 0
        for k in range(N + 1):
            CoeffPolys = (CR.gen(k) ** D).lift(I)
            h = max([g.global_height(prec) for g in CoeffPolys])
            maxh = max(maxh, h)
        L = R((N + 1) * binomial(N + D - d, D - d)).log() + maxh
        C = max(U, L) #height difference dh(P) - L <= h(f(P)) <= dh(P) +U
        return C / (d - 1)

    def multiplier(self, P, n, check=True):
        r"""
        Return the multiplier of the point ``P`` of period ``n`` with
        respect to this dynamical system.

        INPUT:

        - ``P`` -- a point on domain of this map

        - ``n`` -- positive integer, the period of ``P``

        - ``check`` -- boolean (default: ``True``); verify that ``P``
          has period ``n``

        OUTPUT:

        A square matrix of size ``self.codomain().dimension_relative()``
        in the ``base_ring`` of this dynamical system.

        EXAMPLES::

            sage: P.<x,y,z> = ProjectiveSpace(QQ,2)
            sage: f = DynamicalSystem_projective([x^2, y^2, 4*z^2])
            sage: Q = P.point([4,4,1], False)
            sage: f.multiplier(Q,1)
            [2 0]
            [0 2]

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: f = DynamicalSystem_projective([7*x^2 - 28*y^2, 24*x*y])
            sage: f.multiplier(P(2,5), 4)
            [231361/20736]

        ::

            sage: P.<x,y> = ProjectiveSpace(CC,1)
            sage: f = DynamicalSystem_projective([x^3 - 25*x*y^2 + 12*y^3, 12*y^3])
            sage: f.multiplier(P(1,1), 5)
            [0.389017489711934]

        ::

            sage: P.<x,y> = ProjectiveSpace(RR,1)
            sage: f = DynamicalSystem_projective([x^2 - 2*y^2, y^2])
            sage: f.multiplier(P(2,1), 1)
            [4.00000000000000]

        ::

            sage: P.<x,y> = ProjectiveSpace(Qp(13),1)                                   # needs sage.rings.padics
            sage: f = DynamicalSystem_projective([x^2 - 29/16*y^2, y^2])
            sage: f.multiplier(P(5,4), 3)                                               # needs sage.rings.padics
            [6 + 8*13 + 13^2 + 8*13^3 + 13^4 + 8*13^5 + 13^6 + 8*13^7 + 13^8 +
             8*13^9 + 13^10 + 8*13^11 + 13^12 + 8*13^13 + 13^14 + 8*13^15 + 13^16 +
             8*13^17 + 13^18 + 8*13^19 + O(13^20)]

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: f = DynamicalSystem_projective([x^2 - y^2, y^2])
            sage: f.multiplier(P(0,1), 1)
            Traceback (most recent call last):
            ...
            ValueError: (0 : 1) is not periodic of period 1
        """
        if check:
            if self.nth_iterate(P, n) != P:
                raise ValueError("%s is not periodic of period %s" % (P, n))
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
        for i in range(n):
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
        Return the multiplier of the point ``P`` of period ``n`` with
        respect to this dynamical system modulo `p^k`.

        This map must be an endomorphism of projective space defined
        over `\QQ` or `\ZZ`. This function should not be used at the top
        level as it does not perform input checks. It is used primarily
        for the rational preperiodic and periodic point algorithms.

        INPUT:

        - ``P`` -- a point on domain of this map

        - ``n`` -- positive integer, the period of ``P``

        - ``p`` -- positive integer

        - ``k`` -- positive integer

        OUTPUT:

        A square matrix of size ``self.codomain().dimension_relative()``
        in `\ZZ/(p^k)\ZZ`.

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: f = DynamicalSystem_projective([x^2 - 29/16*y^2, y^2])
            sage: f._multipliermod(P(5,4), 3, 11, 1)
            [3]

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: f = DynamicalSystem_projective([x^2 - 29/16*y^2, y^2])
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
        for i in range(n):
            F = []
            R = self(Q, False)
            g = gcd(R._coords)
            R.scale_by(1 / g)
            R_list = list(R)
            for index in range(N + 1):
                R_list[index] = R_list[index] % (p ** k)
            R._coords = tuple(R_list)
            index = N
            while R[index] % p == 0:
                index -= 1
            indexlist.append(index)
            #dehomogenize and compute multiplier
            F = self.dehomogenize((indexlist[i],indexlist[i+1]))
            l = (F.jacobian()(tuple(Q.dehomogenize(indexlist[i])))*l) % (p ** k)
            Q = R
        return l

    def _nth_preimage_tree_helper(self, Q, n, m, **kwds):
        r"""
        A recursive method to fill in ``n``-th preimage tree.

        This helper function is used by ``nth_preimage_tree`` below to actually compute the
        points of the tree and populate the dictionary used to create a ``DiGraph``
        object. Note the addition of an ``m`` parameter, which counts upwards as n counts
        downwards to keep track of what level we are at in the tree for the purposes of
        returning points and displaying the point's level in the tree.
        """
        return_points = kwds.get("return_points", False)
        numerical = kwds.get("numerical", False)
        prec = kwds.get("prec", 100)
        display_complex = kwds.get("display_complex", False)
        digits = kwds.get("digits", 5)
        embed = kwds.get("embed", None)
        D = {}
        if numerical:
            # Solve for preimages numerically
            CR = self.domain().ambient_space().coordinate_ring()
            fn = self.dehomogenize(1)
            poly = (fn[0].numerator()*CR(Q[1]) - fn[0].denominator()*CR(Q[0])).univariate_polynomial()
            K = ComplexField(prec=prec)
            pre = [ProjectiveSpace(K,1)(r) for r in poly.roots(ring=K)]
        else:
            # Solve for preimages algebraically
            pre = self.rational_preimages(Q,1)
        for pt in pre:
            # Fill in dictionary entries of preimage points to Q
            if display_complex:
                pt1 = "(" + str(embed(pt[0]).n(digits=digits)) + ": 1)"
                Q1 = "(" + str(embed(Q[0]).n(digits=digits)) + ": 1)"
                key = pt1 + ", " + str(m)
                D[key] = [Q1 + ", " + str(m-1)]
            else:
                key = str(pt) + ", " + str(m)
                D[key] = [str(Q) + ", " + str(m-1)]
            if return_points:
                # Fill in m-th level preimage points in points list
                kwds["points"][m].append(pt)

        if return_points:
            points = kwds["points"]
            if n == 1:
                # Base case of recursion
                return D, points
            else:
                # For each preimage point of Q, use recursion to find that point's preimages
                # and update the dictionary
                for pt in pre:
                    D.update(self._nth_preimage_tree_helper(pt, n-1, m+1, **kwds)[0])
            return D, points
        else:
            if n == 1:
                # Base case of recursion
                return D
            else:
                # For each preimage point of Q, use recursion to find that point's preimages
                # and update the dictionary
                for pt in pre:
                    D.update(self._nth_preimage_tree_helper(pt, n-1, m+1, **kwds))
            return D

    def nth_preimage_tree(self, Q, n, **kwds):
        r"""
        Return the ``n``-th pre-image tree rooted at ``Q``.

        This map must be an endomorphism of the projective line defined
        over a number field, algebraic field, or finite field.

        INPUT:

        - ``Q`` -- a point in the domain of this map

        - ``n`` -- positive integer, the depth of the pre-image tree

        kwds:

        - ``return_points`` -- boolean (default: ``False``); if ``True``,
          return a list of lists where the index `i` is the level of the tree
          and the elements of the list at that index are the `i`-th preimage
          points as an algebraic element of the splitting field of the
          polynomial `f^n - Q = 0`.

        - ``numerical`` -- boolean (default: ``False``); calculate pre-images
          numerically. Note if this is set to ``True``, preimage points are
          displayed as complex numbers.

        - ``prec`` -- (default: 100) positive integer; the precision of the
          ``ComplexField`` if we compute the preimage points numerically

        - ``display_labels`` -- boolean (default: ``True``); whether to display
          vertex labels. Since labels can be very cluttered, can set
          ``display_labels`` to ``False`` and use ``return_points`` to get a
          hold of the points themselves, either as algebraic or complex numbers.

        - ``display_complex`` -- boolean (default: ``False``); display vertex
          labels as complex numbers. Note if this option is chosen that we must
          choose an embedding from the splitting field ``field_def`` of the
          `n`-th-preimage equation into `\CC`. We make the choice of the first
          embedding returned by ``field_def.embeddings(ComplexField())``.

        - ``digits`` -- positive integer; the number of decimal digits to
          display for complex numbers. This only applies if ``display_complex``
          is set to ``True``.

        OUTPUT:

        If ``return_points`` is ``False``, a :class:`GraphPlot` object representing
        the `n`-th pre-image tree.  If ``return_points`` is ``True``, a tuple
        ``(GP, points)``, where ``GP`` is a :class:`GraphPlot` object, and
        ``points`` is a list of lists as described above under
        ``return_points``.

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: f = DynamicalSystem_projective([x^2 + y^2, y^2])
            sage: Q = P(0,1)
            sage: f.nth_preimage_tree(Q, 2)                                             # needs sage.plot
            GraphPlot object for Digraph on 7 vertices

        ::

            sage: P.<x,y> = ProjectiveSpace(GF(3), 1)
            sage: f = DynamicalSystem_projective([x^2 + x*y + y^2, y^2])
            sage: Q = P(0,1)
            sage: f.nth_preimage_tree(Q, 2, return_points=True)                         # needs sage.plot
            (GraphPlot object for Digraph on 4 vertices,
             [[(0 : 1)], [(1 : 1)], [(0 : 1), (2 : 1)]])
        """
        return_points = kwds.get("return_points", False)
        numerical = kwds.get("numerical", False)
        prec = kwds.get("prec", 100)
        display_labels = kwds.get("display_labels", True)
        display_complex = kwds.get("display_complex", False)
        digits = kwds.get("digits", 5)

        if self.domain().dimension_relative() > 1:
            raise NotImplementedError("only implemented for dimension 1")
        base_ring = self.base_ring()
        if isinstance(base_ring, sage.rings.abc.AlgebraicField):
            if numerical:
                raise ValueError("can't solve numerically over QQbar, no embedding into CC")
            fbar = self
            # No embedding from QQbar into C
            kwds["display_complex"] = False
            display_complex = False
        elif base_ring in NumberFields():
            if numerical:
                field_def = ComplexField(prec=prec)
                embed = base_ring.embeddings(field_def)[0]
                fbar = self.change_ring(embed)
                embed = End(field_def).identity()
                kwds["display_complex"] = True
                display_complex = True
                kwds["embed"] = embed
            else:
                field_def = self.field_of_definition_preimage(Q, n)
                fbar = self.change_ring(field_def)
                if display_complex:
                    embed = field_def.embeddings(ComplexField())[0]
                    kwds["embed"] = embed
        elif base_ring in FiniteFields():
            if numerical:
                raise ValueError("can't solve numerically over a finite field, no embedding into CC")
            field_def = self.field_of_definition_preimage(Q, n)
            fbar = self.change_ring(field_def)
            # No embedding from finite field into C
            kwds["display_complex"] = False
            display_complex = False
        else:
            raise NotImplementedError("only implemented for number fields, algebraic fields, and finite fields")

        Q = fbar.codomain()(Q)
        if return_points:
            # n+1 since we have n levels with root as 0th level
            points = [[] for i in range(n+1)]
            points[0].append(Q)
            kwds["points"] = points
            V, points = fbar._nth_preimage_tree_helper(Q, n, 1, **kwds)
        else:
            V = fbar._nth_preimage_tree_helper(Q, n, 1, **kwds)
        from sage.graphs.digraph import DiGraph
        from sage.graphs.graph_plot import GraphPlot
        G = DiGraph(V)
        if display_complex:
            Q = "(" + str(embed(Q[0]).n(digits=digits)) + ": 1)"
            root = Q + ", " + str(0)
        else:
            root = str(Q) + ", " + str(0)
        options = {'layout':'tree', 'tree_orientation':'up', 'tree_root':root, 'vertex_labels':display_labels}

        if return_points:
            return GraphPlot(G, options), points
        else:
            return GraphPlot(G, options)

    def possible_periods(self, **kwds):
        r"""
        Return the set of possible periods for rational periodic points of
        this dynamical system.

        Must be defined over `\ZZ` or `\QQ`.

        ALGORITHM:

        Calls ``self.possible_periods()`` modulo all primes of good reduction
        in range ``prime_bound``. Return the intersection of those lists.

        INPUT: keyword arguments:

        - ``prime_bound`` -- (default: ``[1, 20]``) a list or tuple of
          two positive integers or an integer for the upper bound

        - ``bad_primes`` -- (optional) a list or tuple of integer primes,
          the primes of bad reduction

        - ``ncpus`` -- (default: all cpus) number of cpus to use in parallel

        OUTPUT: list of positive integers

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: f = DynamicalSystem_projective([x^2 - 29/16*y^2, y^2])
            sage: f.possible_periods(ncpus=1)                                           # needs sage.rings.function_field
            [1, 3]

        ::

            sage: PS.<x,y> = ProjectiveSpace(1,QQ)
            sage: f = DynamicalSystem_projective([5*x^3 - 53*x*y^2 + 24*y^3, 24*y^3])
            sage: f.possible_periods(prime_bound=[1,5])                                 # needs sage.rings.function_field
            Traceback (most recent call last):
            ...
            ValueError: no primes of good reduction in that range
            sage: f.possible_periods(prime_bound=[1,10])                                # needs sage.rings.function_field
            [1, 4, 12]
            sage: f.possible_periods(prime_bound=[1,20])                                # needs sage.rings.function_field
            [1, 4]

        ::

            sage: P.<x,y,z> = ProjectiveSpace(ZZ,2)
            sage: f = DynamicalSystem_projective([2*x^3 - 50*x*z^2 + 24*z^3,
            ....:                                 5*y^3 - 53*y*z^2 + 24*z^3, 24*z^3])
            sage: f.possible_periods(prime_bound=10)                                    # needs sage.rings.function_field
            [1, 2, 6, 20, 42, 60, 140, 420]
            sage: f.possible_periods(prime_bound=20) # long time
            [1, 20]
        """
        if self.domain().base_ring() not in [ZZ, QQ]:
            raise NotImplementedError("must be ZZ or QQ")

        primebound = kwds.pop("prime_bound", [1, 20])
        badprimes = kwds.pop("bad_primes", None)
        num_cpus = kwds.pop("ncpus", ncpus())

        if not isinstance(primebound, (list, tuple)):
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
            if q not in badprimes:
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
            return sorted(periods)

    def _preperiodic_points_to_cyclegraph(self, preper):
        r"""
        Given the complete set of periodic or preperiodic points return the
        digraph representing the orbit.

        If ``preper`` is not the complete set, this function will not fill
        in the gaps.

        INPUT:

        - ``preper`` -- list or tuple of projective points; the complete
          set of rational periodic or preperiodic points

        OUTPUT:

        A digraph representing the orbit the rational preperiodic points
        ``preper`` in projective space.

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: f = DynamicalSystem_projective([x^2 - 2*y^2, y^2])
            sage: preper = [P(-2, 1), P(1, 0), P(0, 1), P(1, 1), P(2, 1), P(-1, 1)]
            sage: f._preperiodic_points_to_cyclegraph(preper)                           # needs sage.graphs
            Looped digraph on 6 vertices
        """
        V = []
        E = []
        #We store the points we encounter is a list, D. Each new point is checked to
        #see if it is in that list (which uses ==) so that equal points with different
        #representations only appear once in the graph.
        D = []
        for val in preper:
            try:
                V.append(D[D.index(val)])
            except ValueError:
                D.append(val)
                V.append(val)
            Q = self(val)
            Q.normalize_coordinates()
            try:
                E.append([D[D.index(Q)]])
            except ValueError:
                D.append(Q)
                E.append([Q])
        from sage.graphs.digraph import DiGraph
        g = DiGraph(dict(zip(V, E)), loops=True)
        return g

    def is_PGL_minimal(self, prime_list=None):
        r"""
        Check if this dynamical system is a minimal model in
        its conjugacy class.

        See [BM2012]_ and [Mol2015]_ for a description of the algorithm.
        For polynomial maps it uses [HS2018]_.

        INPUT:

        - ``prime_list`` -- (optional) list of primes to check minimality

        OUTPUT: boolean

        EXAMPLES::

            sage: PS.<X,Y> = ProjectiveSpace(QQ,1)
            sage: f = DynamicalSystem_projective([X^2 + 3*Y^2, X*Y])
            sage: f.is_PGL_minimal()                                                    # needs sage.rings.function_field
            True

        ::

            sage: PS.<x,y> = ProjectiveSpace(QQ,1)
            sage: f = DynamicalSystem_projective([6*x^2 + 12*x*y + 7*y^2, 12*x*y])
            sage: f.is_PGL_minimal()                                                    # needs sage.rings.function_field
            False

        ::

            sage: PS.<x,y> = ProjectiveSpace(QQ,1)
            sage: f = DynamicalSystem_projective([6*x^2 + 12*x*y + 7*y^2, y^2])
            sage: f.is_PGL_minimal()                                                    # needs sage.rings.function_field
            False
        """
        if self.base_ring() != QQ and self.base_ring() != ZZ:
            raise NotImplementedError("minimal models only implemented over ZZ or QQ")
        if not self.is_morphism():
            raise TypeError("the function is not a morphism")
        if self.degree() == 1:
            raise NotImplementedError("minimality is only for degree 2 or higher")

        f = copy(self)
        f.normalize_coordinates()
        R = f.domain().coordinate_ring()
        F = R(f[0].numerator())
        G = R(f[0].denominator())
        if G.degree() == 0 or F.degree() == 0:
            #can't use BM for polynomial
            from .endPN_minimal_model import HS_minimal
            g, m = HS_minimal(self, return_transformation=True, D=prime_list)
            return m == m.parent().one()

        from .endPN_minimal_model import affine_minimal
        return affine_minimal(self, return_transformation=False, D=prime_list, quick=True)

    def minimal_model(self, return_transformation=False, prime_list=None, algorithm=None, check_primes=True):
        r"""
        Determine if this dynamical system is minimal.

        This dynamical system must be defined over the projective line
        over the rationals. In particular, determine if this map is affine
        minimal, which is enough to decide if it is minimal or not.
        See Proposition 2.10 in [BM2012]_.

        INPUT:

        - ``return_transformation`` -- boolean (default: ``False``); this
          signals a return of the `PGL_2` transformation to conjugate
          this map to the calculated minimal model

        - ``prime_list`` -- (optional) a list of primes, in case one
          only wants to determine minimality at those specific primes

        - ``algorithm`` -- (optional) string; can be one of the following:

        - ``check_primes`` -- (optional) boolean; this signals whether to
            check whether each element in ``prime_list`` is a prime

          * ``'BM'`` -- the Bruin-Molnar algorithm [BM2012]_
          * ``'HS'`` -- the Hutz-Stoll algorithm [HS2018]_

        OUTPUT:

        - a dynamical system on the projective line which is a minimal model
          of this map

        - a `PGL(2,\QQ)` element which conjugates this map to a minimal model

        EXAMPLES::

            sage: PS.<X,Y> = ProjectiveSpace(QQ,1)
            sage: f = DynamicalSystem_projective([X^2 + 3*Y^2, X*Y])
            sage: f.minimal_model(return_transformation=True)                           # needs sage.rings.function_field
            (
            Dynamical System of Projective Space of dimension 1 over Rational
            Field
              Defn: Defined on coordinates by sending (X : Y) to
                    (X^2 + 3*Y^2 : X*Y)
            ,
            [1 0]
            [0 1]
            )

        ::

            sage: PS.<X,Y> = ProjectiveSpace(QQ,1)
            sage: f = DynamicalSystem_projective([7365/2*X^4 + 6282*X^3*Y + 4023*X^2*Y^2
            ....:                                   + 1146*X*Y^3 + 245/2*Y^4,
            ....:                                 -12329/2*X^4 - 10506*X^3*Y - 6723*X^2*Y^2
            ....:                                   - 1914*X*Y^3 - 409/2*Y^4])
            sage: f.minimal_model(return_transformation=True)                           # needs sage.rings.function_field
            (
            Dynamical System of Projective Space of dimension 1 over Rational Field
              Defn: Defined on coordinates by sending (X : Y) to
                    (9847*X^4 + 28088*X^3*Y + 30048*X^2*Y^2 + 14288*X*Y^3 + 2548*Y^4
                    : -12329*X^4 - 35164*X^3*Y - 37614*X^2*Y^2 - 17884*X*Y^3 - 3189*Y^4),
            <BLANKLINE>
            [2 1]
            [0 1]
            )

        ::

            sage: PS.<x,y> = ProjectiveSpace(QQ,1)
            sage: f = DynamicalSystem_projective([6*x^2 + 12*x*y + 7*y^2, 12*x*y])
            sage: f.minimal_model()                                                     # needs sage.rings.function_field
            Dynamical System of Projective Space of dimension 1 over Rational Field
              Defn: Defined on coordinates by sending (x : y) to
                    (x^2 + 12*x*y + 42*y^2 : 2*x*y)

        ::

            sage: PS.<x,y> = ProjectiveSpace(ZZ,1)
            sage: f = DynamicalSystem_projective([6*x^2 + 12*x*y + 7*y^2, 12*x*y + 42*y^2])
            sage: g,M = f.minimal_model(return_transformation=True, algorithm='BM')     # needs sage.rings.function_field
            sage: f.conjugate(M) == g                                                   # needs sage.rings.function_field
            True

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSystem([2*x^2, y^2])
            sage: f.minimal_model(return_transformation=True)                           # needs sage.rings.function_field
            (
            Dynamical System of Projective Space of dimension 1 over Rational Field
              Defn: Defined on coordinates by sending (x : y) to
                    (x^2 : y^2)                                                    ,
            [1 0]
            [0 2]
            )
            sage: f.minimal_model(prime_list=[3])                                       # needs sage.rings.function_field
            Dynamical System of Projective Space of dimension 1 over Rational Field
              Defn: Defined on coordinates by sending (x : y) to
                    (2*x^2 : y^2)

        TESTS::

            sage: PS.<X,Y> = ProjectiveSpace(QQ,1)
            sage: f = DynamicalSystem_projective([X + Y, X - 3*Y])
            sage: f.minimal_model()                                                     # needs sage.rings.function_field
            Traceback (most recent call last):
            ...
            NotImplementedError: minimality is only for degree 2 or higher

        ::

            sage: PS.<X,Y> = ProjectiveSpace(QQ,1)
            sage: f = DynamicalSystem_projective([X^2 - Y^2, X^2 + X*Y])
            sage: f.minimal_model()                                                     # needs sage.rings.function_field
            Traceback (most recent call last):
            ...
            TypeError: the function is not a morphism

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: f = DynamicalSystem([2*x^2, y^2])
            sage: f.minimal_model(algorithm='BM')                                       # needs sage.rings.function_field
            Traceback (most recent call last):
            ...
            TypeError: affine minimality is only considered for
            maps not of the form f or 1/f for a polynomial f

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: f = DynamicalSystem([2*x^2, y^2])
            sage: f.minimal_model(prime_list=[0])                                       # needs sage.rings.function_field
            Traceback (most recent call last):
            ...
            ValueError: prime_list contains 0 which is not prime

        REFERENCES:

        - [BM2012]_
        - [Mol2015]_
        - [HS2018]_
        """
        if self.base_ring() != ZZ and self.base_ring() != QQ:
            raise NotImplementedError("minimal models only implemented over ZZ or QQ")
        if not self.is_morphism():
            raise TypeError("the function is not a morphism")
        if self.degree() == 1:
            raise NotImplementedError("minimality is only for degree 2 or higher")
        if prime_list and check_primes:
            for p in prime_list:
                if not p.is_prime():
                    raise ValueError("prime_list contains " + str(p) + " which is not prime")

        if algorithm == 'BM':
            from .endPN_minimal_model import affine_minimal
            return affine_minimal(self, return_transformation=return_transformation, D=prime_list, quick=False)
        if algorithm == 'HS':
            from .endPN_minimal_model import HS_minimal
            return HS_minimal(self, return_transformation=return_transformation, D=prime_list)
        # algorithm not specified
        f = copy(self)
        f.normalize_coordinates()
        R = f.domain().coordinate_ring()
        F = R(f[0].numerator())
        G = R(f[0].denominator())

        if G.degree() == 0 or F.degree() == 0:
            #can use BM for polynomial
            from .endPN_minimal_model import HS_minimal
            return HS_minimal(self, return_transformation=return_transformation, D=prime_list)

        if prime_list is None:
            prime_list = ZZ(F.resultant().prime_divisors())
        if max(prime_list) > 500:
            from .endPN_minimal_model import affine_minimal
            return affine_minimal(self, return_transformation=return_transformation,
                                  D=prime_list, quick=False)

    def all_minimal_models(self, return_transformation=False, prime_list=None,
                           algorithm=None, check_minimal=True):
        r"""
        Determine a representative in each `SL(2,\ZZ)`-orbit of this map.

        This can be done either with the Bruin-Molnar algorithm or the
        Hutz-Stoll algorithm. The Hutz-Stoll algorithm requires the map
        to have minimal resultant and then finds representatives in orbits
        with minimal resultant. The Bruin-Molnar algorithm finds
        representatives with the same resultant (up to sign) of the given map.

        Bruin-Molnar does not work for polynomials and is more efficient
        for large primes.

        INPUT:

        - ``return_transformation`` -- boolean (default: ``False``); this
          signals a return of the `PGL_2` transformation to conjugate
          this map to the calculated models

        - ``prime_list`` -- (optional) a list of primes, in case one
          only wants to determine minimality at those specific primes

        - ``algorithm`` -- (optional) string; can be one of the following:

          * ``'BM'`` -- the Bruin-Molnar algorithm [BM2012]_
          * ``'HS'`` -- for the Hutz-Stoll algorithm [HS2018]_

          if not specified, properties of the map are utilized to choose

        - ``check_minimal`` -- (optional) boolean; to first check if the map
          is minimal and if not, compute a minimal model before computing
          for orbit representatives

        OUTPUT:

        A list of pairs `(F,m)`, where `F` is dynamical system on the
        projective line and `m` is the associated `PGL(2,\QQ)` element.
        Or just a list of dynamical systems if not returning the conjugation.

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSystem([2*x^2, 3*y^2])
            sage: f.all_minimal_models()                                                # needs sage.rings.function_field
            [Dynamical System of Projective Space of dimension 1 over Rational Field
               Defn: Defined on coordinates by sending (x : y) to
                     (x^2 : y^2)]

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: c = 2*3^6
            sage: f = DynamicalSystem([x^3 - c^2*y^3, x*y^2])
            sage: len(f.all_minimal_models(algorithm='HS'))                             # needs sage.rings.function_field
            14
            sage: len(f.all_minimal_models(prime_list=[2], algorithm='HS'))             # needs sage.rings.function_field
            2

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSystem([237568*x^3 + 1204224*x^2*y + 2032560*x*y^2
            ....:     + 1142289*y^3, -131072*x^3 - 663552*x^2*y - 1118464*x*y^2
            ....:     - 627664*y^3])
            sage: len(f.all_minimal_models(algorithm='BM'))                             # needs sage.rings.function_field
            2

        TESTS::

            sage: # needs sage.rings.function_field
            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: c = 2^2*5^2*11^3
            sage: f = DynamicalSystem([x^3 - c^2*y^3, x*y^2])
            sage: MM = f.all_minimal_models(return_transformation=True, algorithm='BM')
            sage: all(f.conjugate(m) == F for F, m in MM)
            True
            sage: MM = f.all_minimal_models(return_transformation=True, algorithm='HS')
            sage: all(f.conjugate(m) == F for F,m in MM)
            True

        REFERENCES:

        - [BM2012]_
        - [HS2018]_
        """
        if self.base_ring() != ZZ and self.base_ring() != QQ:
            raise NotImplementedError("minimal models only implemented over ZZ or QQ")
        if not self.is_morphism():
            raise TypeError("the function is not a morphism")
        if self.degree() == 1:
            raise NotImplementedError("minimality is only for degree 2 or higher")

        if check_minimal:
            f, m = self.minimal_model(return_transformation=True,
                                      prime_list=prime_list,
                                      algorithm=algorithm)
        else:
            f = self
            m = matrix(ZZ, 2, 2, [1,0,0,1])

        if algorithm == 'BM':
            from .endPN_minimal_model import BM_all_minimal
            models = BM_all_minimal(f, return_transformation=True, D=prime_list)
        elif algorithm == 'HS':
            from .endPN_minimal_model import HS_all_minimal
            models = HS_all_minimal(f, return_transformation=True, D=prime_list)
        else: # algorithm not specified
            f.normalize_coordinates()
            Aff_f = f.dehomogenize(1)
            R = Aff_f.domain().coordinate_ring()
            F = R(Aff_f[0].numerator())
            G = R(Aff_f[0].denominator())
            if G.degree() == 0 or F.degree() == 0:
                #can use BM for polynomial
                from .endPN_minimal_model import HS_all_minimal
                models = HS_all_minimal(f, return_transformation=True, D=prime_list)
            elif prime_list is None:
                prime_list = ZZ(f.resultant()).prime_divisors()
                if prime_list == []:
                    models = [[f,m]]
                elif max(prime_list) > 500:
                    from .endPN_minimal_model import BM_all_minimal
                    models = BM_all_minimal(f, return_transformation=True, D=prime_list)
                else:
                    from .endPN_minimal_model import HS_all_minimal
                    models = HS_all_minimal(f, return_transformation=True, D=prime_list)

        if return_transformation:
            models = [[g, t*m] for g,t in models]
        else:
            models = [g for g,t in models]
        return models

    def affine_preperiodic_model(self, m, n, return_conjugation=False):
        r"""
        Return a dynamical system conjugate to this one with affine (n, m) preperiodic points.

        If the base ring of this dynamical system is finite, there may
        not be a model with affine preperiodic points, in which case a
        :exc:`ValueError` is raised.

        INPUT:

        - ``m`` -- the preperiod of the preperiodic points to make affine

        - ``n`` -- the period of the preperiodic points to make affine

        - ``return_conjugation`` -- boolean (default: ``False``); if ``True``, return a tuple
          ``(g, phi)`` where ``g`` is a model with affine (n, m) preperiodic points
          and ``phi`` is the matrix that moves ``f`` to ``g``.

        OUTPUT: a dynamical system conjugate to this one

        EXAMPLES::

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: f = DynamicalSystem_projective([x^2, y^2, z^2])
            sage: g = f.affine_preperiodic_model(0, 1); g                               # needs sage.rings.function_field
            Dynamical System of Projective Space of dimension 2 over Rational Field
              Defn: Defined on coordinates by sending (x : y : z) to
                    (-x^2 : -2*x^2 + 2*x*y - y^2 : 2*x^2 - 2*x*y + 2*y^2 + 2*y*z + z^2)

        We can check that ``g`` has affine fixed points::

            sage: g.periodic_points(1)                                                  # needs sage.rings.function_field
            [(-1 : -1 : 1),
             (-1/2 : -1 : 1),
             (-1/2 : -1/2 : 1),
             (-1/3 : -2/3 : 1),
             (0 : -1 : 1),
             (0 : -1/2 : 1),
             (0 : 0 : 1)]

        ::

            sage: # needs sage.rings.finite_rings
            sage: P.<x,y,z> = ProjectiveSpace(GF(9), 2)
            sage: f = DynamicalSystem_projective([x^2, y^2, z^2])
            sage: f.affine_preperiodic_model(0, 1)                                      # needs sage.rings.function_field
            Dynamical System of Projective Space of dimension 2
             over Finite Field in z2 of size 3^2
              Defn: Defined on coordinates by sending (x : y : z) to
                    ((-z2)*x^2 : z2*x^2 + (-z2)*x*y + (-z2)*y^2 :
                     (-z2)*x^2 + z2*x*y + (z2 + 1)*y^2 - y*z + z^2)

        ::

            sage: R.<c> = GF(3)[]
            sage: P.<x,y,z> = ProjectiveSpace(R, 2)
            sage: f = DynamicalSystem_projective([x^2, y^2, z^2])
            sage: f.affine_preperiodic_model(0, 1)  # long time
            Dynamical System of Projective Space of dimension 2 over
             Univariate Polynomial Ring in c over Finite Field of size 3
              Defn: Defined on coordinates by sending (x : y : z) to
                    (2*c^3*x^2 : c^3*x^2 + 2*c^3*x*y + 2*c^3*y^2 :
                     2*c^3*x^2 + c^3*x*y + (c^3 + c^2)*y^2 + 2*c^2*y*z + c^2*z^2)

        ::

            sage: # needs sage.rings.number_field
            sage: K.<k> = CyclotomicField(3)
            sage: P.<x,y,z> = ProjectiveSpace(K, 2)
            sage: f = DynamicalSystem_projective([x^2 + k*x*y + y^2, z^2, y^2])
            sage: f.affine_preperiodic_model(1, 1)                                      # needs sage.rings.function_field
            Dynamical System of Projective Space of dimension 2
             over Cyclotomic Field of order 3 and degree 2
              Defn: Defined on coordinates by sending (x : y : z) to
                    (-y^2 : x^2 : x^2 + (-k)*x*z + z^2)

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSystem_projective([x^2 + y^2, y^2])
            sage: g, mat = f.affine_preperiodic_model(0, 1, return_conjugation=True)    # needs sage.rings.function_field
            sage: g == f.conjugate(mat)                                                 # needs sage.rings.function_field
            True

        ::

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: X = P.subscheme(2*y - z)
            sage: f = DynamicalSystem_projective([x^2 + y^2, z^2 + y^2, z^2], domain=X)
            sage: f.affine_preperiodic_model(0, 1)                                      # needs sage.rings.function_field
            Dynamical System of Closed subscheme of Projective Space of dimension 2
             over Rational Field defined by: 2*y - z
              Defn: Defined on coordinates by sending (x : y : z) to
                    (-x^2 - y^2 : y^2 : x^2 + z^2)

        TESTS::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSystem_projective([x^2 + 2*y^2, x^2])
            sage: g, mat = f.affine_preperiodic_model(0, 1, return_conjugation=True)    # needs sage.rings.function_field
            sage: f.conjugate(mat) == g                                                 # needs sage.rings.function_field
            True
        """
        n = ZZ(n)
        if n < 1:
            raise ValueError('period must be positive')
        m = ZZ(m)
        if m < 0:
            raise ValueError('preperiod must be nonnegative')
        f = self
        CR = f.coordinate_ring()
        dom = f.domain()
        PS = f.codomain().ambient_space()
        N = PS.dimension_relative() + 1
        R = f.base_ring()
        F_1 = f.nth_iterate_map(n+m)
        F_2 = f.nth_iterate_map(m)
        L = [F_1[i]*F_2[j] - F_1[j]*F_2[i] for i in range(N)
            for j in range(i+1, N)]
        X = PS.subscheme(L + list(dom.defining_polynomials()))
        hyperplane_at_infinity = PS.subscheme(CR.gens()[-1])
        if R.is_field():
            F = R
        else:
            F = FractionField(R)
        if X.intersection(hyperplane_at_infinity).change_ring(F).dimension() >= 0:
            hyperplane_found = False
            attempted_combinations = {}
            # in order to find a hyperplane to move to infinity,
            # we need to enumerate an infinite number of hyperplanes
            if R.is_finite():
                # when R is finite, we try all hyperplanes
                for tup in product(R, repeat=N):
                    if list(tup) != [0]*N:
                        if PS(tup) not in attempted_combinations:
                            hyperplane = PS.subscheme(sum([tup[i]*PS.gens()[i] for i in range(N)]))
                            if X.intersection(hyperplane).change_ring(F).dimension() < 0:
                                hyperplane_found = True
                                break
                if not hyperplane_found:
                    raise ValueError('no possible conjugation over %s makes all preperiodic points affine' % R)
            else:
                # if the characteristic is 0, R contains Z
                if R.characteristic() == 0:
                    for height_bound in count(1):
                        terms = ZZ.range(height_bound)
                        for tup in product(terms, repeat=N):
                            if list(tup) != [0]*N:
                                if PS(tup) not in attempted_combinations:
                                    hyperplane = PS.subscheme(sum([tup[i]*PS.gens()[i] for i in range(N)]))
                                    if X.intersection(hyperplane).change_ring(F).dimension() < 0:
                                        hyperplane_found = True
                                        break
                        if hyperplane_found:
                            break
                else:
                    if isinstance(R, (PolynomialRing_generic,
                                      MPolynomialRing_base,
                                      FractionField_generic)):
                        # for polynomial rings, we can get an infinite family of hyperplanes
                        # by increasing the degree
                        var = R.gen()
                        for degree in count(0):
                            ZZ_terms = ZZ.range(R.characteristic())
                            terms = ZZ_terms[:]
                            for i in ZZ_terms:
                                terms.append(i*var**degree)
                            for tup in product(terms, repeat=N):
                                if list(tup) != [0]*N:
                                    if PS(tup) not in attempted_combinations:
                                        hyperplane = PS.subscheme(sum([tup[i]*PS.gens()[i] for i in range(N)]))
                                        if X.intersection(hyperplane).change_ring(F).dimension() < 0:
                                            hyperplane_found = True
                                            break
                            if hyperplane_found:
                                break
                    else:
                        raise NotImplementedError('cannot find affine periodic model over %s' % (R))
            source = PS.subscheme(CR.gens()[-1])
            mat = PS.hyperplane_transformation_matrix(source, hyperplane)
            if R.is_field():
                g = f.conjugate(mat)
            else:
                g = f.conjugate(mat, adjugate=True)
            if return_conjugation:
                return (g, mat)
            return g
        if return_conjugation:
            return (f, matrix.identity(N))
        return f

    def automorphism_group(self, **kwds):
        r"""
        Calculate the subgroup of `PGL2` that is the automorphism group
        of this dynamical system.

        The automorphism group is the set of `PGL(2)` elements that fixes
        this map under conjugation.

        INPUT:

        The following keywords are used in most cases:

        - ``num_cpus`` -- (default: 2) the number of threads to use. Setting to a
          larger number can greatly speed up this function

        The following keywords are used only when the dimension of the domain is 1 and
        the base ring is the rationals, but ignored in all other cases:

        - ``starting_prime`` -- (default: 5) the first prime to use for CRT

        - ``algorithm`` -- (optional) can be one of the following:

          * ``'CRT'`` -- Chinese Remainder Theorem
          * ``'fixed_points'`` -- fixed points algorithm

        - ``return_functions`` -- boolean (default: ``False``); ``True``
          returns elements as linear fractional transformations and
          ``False`` returns elements as `PGL2` matrices

        - ``iso_type`` -- boolean (default: ``False``); ``True`` returns the
          isomorphism type of the automorphism group

        OUTPUT: list of elements in the automorphism group

        AUTHORS:

        - Original algorithm written by Xander Faber, Michelle Manes,
          Bianca Viray

        - Modified by Joao Alberto de Faria, Ben Hutz, Bianca Thompson

        REFERENCES:

        - [FMV2014]_

        EXAMPLES::

            sage: R.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSystem_projective([x^2 - y^2, x*y])
            sage: f.automorphism_group(return_functions=True)                           # needs sage.libs.pari
            [x, -x]

        ::

            sage: R.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSystem_projective([x^2 + 5*x*y + 5*y^2, 5*x^2 + 5*x*y + y^2])
            sage: f.automorphism_group()                                                # needs sage.libs.pari
            [
            [1 0]  [0 2]
            [0 1], [2 0]
            ]

        ::

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: f = DynamicalSystem([x^3, y^3, z^3])
            sage: len(f.automorphism_group())                                           # needs sage.rings.function_field
            24

        ::

            sage: R.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSystem_projective([x^2 - 2*x*y - 2*y^2, -2*x^2 - 2*x*y + y^2])
            sage: f.automorphism_group(return_functions=True)                           # needs sage.libs.pari
            [x, 1/x, -x - 1, -x/(x + 1), (-x - 1)/x, -1/(x + 1)]

        ::

            sage: R.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSystem_projective([3*x^2*y - y^3, x^3 - 3*x*y^2])
            sage: lst, label = f.automorphism_group(algorithm='CRT',                    # needs sage.libs.pari
            ....:                                   return_functions=True,
            ....:                                   iso_type=True)
            sage: sorted(lst), label                                                    # needs sage.libs.pari
            ([-1/x, 1/x, (-x - 1)/(x - 1), (-x + 1)/(x + 1), (x - 1)/(x + 1),
              (x + 1)/(x - 1), -x, x],
             'Dihedral of order 8')

        ::

            sage: A.<z> = AffineSpace(QQ, 1)
            sage: f = DynamicalSystem_affine([1/z^3])
            sage: F = f.homogenize(1)
            sage: F.automorphism_group()                                                # needs sage.libs.pari
            [
            [1 0]  [0 2]  [-1  0]  [ 0 -2]
            [0 1], [2 0], [ 0  1], [ 2  0]
            ]

        ::

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: f = DynamicalSystem_projective([x**2 + x*z, y**2, z**2])
            sage: f.automorphism_group()                                                # needs sage.rings.function_field
            [
            [1 0 0]
            [0 1 0]
            [0 0 1]
            ]

        ::

            sage: # needs sage.rings.number_field
            sage: K.<w> = CyclotomicField(3)
            sage: P.<x,y> = ProjectiveSpace(K, 1)
            sage: D6 = DynamicalSystem_projective([y^2, x^2])
            sage: sorted(D6.automorphism_group())
            [
            [-w - 1      0]  [     0 -w - 1]  [w 0]  [0 w]  [0 1]  [1 0]
            [     0      1], [     1      0], [0 1], [1 0], [1 0], [0 1]
            ]
        """
        alg = kwds.get('algorithm', None)
        p = kwds.get('starting_prime', 5)
        return_functions = kwds.get('return_functions', False)
        iso_type = kwds.get('iso_type', False)
        num_cpus = kwds.get('num_cpus', 2)
        if self.domain().dimension_relative() != 1:
            return self.conjugating_set(self, num_cpus)
        if self.base_ring() != QQ and self.base_ring() != ZZ:
            return self.conjugating_set(self, num_cpus)
        self.normalize_coordinates()
        if (self.degree() == 1) or (self.degree() == 0):
            raise NotImplementedError("rational function of degree 1 not implemented")
        f = self.dehomogenize(1)
        R = PolynomialRing(f.base_ring(),'x')
        if isinstance(f[0], FractionFieldElement):
            F = (f[0].numerator().univariate_polynomial(R))/f[0].denominator().univariate_polynomial(R)
        else:
            F = f[0].univariate_polynomial(R)
        if alg is None:
            if self.degree() <= 12:
                return automorphism_group_QQ_fixedpoints(F, return_functions, iso_type)
            return automorphism_group_QQ_CRT(F, p, return_functions, iso_type)
        elif alg == 'CRT':
            return automorphism_group_QQ_CRT(F, p, return_functions, iso_type)
        return automorphism_group_QQ_fixedpoints(F, return_functions, iso_type)

    def critical_subscheme(self):
        r"""
        Return the critical subscheme of this dynamical system.

        OUTPUT: projective subscheme

        EXAMPLES::

            sage: set_verbose(None)
            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: f = DynamicalSystem_projective([x^3 - 2*x*y^2 + 2*y^3, y^3])
            sage: f.critical_subscheme()                                                # needs sage.rings.function_field
            Closed subscheme of Projective Space of dimension 1 over Rational Field
            defined by:
              9*x^2*y^2 - 6*y^4

        ::

            sage: set_verbose(None)
            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: f = DynamicalSystem_projective([2*x^2 - y^2, x*y])
            sage: f.critical_subscheme()                                                # needs sage.rings.function_field
            Closed subscheme of Projective Space of dimension 1 over Rational Field
            defined by:
              4*x^2 + 2*y^2

        ::

            sage: P.<x,y,z> = ProjectiveSpace(QQ,2)
            sage: f = DynamicalSystem_projective([2*x^2 - y^2, x*y, z^2])
            sage: f.critical_subscheme()                                                # needs sage.rings.function_field
            Closed subscheme of Projective Space of dimension 2 over Rational Field
            defined by:
              8*x^2*z + 4*y^2*z

        ::

            sage: # needs sage.rings.finite_rings
            sage: P.<x,y,z,w> = ProjectiveSpace(GF(81), 3)
            sage: g = DynamicalSystem_projective([x^3 + y^3, y^3 + z^3, z^3 + x^3, w^3])
            sage: g.critical_subscheme()
            Closed subscheme of Projective Space of dimension 3 over Finite Field in
            z4 of size 3^4 defined by:
              0

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: f = DynamicalSystem_projective([x^2, x*y])
            sage: f.critical_subscheme()                                                # needs sage.rings.function_field
            Traceback (most recent call last):
            ...
            TypeError: the function is not a morphism
        """
        PS = self.domain()
        if not isinstance(PS, ProjectiveSpace_ring):
            raise NotImplementedError("not implemented for subschemes")
        if not self.is_morphism():
            raise TypeError("the function is not a morphism")
        wr = self.wronskian_ideal()
        crit_subscheme = self.codomain().subscheme(wr)
        return crit_subscheme

    def critical_points(self, R=None):
        r"""
        Return the critical points of this dynamical system defined over
        the ring `R` or the base ring of this map.

        Must be dimension 1.

        INPUT:

        - ``R`` -- (optional) a ring

        OUTPUT: list of projective space points defined over `R`

        EXAMPLES::

            sage: set_verbose(None)
            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: f = DynamicalSystem_projective([x^3 - 2*x*y^2 + 2*y^3, y^3])
            sage: f.critical_points()                                                   # needs sage.rings.function_field
            [(1 : 0)]
            sage: K.<w> = QuadraticField(6)                                             # needs sage.rings.number_field
            sage: f.critical_points(K)                                                  # needs sage.rings.number_field
            [(-1/3*w : 1), (1/3*w : 1), (1 : 0)]

        ::

            sage: set_verbose(None)
            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: f = DynamicalSystem_projective([2*x^2 - y^2, x*y])
            sage: f.critical_points(QQbar)                                              # needs sage.rings.number_field
            [(-0.7071067811865475?*I : 1), (0.7071067811865475?*I : 1)]
        """
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

    def ramification_type(self, R=None, stable=True):
        r"""
        Return the ramification type of endomorphisms of `\mathbb{P}^1`.

        Only branch points defined over the ring ``R`` contribute to
        the ramification type if specified, otherwise ``R`` is the
        ring of definition for ``self``.

        Note that branch points defined over ``R`` may not be
        geometric points if stable not set to ``True``.

        If ``R`` is specified, ``stable`` is ignored.

        If ``stable``, then this will return the ramification type
        over an extension which splits the Galois orbits of critical
        points.

        INPUT:

        - ``R`` -- ring or morphism (optional)
        - ``split`` -- boolean (optional)

        OUTPUT:

        list of lists, each term being the list of ramification indices
        in the pre-images of one critical value

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: F = DynamicalSystem_projective([x^4, y^4])
            sage: F.ramification_type()                                                 # needs sage.rings.number_field
            [[4], [4]]

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: F = DynamicalSystem_projective([x^3, 4*y^3 - 3*x^2*y])
            sage: F.ramification_type()                                                 # needs sage.rings.number_field
            [[2], [2], [3]]

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: F = DynamicalSystem_projective([(x + y)^4, 16*x*y*(x-y)^2])
            sage: F.ramification_type()                                                 # needs sage.rings.number_field
            [[2], [2, 2], [4]]

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: F = DynamicalSystem_projective([(x + y)*(x - y)^3, y*(2*x+y)^3])
            sage: F.ramification_type()                                                 # needs sage.rings.number_field
            [[3], [3], [3]]

            sage: F = DynamicalSystem_projective([x^3 - 2*x*y^2 + 2*y^3, y^3])
            sage: F.ramification_type()                                                 # needs sage.rings.number_field
            [[2], [2], [3]]
            sage: F.ramification_type(R=F.base_ring())                                  # needs sage.rings.function_field
            [[2], [3]]
        """
        # Change base ring if specified.
        if R is None:
            if stable:
                L,phi = self.field_of_definition_critical(return_embedding=True)
                F = self.change_ring(phi)
            else:
                F = self
        else:
            F = self.change_ring(R)

        C = F.critical_subscheme()
        ram_type = {}
        fc = C.defining_ideal().gens()[0]
        for f, m in fc.factor():
            c = F(F.domain().subscheme(f))  # critical value
            if c in ram_type:
                ram_type[c].append(m + 1)
            else:
                ram_type[c] = [m + 1]
        return sorted(ram_type.values())

    def is_postcritically_finite(self, err=0.01, use_algebraic_closure=True):
        r"""
        Determine if this dynamical system is post-critically finite.

        Only for endomorphisms of `\mathbb{P}^1`. It checks if each critical
        point is preperiodic. The optional parameter ``err`` is passed into
        :meth:`is_preperiodic` as part of the preperiodic check.

        The computations can be done either over the algebraic closure of the
        base field or over the minimal extension of the base field that
        contains the critical points.

        INPUT:

        - ``err`` -- (default: 0.01) positive real number

        - ``use_algebraic_closure`` -- boolean (default: ``True``); if ``True``, uses the
          algebraic closure. If ``False``, uses the smallest extension of the base field
          containing all the critical points.

        OUTPUT: boolean

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: f = DynamicalSystem_projective([x^2 - y^2, y^2])
            sage: f.is_postcritically_finite()                                          # needs sage.rings.number_field
            True

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: f = DynamicalSystem_projective([x^3 - y^3, y^3])
            sage: f.is_postcritically_finite()                                          # needs sage.rings.number_field
            False

        ::

            sage: # needs sage.rings.number_field
            sage: R.<z> = QQ[]
            sage: K.<v> = NumberField(z^8 + 3*z^6 + 3*z^4 + z^2 + 1)
            sage: PS.<x,y> = ProjectiveSpace(K,1)
            sage: f = DynamicalSystem_projective([x^3 + v*y^3, y^3])
            sage: f.is_postcritically_finite()  # long time
            True

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: f = DynamicalSystem_projective([6*x^2 + 16*x*y + 16*y^2,
            ....:                                 -3*x^2 - 4*x*y - 4*y^2])
            sage: f.is_postcritically_finite()                                          # needs sage.rings.number_field
            True

        ::

            sage: # needs sage.libs.gap sage.rings.number_field
            sage: K = UniversalCyclotomicField()
            sage: P.<x,y> = ProjectiveSpace(K,1)
            sage: F = DynamicalSystem_projective([x^2 - y^2, y^2], domain=P)
            sage: F.is_postcritically_finite()
            True

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: f = DynamicalSystem_projective([8*x^4 - 8*x^2*y^2 + y^4, y^4])
            sage: f.is_postcritically_finite(use_algebraic_closure=False)  # long time
            True

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: f = DynamicalSystem_projective([x^4 - x^2*y^2 + y^4, y^4])
            sage: f.is_postcritically_finite(use_algebraic_closure=False)               # needs sage.rings.number_field
            False

        ::

            sage: # needs sage.rings.number_field
            sage: P.<x,y> = ProjectiveSpace(QQbar,1)
            sage: f = DynamicalSystem_projective([x^4 - x^2*y^2, y^4])
            sage: f.is_postcritically_finite()
            False
        """
        #iteration of subschemes not yet implemented
        if self.domain().dimension_relative() > 1:
            raise NotImplementedError("only implemented in dimension 1")

        K = FractionField(self.codomain().base_ring())
        if use_algebraic_closure:
            Kbar = K.algebraic_closure()
            if Kbar.has_coerce_map_from(K):
                F = self.change_ring(Kbar)
            else:
                embeds = K.embeddings(Kbar)
                if embeds:
                    F = self.change_ring(embeds[0])
                else:
                    raise ValueError("no embeddings of base field to algebraic closure")
        else:
            embedding = self.field_of_definition_critical(return_embedding=True)[1]
            F = self.change_ring(embedding)

        crit_points = F.critical_points()
        pcf = True
        i = 0
        while pcf and i < len(crit_points):
            if not crit_points[i].is_preperiodic(F, err):
                pcf = False
            i += 1
        return pcf

    def is_dynamical_belyi_map(self):
        r"""
        Return if this dynamical system is a dynamical Belyi map.

        We define a dynamical Belyi map to be a map conjugate to a
        dynamical system `f: \mathbb{P}^1  \to \mathbb{P}^1`
        where the branch points are contained in `\{0, 1, \infty \}`
        and the postcritical set is contained in `\{0, 1, \infty \}`.

        Output: Boolean

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSystem_projective([-2*x^3 - 9*x^2*y - 12*x*y^2 - 6*y^3, y^3])
            sage: f.is_dynamical_belyi_map()                                            # needs sage.rings.number_field
            True

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSystem_projective([5*x^7 - 7*x^6*y, -7*x*y^6 + 5*y^7])
            sage: f.is_dynamical_belyi_map()                                            # needs sage.rings.number_field
            True

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSystem_projective([x^2 + y^2, y^2])
            sage: f.is_dynamical_belyi_map()                                            # needs sage.rings.number_field
            False

        ::

            sage: # needs sage.rings.number_field
            sage: F = QuadraticField(-7)
            sage: P.<x,y> = ProjectiveSpace(F, 1)
            sage: f = DynamicalSystem_projective([5*x^7 - 7*x^6*y, -7*x*y^6 + 5*y^7])
            sage: f.is_dynamical_belyi_map()
            True

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSystem_projective([2*x^3 + 3*x^2*y - 3*x*y^2 + 2*y^3,
            ....:                                 x^3 + y^3])
            sage: f.is_dynamical_belyi_map()                                            # needs sage.rings.number_field
            False

        ::

            sage: # needs sage.rings.number_field
            sage: R.<t> = PolynomialRing(QQ)
            sage: N.<c> = NumberField(t^3 - 2)
            sage: P.<x,y> = ProjectiveSpace(N, 1)
            sage: f=DynamicalSystem_projective([x^2 + c*y^2, x*y])
            sage: f.is_dynamical_belyi_map()
            False

        ::

            sage: P.<x,y> = ProjectiveSpace(GF(7), 1)
            sage: f = DynamicalSystem_projective([x^3 + 6*y^3, y^3])
            sage: f.is_dynamical_belyi_map()                                            # needs sage.libs.pari
            False
        """
        P = self.codomain()
        if not isinstance(P, ProjectiveSpace_ring):
            raise NotImplementedError('only implemented for dynamical systems on projective space')
        if P.dimension_relative() != 1:
            raise NotImplementedError('only implemented for maps on projective space of dimension 1')
        embed = self.field_of_definition_critical(return_embedding=True)[1]
        f = self.change_ring(embed)
        crit_list = f.critical_points()
        crit_orbit = []
        for i in crit_list:
            crit_orbit += f.orbit(i, 4)
        if len(set(crit_orbit)) > 3:
            return False
        return True

    def critical_point_portrait(self, check=True, use_algebraic_closure=True):
        r"""
        If this dynamical system  is post-critically finite, return its
        critical point portrait.

        This is the directed graph of iterates starting with the critical
        points. Must be dimension 1. If ``check`` is ``True``, then the
        map is first checked to see if it is postcritically finite.

        The computations can be done either over the algebraic closure of the
        base field or over the minimal extension of the base field that
        contains the critical points.

        INPUT:

        - ``check`` -- boolean (default: ``True``)

        - ``use_algebraic_closure`` -- boolean (default: ``True``); if ``True``, uses the
          algebraic closure. If ``False``, uses the smallest extension of the base field
          containing all the critical points.

        OUTPUT: a digraph

        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: R.<z> = QQ[]
            sage: K.<v> = NumberField(z^6 + 2*z^5 + 2*z^4 + 2*z^3 + z^2 + 1)
            sage: PS.<x,y> = ProjectiveSpace(K,1)
            sage: f = DynamicalSystem_projective([x^2 + v*y^2, y^2])
            sage: f.critical_point_portrait(check=False)        # long time             # needs sage.graphs
            Looped digraph on 6 vertices

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: f = DynamicalSystem_projective([x^5 + 5/4*x*y^4, y^5])
            sage: f.critical_point_portrait(check=False)                                # needs sage.graphs sage.rings.number_field
            Looped digraph on 5 vertices

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: f = DynamicalSystem_projective([x^2 + 2*y^2, y^2])
            sage: f.critical_point_portrait()                                           # needs sage.rings.number_field
            Traceback (most recent call last):
            ...
            TypeError: map must be post-critically finite

        ::

            sage: # needs sage.rings.number_field
            sage: R.<t> = QQ[]
            sage: K.<v> = NumberField(t^3 + 2*t^2 + t + 1)
            sage: phi = K.embeddings(QQbar)[0]
            sage: P.<x, y> = ProjectiveSpace(K, 1)
            sage: f = DynamicalSystem_projective([x^2 + v*y^2, y^2])
            sage: f.change_ring(phi).critical_point_portrait()                          # needs sage.graphs
            Looped digraph on 4 vertices

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: f = DynamicalSystem_projective([8*x^4 - 8*x^2*y^2 + y^4, y^4])
            sage: f.critical_point_portrait(use_algebraic_closure=False)  # long time
            Looped digraph on 6 vertices

        ::

            sage: # needs sage.rings.number_field
            sage: P.<x,y> = ProjectiveSpace(QQbar,1)
            sage: f = DynamicalSystem_projective([8*x^4 - 8*x^2*y^2 + y^4, y^4])
            sage: f.critical_point_portrait()   # long time                             # needs sage.graphs
            Looped digraph on 6 vertices

        ::

            sage: P.<x,y> = ProjectiveSpace(GF(3),1)
            sage: f = DynamicalSystem_projective([x^2 + x*y - y^2, x*y])
            sage: f.critical_point_portrait(use_algebraic_closure=False)                # needs sage.libs.pari
            Looped digraph on 6 vertices
            sage: f.critical_point_portrait() #long time
            Looped digraph on 6 vertices
        """
        #input checking done in is_postcritically_finite
        if check:
            if not self.is_postcritically_finite():
                raise TypeError("map must be post-critically finite")
        K = FractionField(self.base_ring())
        if use_algebraic_closure:
            Kbar = K.algebraic_closure()
            if Kbar.has_coerce_map_from(K):
                F = self.change_ring(Kbar)
            else:
                embeds = K.embeddings(Kbar)
                if embeds:
                    F = self.change_ring(embeds[0])
                else:
                    raise ValueError("no embeddings of base field to algebraic closure")
        else:
            embedding = self.field_of_definition_critical(return_embedding=True)[1]
            F = self.change_ring(embedding)
        crit_points = F.critical_points()
        N = len(crit_points)
        for i in range(N):
            done = False
            Q = F(crit_points[i])
            while not done:
                if Q in crit_points:
                    done = True
                else:
                    crit_points.append(Q)
                Q = F(Q)
        return F._preperiodic_points_to_cyclegraph(crit_points)

    def critical_height(self, **kwds):
        r"""
        Compute the critical height of this dynamical system.

        The critical height is defined by J. Silverman as
        the sum of the canonical heights of the critical points.
        This must be dimension 1 and defined over a number field
        or number field order.

        The computations can be done either over the algebraic closure of the
        base field or over the minimal extension of the base field that
        contains the critical points.

        INPUT: keyword arguments:

        - ``badprimes`` -- (optional) a list of primes of bad reduction

        - ``N`` -- (default: 10) positive integer; number of terms of
          the series to use in the local green functions

        - ``prec`` -- (default: 100) positive integer, float point
          or `p`-adic precision

        - ``error_bound`` -- (optional) a positive real number

        - ``use_algebraic_closure`` -- boolean (default: ``True``); if
          ``True``, uses the algebraic closure. If ``False``, uses the
          smallest extension of the base field containing all the critical
          points.

        OUTPUT: real number

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: f = DynamicalSystem_projective([x^3 + 7*y^3, 11*y^3])
            sage: f.critical_height()                                                   # needs sage.rings.number_field
            1.1989273321156851418802151128

        ::

            sage: # needs sage.rings.number_field
            sage: K.<w> = QuadraticField(2)
            sage: P.<x,y> = ProjectiveSpace(K,1)
            sage: f = DynamicalSystem_projective([x^2 + w*y^2, y^2])
            sage: f.critical_height()
            0.16090842452312941163719755472

        Postcritically finite maps have critical height 0::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: f = DynamicalSystem_projective([x^3 - 3/4*x*y^2 + 3/4*y^3, y^3])
            sage: f.critical_height(error_bound=0.0001)                                 # needs sage.rings.number_field
            0.00000000000000000000000000000

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: f = DynamicalSystem_projective([x^3 + 3*x*y^2, y^3])
            sage: f.critical_height(use_algebraic_closure=False)                        # needs sage.rings.number_field
            0.000023477016733897112886491967991
            sage: f.critical_height()                                                   # needs sage.rings.number_field
            0.000023477016733897112886491967991
        """
        PS = self.codomain()
        if PS.dimension_relative() > 1:
            raise NotImplementedError("only implemented in dimension 1")

        K = FractionField(PS.base_ring())
        use_algebraic_closure = kwds.get("use_algebraic_closure", True)
        if use_algebraic_closure:
            Kbar = K.algebraic_closure()
            if Kbar.has_coerce_map_from(K):
                F = self.change_ring(Kbar)
            else:
                embeds = K.embeddings(Kbar)
                if embeds:
                    F = self.change_ring(embeds[0])
                else:
                    raise ValueError("no embeddings of base field to algebraic closure")
        else:
            embedding = self.field_of_definition_critical(return_embedding=True)[1]
            F = self.change_ring(embedding)
        crit_points = F.critical_points()
        n = len(crit_points)
        err_bound = kwds.get("error_bound", None)
        if err_bound is not None:
            kwds["error_bound"] = err_bound / n
        ch = 0
        for P in crit_points:
            ch += F.canonical_height(P, **kwds)
        return ch

    def preperiodic_points(self, m, n, **kwds):
        r"""
        Compute the preperiodic points of period ``m, n`` of this dynamical system
        defined over the ring ``R`` or the base ring of the map.

        This is done by finding the rational points on the variety
        defining the points of period ``m, n``.

        For rational maps, where there are potentially infinitely many periodic
        points of a given period, you must use the ``return_scheme`` option.
        Note that this scheme will include the indeterminacy locus.

        INPUT:

        - ``n`` -- positive integer; the period

        - ``m`` -- nonnegative integer; the preperiod

        kwds:

        - ``minimal`` -- boolean (default: ``True``); ``True`` specifies to
          find only the preperiodic points of minimal period ``m``,``n`` and
          ``False`` specifies to find all preperiodic points of period
          ``m``, ``n``

        - ``formal`` -- boolean (default: ``False``); ``True`` specifies to
          find the formal periodic points only. The formal periodic points
          are the points in the support of the dynatomic cycle.

        - ``R`` -- (default: the base ring of the dynamical system) a
          commutative ring over which to find the preperiodic points

        - ``return_scheme`` -- boolean (default: ``False``); return a
          subscheme of the ambient space that defines the ``m``,``n`` th
          preperiodic points

        OUTPUT:

        A list of preperiodic points of this map or the subscheme defining
        the preperiodic points.

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQbar, 1)                                   # needs sage.rings.number_field
            sage: f = DynamicalSystem_projective([x^2 - y^2, y^2])                      # needs sage.rings.number_field
            sage: f.preperiodic_points(0, 1)                                            # needs sage.rings.number_field
            [(-0.618033988749895? : 1), (1 : 0), (1.618033988749895? : 1)]

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSystem_projective([x^2 - 29/16*y^2, y^2])
            sage: f.preperiodic_points(1, 3)                                            # needs sage.rings.function_field
            [(-5/4 : 1), (1/4 : 1), (7/4 : 1)]

        ::

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: f = DynamicalSystem_projective([x^2 - 3/4*y^2, y^2 , z^2])
            sage: f.preperiodic_points(0, 2, formal=True)                               # needs sage.rings.function_field
            [(-1/2 : 1 : 0), (-1/2 : 1 : 1)]

        ::

            sage: # needs sage.rings.number_field
            sage: P.<x,y> = ProjectiveSpace(QQbar, 1)
            sage: f = DynamicalSystem_projective([x^2 - x*y + 2*y^2, x^2 - y^2])
            sage: f.preperiodic_points(1, 2, minimal=False)                             # needs sage.rings.function_field
            [(-3.133185666641252? : 1),
            (-1 : 1),
            (-0.3478103847799310? - 1.028852254136693?*I : 1),
            (-0.3478103847799310? + 1.028852254136693?*I : 1),
            (0.8165928333206258? - 0.6710067557437100?*I : 1),
            (0.8165928333206258? + 0.6710067557437100?*I : 1),
            (1 : 0),
            (1 : 1),
            (1.695620769559862? : 1),
            (3 : 1)]

        ::

            sage: # needs sage.rings.number_field
            sage: R.<w> = QQ[]
            sage: K.<s> = NumberField(w^6 - 3*w^5 + 5*w^4 - 5*w^3 + 5*w^2 - 3*w + 1)
            sage: P.<x,y,z> = ProjectiveSpace(K, 2)
            sage: f = DynamicalSystem_projective([x^2 + z^2, y^2 + x^2, z^2 + y^2])
            sage: sorted(f.preperiodic_points(0, 1), key=str)                           # needs sage.rings.function_field
            [(-2*s^5 + 4*s^4 - 5*s^3 + 3*s^2 - 4*s : -2*s^5 + 5*s^4 - 7*s^3 + 6*s^2 - 7*s + 3 : 1),
             (-s^5 + 3*s^4 - 4*s^3 + 4*s^2 - 4*s + 2 : -s^5 + 2*s^4 - 2*s^3 + s^2 - s : 1),
             (-s^5 + 3*s^4 - 5*s^3 + 4*s^2 - 3*s + 1 : s^5 - 2*s^4 + 3*s^3 - 3*s^2 + 4*s - 1 : 1),
             (1 : 1 : 1),
             (2*s^5 - 6*s^4 + 9*s^3 - 8*s^2 + 7*s - 4 : 2*s^5 - 5*s^4 + 7*s^3 - 5*s^2 + 6*s - 2 : 1),
             (s^5 - 2*s^4 + 2*s^3 + s : s^5 - 3*s^4 + 4*s^3 - 3*s^2 + 2*s - 1 : 1),
             (s^5 - 2*s^4 + 3*s^3 - 3*s^2 + 3*s - 1 : -s^5 + 3*s^4 - 5*s^3 + 4*s^2 - 4*s + 2 : 1)]

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSystem_projective([x^2 + 1/4*y^2, y^2])
            sage: f.preperiodic_points(1, 1, formal=True)
            [(-1/2 : 1), (1 : 0)]

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSystem_projective([x^2 - 3/4*y^2, y^2])
            sage: f.preperiodic_points(0, 2, formal=True)                               # needs sage.libs.pari
            [(-1/2 : 1)]

        ::

            sage: # needs sage.rings.number_field
            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: K.<v> = QuadraticField(5)
            sage: phi = QQ.embeddings(K)[0]
            sage: f = DynamicalSystem_projective([x^2 - y^2, y^2])
            sage: f.preperiodic_points(1, 1, R=phi)
            [(-1/2*v - 1/2 : 1), (1/2*v - 1/2 : 1)]

        ::

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: X = P.subscheme(2*x - y)
            sage: f = DynamicalSystem_projective([x^2 - y^2, 2*(x^2 - y^2), y^2 - z^2],
            ....:                                domain=X)
            sage: f.preperiodic_points(1, 1)                                            # needs sage.rings.function_field
            [(-1/4 : -1/2 : 1), (1 : 2 : 1)]

        ::

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: f = DynamicalSystem_projective([x^2 - 3/4*y^2, z^2, y^2])
            sage: f.preperiodic_points(1, 1)                                            # needs sage.rings.function_field
            [(-3/2 : -1 : 1), (-3/2 : 1 : 1), (-1/2 : -1 : 1), (1/2 : -1 : 1),
             (1/2 : 1 : 1), (3/2 : -1 : 1)]

        ::

            sage: P.<x,y,z> = ProjectiveSpace(GF(5), 2)
            sage: f = DynamicalSystem_projective([x^2, y^2, z^2])
            sage: sorted(f.preperiodic_points(2, 1))                                    # needs sage.rings.function_field
            [(0 : 2 : 1), (0 : 3 : 1), (1 : 2 : 1), (1 : 3 : 1), (2 : 0 : 1), (2 : 1 : 0),
             (2 : 1 : 1), (2 : 2 : 1), (2 : 3 : 1), (2 : 4 : 1), (3 : 0 : 1), (3 : 1 : 0),
             (3 : 1 : 1), (3 : 2 : 1), (3 : 3 : 1), (3 : 4 : 1), (4 : 2 : 1), (4 : 3 : 1)]

        ::

            sage: P.<x,y,z> = ProjectiveSpace(GF(5), 2)
            sage: f = DynamicalSystem_projective([x^2, x*y, z^2])
            sage: f.preperiodic_points(2, 1, return_scheme=True, minimal=False)
            Closed subscheme of
             Projective Space of dimension 2 over Finite Field of size 5 defined by:
              0,
              x^8*z^4 - x^4*z^8,
              x^7*y*z^4 - x^3*y*z^8

        When the ring over which to find the preperiodic points is a number
        field, the ordering of the preperiodic points might depend on the
        architecture (32 or 64 bits)::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: R.<z> = QQ[]
            sage: K.<v> = NumberField(z^4 - z^2 - 1)                                    # needs sage.rings.number_field
            sage: f = DynamicalSystem_projective([x^2 - y^2, y^2])
            sage: sorted(f.preperiodic_points(2, 1, R=K), key=str)                      # needs sage.rings.number_field
            [(-v : 1), (v : 1)]

        ::

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: f = DynamicalSystem_projective([x^2 - 3/4*y^2, y^2, z^2])
            sage: f.preperiodic_points(0, 2, formal=True)                               # needs sage.rings.function_field
            [(-1/2 : 1 : 0), (-1/2 : 1 : 1)]

        ::

            sage: S.<c> = QQ[]
            sage: R.<x,y> = PolynomialRing(S, 2)
            sage: P = ProjectiveSpace(R)
            sage: f = DynamicalSystem_projective([x^2 + c*y^2, y^2])
            sage: f.preperiodic_points(1, 2, return_scheme=True)                        # needs sage.rings.function_field
            Closed subscheme of Projective Space of dimension 1 over Univariate
             Polynomial Ring in c over Rational Field defined by:
              x^2 - x*y + (c + 1)*y^2

        TESTS::

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: f = DynamicalSystem_projective([x^2, x*y, z^2])
            sage: f.preperiodic_points(2, 1, minimal=False)                             # needs sage.rings.function_field
            Traceback (most recent call last):
            ...
            TypeError: use return_scheme=True

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSystem_projective([x^2 - 29/16*y^2, y^2])
            sage: f.preperiodic_points(1.2, 3)
            Traceback (most recent call last):
            ...
            TypeError: Attempt to coerce non-integral RealNumber to Integer
            sage: f.preperiodic_points(1, 3.1)
            Traceback (most recent call last):
            ...
            TypeError: Attempt to coerce non-integral RealNumber to Integer

        ::

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: f = DynamicalSystem([x^2 - z^2, y^2 - 21/16*z^2, z^2])
            sage: len(f.preperiodic_points(1, 2, minimal=True, formal=False)) == 16     # needs sage.rings.function_field
            True

        ::

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: f = DynamicalSystem_projective([x^2 - y^2, 2*(x^2 - y^2), y^2 - z^2])
            sage: f.preperiodic_points(2, 2)                                            # needs sage.rings.function_field
            Traceback (most recent call last):
            ...
            ValueError: dynamical system is not a morphism,
            cannot calculate minimal or formal preperiodic points
        """
        n = ZZ(n)
        m = ZZ(m)
        if n <= 0:
            raise ValueError("a positive integer period must be specified")
        if m < 0:
            raise ValueError("a nonnegative preperiod must be specified")
        R = kwds.pop('R', None)
        if R is None:
            f_sub = self
            R = self.base_ring()
        else:
            f_sub = self.change_ring(R)
            R = f_sub.base_ring() #in the case when R is an embedding
        if isinstance(R, FractionField_1poly_field) or R in FunctionFields():
            raise NotImplementedError('Periodic points not implemented for function fields; '
                'clear denominators and use the polynomial ring instead')
        CR = f_sub.coordinate_ring()
        dom = f_sub.domain()
        PS = f_sub.codomain().ambient_space()
        if dom != PS:
            f = DynamicalSystem(f_sub.defining_polynomials())
        else:
            f = f_sub
        N = PS.dimension_relative() + 1
        formal = kwds.pop('formal', False)
        minimal = kwds.pop('minimal', True)
        return_scheme = kwds.pop('return_scheme', False)
        if formal and N == 2 and dom == PS:
            X = PS.subscheme([f.dynatomic_polynomial([m,n])])
        else:
            F_1 = f.nth_iterate_map(n+m)
            F_2 = f.nth_iterate_map(m)
            L = [F_1[i]*F_2[j] - F_1[j]*F_2[i] for i in range(N)
                    for j in range(i+1, N)]
            X = PS.subscheme(L + list(dom.defining_polynomials()))
            if (minimal or formal) and (n != 1 or m != 0):
                if not f_sub.is_morphism():
                    raise ValueError('dynamical system is not a morphism, cannot calculate minimal or formal preperiodic points')
                if formal:
                    d = f.degree()
                    # we need a model with no preperiodic points at infinity
                    new_f, mat = f.affine_preperiodic_model(m, n, return_conjugation=True)
                    new_f.normalize_coordinates()
                    # we now deform by a parameter t
                    T = R['t']
                    t = T.gens()[0]
                    Pt = ProjectiveSpace(N-1, R=T, names=[str(i) for i in CR.gens()])
                    deformed_polys = [poly + t*Pt.gens()[-1]**d for poly in new_f.defining_polynomials()[:-1]]
                    deformed_polys += [new_f.defining_polynomials()[-1]]
                    f_deformed = DynamicalSystem(deformed_polys)

                    # after deforming by the parameter, the preperiodic points with multiplicity
                    # will separate into different points. we can now calculate the minimal preperiodic
                    # points with the parameter, and then specialize to get the formal preperiodic points
                    ideal = f_deformed.preperiodic_points(m, n, return_scheme=True).defining_ideal()
                    L = [poly.specialization({t:0}) for poly in ideal.gens()]
                    X = PS.subscheme(L)
                    subs_list = mat.inverse()*vector(CR.gens())
                    subs = {}
                    for i in range(len(subs_list)):
                        subs[PS.gens()[i]] = subs_list[i]
                    if R.is_field():
                        X = PS.subscheme([poly.subs(subs) for poly in L])
                    else:
                        K = [poly.subs(subs) for poly in L]
                        K = [poly*poly.denominator() for poly in K]
                        X = PS.subscheme(K)
                elif minimal:
                    Sn = []
                    for k in ZZ(n).divisors():
                        if ZZ(n/k).is_prime():
                            Sn.append(k)
                    if isinstance(R, (PolynomialRing_generic,
                                      MPolynomialRing_base)):
                        phi = FlatteningMorphism(CR)
                        flatCR = phi.codomain()
                        Ik = flatCR.ideal(1)
                        for k in Sn:
                            Ik *= f.preperiodic_points(m, k, return_scheme=True, minimal=False).defining_ideal()
                        if m != 0:
                            Ik *= f.preperiodic_points(m-1, n, return_scheme=True, minimal=False).defining_ideal()
                        psi = UnflatteningMorphism(flatCR, CR)
                        In = flatCR.ideal([phi(i) for i in X.defining_polynomials()])
                        X = PS.subscheme([psi(i) for i in In.saturation(Ik)[0].gens()])
                    else:
                        Ik = CR.ideal(1)
                        for k in Sn:
                            Ik *= f.preperiodic_points(m, k, return_scheme=True, minimal=False).defining_ideal()
                        if m != 0:
                            Ik *= f.preperiodic_points(m-1, n, return_scheme=True, minimal=False).defining_ideal()
                        In = X.defining_ideal()
                        X = PS.subscheme(In.saturation(Ik)[0])
        if dom != PS:
            X = PS.subscheme(list(X.defining_polynomials()) + list(dom.defining_polynomials()))
        if return_scheme:  # this includes the indeterminacy locus points!
            return X
        if X.dimension() <= 0:
            if R in NumberFields() or isinstance(R, sage.rings.abc.AlgebraicField) or R in FiniteFields():
                Z = f.base_indeterminacy_locus()
                points = [dom(Q) for Q in X.rational_points()]
                good_points = []
                for Q in points:
                    try:
                        Z(list(Q))
                    except TypeError:
                        good_points.append(Q)
                good_points.sort()
                return good_points
            else:
                raise NotImplementedError("ring must a number field or finite field")
        else: #a higher dimensional scheme
            raise TypeError("use return_scheme=True")

    def periodic_points(self, n, minimal=True, formal=False, R=None, algorithm='variety',
                        return_scheme=False):
        r"""
        Compute the periodic points of period ``n`` of this dynamical system
        defined over the ring ``R`` or the base ring of the map.

        This can be done either by finding the rational points on the variety
        defining the points of period ``n``, or, for finite fields,
        finding the cycle of appropriate length in the cyclegraph. For small
        cardinality fields, the cyclegraph algorithm is effective for any
        map and length cycle, but is slow when the cyclegraph is large.
        The variety algorithm is good for small period, degree, and dimension,
        but is slow as the defining equations of the variety get more
        complicated.

        For rational maps, where there are potentially infinitely many periodic
        points of a given period, you must use the ``return_scheme`` option.
        Note that this scheme will include the indeterminacy locus.

        INPUT:

        - ``n`` -- positive integer

        - ``minimal`` -- boolean (default: ``True``); ``True`` specifies to
          find only the periodic points of minimal period ``n`` and ``False``
          specifies to find all periodic points of period ``n``

        - ``formal`` -- boolean (default: ``False``); ``True`` specifies to
          find the formal periodic points only. The formal periodic points
          are the points in the support of the dynatomic cycle.

        - ``R`` -- (optional) a commutative ring. Defaults to the base ring of
          this map

        - ``algorithm`` -- (default: ``'variety'``) must be one of
          the following:

          * ``'variety'`` -- find the rational points on the appropriate variety
          * ``'cyclegraph'`` -- find the cycles from the cycle graph

        - ``return_scheme`` -- return a subscheme of the ambient space
          that defines the ``n`` th periodic points

        OUTPUT:

        A list of periodic points of this map or the subscheme defining
        the periodic points.

        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: set_verbose(None)
            sage: P.<x,y> = ProjectiveSpace(QQbar, 1)
            sage: f = DynamicalSystem_projective([x^2 - x*y + y^2, x^2 - y^2 + x*y])
            sage: f.periodic_points(1)
            [(-0.50000000000000000? - 0.866025403784439?*I : 1),
             (-0.50000000000000000? + 0.866025403784439?*I : 1),
             (1 : 1)]

        ::

            sage: # needs sage.rings.number_field
            sage: P.<x,y,z> = ProjectiveSpace(QuadraticField(5,'t'), 2)
            sage: f = DynamicalSystem_projective([x^2 - 21/16*z^2, y^2 - z^2, z^2])
            sage: f.periodic_points(2)
            [(-5/4 : -1 : 1), (-5/4 : -1/2*t + 1/2 : 1), (-5/4 : 0 : 1),
             (-5/4 : 1/2*t + 1/2 : 1), (-3/4 : -1 : 1), (-3/4 : 0 : 1),
             (1/4 : -1 : 1), (1/4 : -1/2*t + 1/2 : 1), (1/4 : 0 : 1),
             (1/4 : 1/2*t + 1/2 : 1), (7/4 : -1 : 1), (7/4 : 0 : 1)]

        ::

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: f = DynamicalSystem_projective([x^2 - 3/4*y^2, y^2 , z^2])
            sage: f.periodic_points(2, formal=True)                                     # needs sage.rings.function_field
            [(-1/2 : 1 : 0), (-1/2 : 1 : 1)]

        ::

            sage: # needs sage.rings.number_field
            sage: w = QQ['w'].0
            sage: K = NumberField(w^6 - 3*w^5 + 5*w^4 - 5*w^3 + 5*w^2 - 3*w + 1,'s')
            sage: P.<x,y,z> = ProjectiveSpace(K, 2)
            sage: f = DynamicalSystem_projective([x^2 + z^2, y^2 + x^2, z^2 + y^2])
            sage: sorted(f.periodic_points(1), key=str)                                 # needs sage.rings.function_field
            [(-2*s^5 + 4*s^4 - 5*s^3 + 3*s^2 - 4*s : -2*s^5 + 5*s^4 - 7*s^3 + 6*s^2 - 7*s + 3 : 1),
             (-s^5 + 3*s^4 - 4*s^3 + 4*s^2 - 4*s + 2 : -s^5 + 2*s^4 - 2*s^3 + s^2 - s : 1),
             (-s^5 + 3*s^4 - 5*s^3 + 4*s^2 - 3*s + 1 : s^5 - 2*s^4 + 3*s^3 - 3*s^2 + 4*s - 1 : 1),
             (1 : 1 : 1),
             (2*s^5 - 6*s^4 + 9*s^3 - 8*s^2 + 7*s - 4 : 2*s^5 - 5*s^4 + 7*s^3 - 5*s^2 + 6*s - 2 : 1),
             (s^5 - 2*s^4 + 2*s^3 + s : s^5 - 3*s^4 + 4*s^3 - 3*s^2 + 2*s - 1 : 1),
             (s^5 - 2*s^4 + 3*s^3 - 3*s^2 + 3*s - 1 : -s^5 + 3*s^4 - 5*s^3 + 4*s^2 - 4*s + 2 : 1)]

        ::

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: f = DynamicalSystem_projective([x^2 - 21/16*z^2, y^2 - 2*z^2, z^2])
            sage: f.periodic_points(2, False)                                           # needs sage.rings.function_field
            [(-5/4 : -1 : 1), (-5/4 : 2 : 1), (-3/4 : -1 : 1),
             (-3/4 : 2 : 1), (0 : 1 : 0), (1/4 : -1 : 1), (1/4 : 2 : 1),
             (1 : 0 : 0), (1 : 1 : 0), (7/4 : -1 : 1), (7/4 : 2 : 1)]

        ::

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: f = DynamicalSystem_projective([x^2 - 21/16*z^2, y^2 - 2*z^2, z^2])
            sage: f.periodic_points(2)                                                  # needs sage.rings.function_field
            [(-5/4 : -1 : 1), (-5/4 : 2 : 1), (1/4 : -1 : 1), (1/4 : 2 : 1)]

        ::

            sage: set_verbose(None)
            sage: P.<x,y> = ProjectiveSpace(ZZ, 1)
            sage: f = DynamicalSystem_projective([x^2 + y^2, y^2])
            sage: f.periodic_points(2, R=QQbar, minimal=False)                          # needs sage.rings.number_field
            [(-0.50000000000000000? - 1.322875655532296?*I : 1),
             (-0.50000000000000000? + 1.322875655532296?*I : 1),
             (0.50000000000000000? - 0.866025403784439?*I : 1),
             (0.50000000000000000? + 0.866025403784439?*I : 1),
             (1 : 0)]

        ::

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: f = DynamicalSystem_projective([x^2 - 3/4*z^2, y^2 - 3/4*z^2, z^2])
            sage: f.periodic_points(2, formal=True)                                     # needs sage.rings.function_field
            [(-1/2 : -1/2 : 1), (-1/2 : 3/2 : 1), (3/2 : -1/2 : 1)]

        ::

            sage: P.<x,y> = ProjectiveSpace(GF(307), 1)
            sage: f = DynamicalSystem_projective([x^10 + y^10, y^10])
            sage: f.periodic_points(16, minimal=True, algorithm='cyclegraph')           # needs sage.graphs
            [(69 : 1), (185 : 1), (120 : 1), (136 : 1), (97 : 1), (183 : 1),
             (170 : 1), (105 : 1), (274 : 1), (275 : 1), (154 : 1), (156 : 1),
             (87 : 1), (95 : 1), (161 : 1), (128 : 1)]

        ::

            sage: # needs sage.rings.finite_rings
            sage: P.<x,y> = ProjectiveSpace(GF(13^2, 't'), 1)
            sage: f = DynamicalSystem_projective([x^3 + 3*y^3, x^2*y])
            sage: f.periodic_points(30, minimal=True, algorithm='cyclegraph')           # needs sage.graphs
            [(t + 3 : 1), (6*t + 6 : 1), (7*t + 1 : 1), (2*t + 8 : 1),
             (3*t + 4 : 1), (10*t + 12 : 1), (8*t + 10 : 1), (5*t + 11 : 1),
             (7*t + 4 : 1), (4*t + 8 : 1), (9*t + 1 : 1), (2*t + 2 : 1),
             (11*t + 9 : 1), (5*t + 7 : 1), (t + 10 : 1), (12*t + 4 : 1),
             (7*t + 12 : 1), (6*t + 8 : 1), (11*t + 10 : 1), (10*t + 7 : 1),
             (3*t + 9 : 1), (5*t + 5 : 1), (8*t + 3 : 1), (6*t + 11 : 1),
             (9*t + 12 : 1), (4*t + 10 : 1), (11*t + 4 : 1), (2*t + 7 : 1),
             (8*t + 12 : 1), (12*t + 11 : 1)]

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSystem_projective([3*x^2 + 5*y^2, y^2])
            sage: f.periodic_points(2, R=GF(3), minimal=False)                          # needs sage.rings.function_field
            [(2 : 1)]
            sage: f.periodic_points(2, R=GF(7))                                         # needs sage.rings.function_field
            []

        ::

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: f = DynamicalSystem_projective([x^2, x*y, z^2])
            sage: f.periodic_points(1)                                                  # needs sage.rings.function_field
            Traceback (most recent call last):
            ...
            TypeError: use return_scheme=True

        ::

            sage: # needs sage.rings.number_field
            sage: R.<x> = QQ[]
            sage: K.<u> = NumberField(x^2 - x + 3)
            sage: P.<x,y,z> = ProjectiveSpace(K, 2)
            sage: X = P.subscheme(2*x - y)
            sage: f = DynamicalSystem_projective([x^2 - y^2, 2*(x^2 - y^2), y^2 - z^2],
            ....:                                domain=X)
            sage: f.periodic_points(2)                                                  # needs sage.rings.function_field
            [(-1/5*u - 1/5 : -2/5*u - 2/5 : 1), (1/5*u - 2/5 : 2/5*u - 4/5 : 1)]

        ::

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: f = DynamicalSystem_projective([x^2 - y^2, x^2 - z^2, y^2 - z^2])
            sage: f.periodic_points(1)                                                  # needs sage.rings.function_field
            [(-1 : 0 : 1)]
            sage: f.periodic_points(1, return_scheme=True)
            Closed subscheme of Projective Space of dimension 2 over Rational Field
            defined by:
              -x^3 + x^2*y - y^3 + x*z^2,
              -x*y^2 + x^2*z - y^2*z + x*z^2,
              -y^3 + x^2*z + y*z^2 - z^3

        ::

            sage: P.<x,y>=ProjectiveSpace(GF(3), 1)
            sage: f = DynamicalSystem_projective([x^2 - 2*y^2, y^2])
            sage: f.periodic_points(2, R=GF(3^2,'t'))                                   # needs sage.rings.finite_rings
            [(t + 2 : 1), (2*t : 1)]

        ::

            sage: S.<c> = QQ[]
            sage: R.<x,y> = PolynomialRing(S, 2)
            sage: P = ProjectiveSpace(R)
            sage: f = DynamicalSystem_projective([x^2 + c*y^2, y^2])
            sage: f.periodic_points(2, return_scheme=True)                              # needs sage.rings.function_field
            Closed subscheme of Projective Space of dimension 1 over Univariate
            Polynomial Ring in c over Rational Field defined by:
              x^2 + x*y + (c + 1)*y^2

        ::

            sage: P.<x,y,z> = ProjectiveSpace(ZZ, 2)
            sage: f = DynamicalSystem([x^2 - 2*y^2, y^2, z^2])
            sage: X = f.periodic_points(2, minimal=False, formal=True,  # long time
            ....:                       return_scheme=True)
            sage: len(X.defining_polynomials())                         # long time
            19

        TESTS::

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: f = DynamicalSystem_projective([x^2 - y^2, 2*(x^2 - y^2), y^2 - z^2])
            sage: f.periodic_points(2, minimal=True)                                    # needs sage.rings.function_field
            Traceback (most recent call last):
            ...
            ValueError: dynamical system is not a morphism, cannot calculate minimal or formal periodic points
        """
        if n <= 0:
            raise ValueError("a positive integer period must be specified")
        if R is None:
            f_sub = self
            R = self.base_ring()
        else:
            f_sub = self.change_ring(R)
            R = f_sub.base_ring() #in the case when R is an embedding
        if isinstance(R, FractionField_1poly_field) or R in FunctionFields():
            raise NotImplementedError('periodic points not implemented for fraction function fields; '
                'clear denominators and use the polynomial ring instead')
        if isinstance(R, FractionField_generic):
            if isinstance(R.ring(), MPolynomialRing_base):
                raise NotImplementedError('periodic points not implemented for fraction function fields; '
                    'clear denominators and use the polynomial ring instead')
        CR = f_sub.coordinate_ring()
        dom = f_sub.domain()
        PS = f_sub.codomain().ambient_space()
        N = PS.dimension_relative() + 1
        FF = FractionField(R)
        if dom != PS:
            f = DynamicalSystem(f_sub.defining_polynomials())
        else:
            f = f_sub
        if algorithm == 'cyclegraph':
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
                return points
            else:
                raise TypeError("ring must be finite to generate cyclegraph")
        elif algorithm == 'variety':
            if formal and N == 2 and dom == PS:
                X = PS.subscheme([f.dynatomic_polynomial(n)])
            else:
                F = f.nth_iterate_map(n)
                L = [F[i]*CR.gen(j) - F[j]*CR.gen(i) for i in range(N)
                    for j in range(i+1, N)]
                L = [t for t in L if t != 0]
                X = PS.subscheme(L)
                if (minimal or formal) and n != 1:
                    if not f_sub.is_morphism():
                        raise ValueError('dynamical system is not a morphism, cannot calculate minimal or formal periodic points')
                    if formal:
                        if N > 2:
                            d = f.degree()
                            # we need a model with no periodic points at infinity
                            new_f, mat = f.affine_preperiodic_model(0, n, return_conjugation=True)
                            new_f.normalize_coordinates()
                            # we now deform by a parameter t
                            T = R['t']
                            t = T.gens()[0]
                            Pt = ProjectiveSpace(N-1, R=T, names=[str(i) for i in CR.gens()])
                            deformed_polys = [poly + t*Pt.gens()[-1]**d for poly in new_f.defining_polynomials()[:-1]]
                            deformed_polys += [new_f.defining_polynomials()[-1]]
                            f_deformed = DynamicalSystem(deformed_polys)

                            # after deforming by the parameter, the preperiodic points with multiplicity
                            # will separate into different points. we can now calculate the minimal preperiodic
                            # points with the parameter, and then specialize to get the formal periodic points
                            ideal = f_deformed.periodic_points(n, return_scheme=True).defining_ideal()
                            L = [poly.specialization({t:0}) for poly in ideal.gens()]
                            subs_list = mat.inverse()*vector(CR.gens())
                            subs = {}
                            for i in range(len(subs_list)):
                                subs[PS.gens()[i]] = subs_list[i]
                            if R.is_field():
                                X = PS.subscheme([poly.subs(subs) for poly in L])
                            else:
                                K = [poly.subs(subs) for poly in L]
                                K = [poly*poly.denominator() for poly in K]
                                X = PS.subscheme(K)
                    elif minimal:
                        Sn = []
                        for k in ZZ(n).divisors():
                            if ZZ(n//k).is_prime():
                                Sn.append(k)
                        if isinstance(R, (PolynomialRing_generic,
                                          MPolynomialRing_base)):
                            phi = FlatteningMorphism(CR)
                            flatCR = phi.codomain()
                            Ik = flatCR.ideal(1)
                            for k in Sn:
                                Ik *= f.periodic_points(k, return_scheme=True, minimal=False).defining_ideal()
                            psi = UnflatteningMorphism(flatCR, CR)
                            In = flatCR.ideal([phi(i) for i in X.defining_polynomials()])
                            X = PS.subscheme([psi(i) for i in In.saturation(Ik)[0].gens()])
                        else:
                            Ik = CR.ideal(1)
                            for k in Sn:
                                Ik *= f.periodic_points(k, return_scheme=True, minimal=False).defining_ideal()
                            In = X.defining_ideal()
                            X = PS.subscheme(In.saturation(Ik)[0])
            if dom != PS:
                X = PS.subscheme(list(X.defining_polynomials()) + list(dom.defining_polynomials()))
            if return_scheme:  # this includes the indeterminacy locus points!
                return X
            if X.change_ring(FF).dimension() <= 0:
                if R in NumberFields() or isinstance(R, sage.rings.abc.AlgebraicField) or R in FiniteFields():
                    Z = f.base_indeterminacy_locus()
                    points = [dom(Q) for Q in X.rational_points()]
                    good_points = []
                    for Q in points:
                        try:
                            Z(list(Q))
                        except TypeError:
                            good_points.append(Q)
                    good_points.sort()
                    return good_points
                else:
                    raise NotImplementedError("ring must be a number field or finite field")
            else: #a higher dimensional scheme
                raise TypeError("use return_scheme=True")
        else:
            raise ValueError("algorithm must be either 'variety' or 'cyclegraph'")

    def multiplier_spectra(self, n, formal=False, type='point', use_algebraic_closure=True, check=True):
        r"""
        Compute the ``n`` multiplier spectra of this dynamical system.

        This is the set of multipliers of all peroidic points of
        period ``n`` included with the appropriate multiplicity.
        User can also specify to compute the formal ``n`` multiplier spectra
        instead which includes the multipliers of all formal periodic points
        of period ``n`` with appropriate multiplicity. The map must be defined over
        projective space over a number field or finite field.

        By default, the computations are done over the algebraic closure of the
        base field. If the map is defined over projective space of dimension 1,
        the computation can be done over the minimal extension of the base field that
        contains the periodic points. Otherwise, it will be done over the base ring
        of the map.

        INPUT:

        - ``n`` -- positive integer, the period

        - ``formal`` -- boolean (default: ``False``); ``True`` specifies
          to find the formal ``n`` multiplier spectra of this map and
          ``False`` specifies to find the ``n`` multiplier spectra

        - ``type`` -- (default: ``'point'``) string; either ``'point'``
          or ``'cycle'`` depending on whether you compute one multiplier
          per point or one per cycle

        - ``use_algebraic_closure`` -- boolean (default: ``True``); if ``True`` uses the
          algebraic closure. Using the algebraic closure can sometimes lead to numerical instability
          and extraneous errors. For most accurate results in dimension 1, set to ``False``.
          If ``False``, and the map is defined over projective space of
          dimension 1, uses the smallest extension of the base field
          containing all the periodic points. If the map is defined over projective space
          of dimension greater than 1, then the base ring of the map is used.

        - ``check`` -- boolean (default: ``True``); whether to check if the
          full multiplier spectra was computed. If ``False``, can lead to
          mathematically incorrect answers in dimension greater than 1. Ignored
          if ``use_algebraic_closure`` is ``True`` or if this dynamical system is defined
          over projective space of dimension 1.

        OUTPUT:

        A list of field elements if the domain of the map is projective space of
        dimension 1. If the domain of the map is projective space of dimension
        greater than 1, a list of matrices

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSystem_projective([x^2 - 3/4*y^2, y^2])
            sage: sorted(f.multiplier_spectra(2, type='point'))                         # needs sage.rings.number_field
            [0, 1, 1, 1, 9]
            sage: sorted(f.multiplier_spectra(2, type='cycle'))                         # needs sage.rings.number_field
            [0, 1, 1, 9]

        ::

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: f = DynamicalSystem_projective([x^2, z^2, y^2])
            sage: f.multiplier_spectra(1)                                               # needs sage.rings.number_field
            [
            [                       2 1 - 1.732050807568878?*I]
            [                       0                       -2],
            [                       2 1 + 1.732050807568878?*I]  [ 0  0]  [ 0  0]
            [                       0                       -2], [ 0 -2], [ 0 -2],
            [ 0  0]  [0 0]  [ 2 -2]
            [ 0 -2], [0 0], [ 0 -2]
            ]

        ::

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: f = DynamicalSystem_projective([x^2, z^2, y^2])
            sage: f.multiplier_spectra(2, formal=True)  # long time
            [
            [4 0]  [4 0]  [4 0]  [4 0]  [4 0]  [4 0]  [4 0]  [4 0]  [0 0]  [0 0]
            [0 4], [0 0], [0 0], [0 4], [0 4], [0 0], [0 0], [0 4], [0 0], [0 0],
            [4 0]  [4 0]  [4 0]  [4 0]
            [0 4], [0 4], [0 0], [0 0]
            ]

        ::

            sage: # needs sage.rings.number_field
            sage: set_verbose(None)
            sage: z = QQ['z'].0
            sage: K.<w> = NumberField(z^4 - 4*z^2 + 1,'z')
            sage: P.<x,y> = ProjectiveSpace(K, 1)
            sage: f = DynamicalSystem_projective([x^2 - w/4*y^2, y^2])
            sage: sorted(f.multiplier_spectra(2, formal=False, type='cycle'))
            [0,
             0.0681483474218635? - 1.930649271699173?*I,
             0.0681483474218635? + 1.930649271699173?*I,
             5.931851652578137? + 0.?e-49*I]

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSystem_projective([4608*x^10 - 2910096*x^9*y + 325988068*x^8*y^2
            ....:         + 31825198932*x^7*y^3 - 4139806626613*x^6*y^4 - 44439736715486*x^5*y^5
            ....:         + 2317935971590902*x^4*y^6 - 15344764859590852*x^3*y^7
            ....:         + 2561851642765275*x^2*y^8 + 113578270285012470*x*y^9
            ....:         - 150049940203963800*y^10, 4608*y^10])
            sage: sorted(f.multiplier_spectra(1))                                       # needs sage.rings.number_field
            [-119820502365680843999,
             -7198147681176255644585/256,
             -3086380435599991/9,
             -3323781962860268721722583135/35184372088832,
             -4290991994944936653/2097152,
             0,
             529278480109921/256,
             1061953534167447403/19683,
             848446157556848459363/19683,
             82911372672808161930567/8192,
             3553497751559301575157261317/8192]

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSystem_projective([x^2 - 7/4*y^2, y^2])
            sage: f.multiplier_spectra(3, formal=True, type='cycle')                    # needs sage.rings.number_field
            [1, 1]
            sage: f.multiplier_spectra(3, formal=True, type='point')                    # needs sage.rings.number_field
            [1, 1, 1, 1, 1, 1]

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSystem_projective([x^4 + 3*y^4, 4*x^2*y^2])
            sage: f.multiplier_spectra(1, use_algebraic_closure=False)                  # needs sage.rings.number_field
            [0,
             -1,
             1/128*a^5 - 13/384*a^4 + 5/96*a^3 + 1/16*a^2 + 43/128*a + 303/128,
             -1/288*a^5 + 1/96*a^4 + 1/24*a^3 - 1/3*a^2 + 5/32*a - 115/32,
             -5/1152*a^5 + 3/128*a^4 - 3/32*a^3 + 13/48*a^2 - 63/128*a - 227/128]
            sage: f.multiplier_spectra(1)                                               # needs sage.rings.number_field
            [0,
             -1,
             1.951373035591442?,
             -2.475686517795721? - 0.730035681602057?*I,
             -2.475686517795721? + 0.730035681602057?*I]

        ::

            sage: P.<x,y> = ProjectiveSpace(GF(5), 1)
            sage: f = DynamicalSystem_projective([x^4 + 2*y^4, 4*x^2*y^2])
            sage: f.multiplier_spectra(1, use_algebraic_closure=False)                  # needs sage.rings.finite_rings
            [0, 3*a + 3, 2*a + 1, 1, 1]
            sage: f.multiplier_spectra(1)
            [0, 2*z2 + 1, 3*z2 + 3, 1, 1]

        ::

            sage: # needs sage.rings.number_field
            sage: P.<x,y> = ProjectiveSpace(QQbar, 1)
            sage: f = DynamicalSystem_projective([x^5 + 3*y^5, 4*x^3*y^2])
            sage: f.multiplier_spectra(1)
            [0,
             -4.106544657178796?,
             -7/4,
             1.985176555073911?,
             -3.064315948947558? - 1.150478041113253?*I,
             -3.064315948947558? + 1.150478041113253?*I]

        ::

            sage: K = GF(3).algebraic_closure()
            sage: P.<x,y> = ProjectiveSpace(K, 1)
            sage: f = DynamicalSystem_projective([x^5 + 2*y^5, 4*x^3*y^2])
            sage: f.multiplier_spectra(1)
            [0, z3 + 2, z3 + 1, z3, 1, 1]

        TESTS::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSystem_projective([x^2 + y^2, x*y])
            sage: f.multiplier_spectra(1)                                               # needs sage.rings.number_field
            [1, 1, 1]

        ::

            sage: K = GF(3).algebraic_closure()
            sage: P.<x,y,z> = ProjectiveSpace(K, 2)
            sage: f = DynamicalSystem_projective([x^2 + 2*y^2, 4*x*y, z^2])
            sage: f.multiplier_spectra(1)                                               # needs sage.rings.number_field
            [
            [0 0]  [1 0]  [1 0]  [1 0]  [2 0]  [2 0]  [2 0]
            [0 0], [0 0], [0 0], [0 0], [0 1], [0 1], [0 1]
            ]

        ::

            sage: F.<a> = GF(7)
            sage: P.<x,y> = ProjectiveSpace(F, 1)
            sage: f = DynamicalSystem_projective([x^2 + y^2, y^2])
            sage: sorted(f.multiplier_spectra(1))
            [0, 3, 6]

        ::

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: f = DynamicalSystem_projective([x^2, z^2, y^2])
            sage: g = f.change_ring(QQbar)                                              # needs sage.rings.number_field
            sage: f.multiplier_spectra(1) == g.multiplier_spectra(1)    # long time, needs sage.rings.number_field
            True

        ::

            sage: # needs sage.rings.number_field
            sage: K.<w> = QuadraticField(5)
            sage: P.<x,y,z> = ProjectiveSpace(K, 2)
            sage: f = DynamicalSystem_projective([x^2 + w*x*y + y^2, y^2, z^2])
            sage: f.multiplier_spectra(1)                               # long time
            [
            [1.000000000000000? - 1.572302755514847?*I                                         0]
            [1.000000000000000? - 1.572302755514847?*I 0.618033988749895? - 1.757887921270715?*I]
            [1.000000000000000? + 1.572302755514847?*I                                         0]
            [1.000000000000000? + 1.572302755514847?*I 0.618033988749895? + 1.757887921270715?*I]
            [                                        0                                         0],
            [                                        0                                         2],
            [                                        0                                         0],
            [                                        0                                         2],
            [0 0]  [0 0]  [                 2 2.236067977499790?]
            [0 0], [0 0], [                 0                  0]
            ]

        ::

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: f = DynamicalSystem_projective([x^2, z^2, y^2])
            sage: f.multiplier_spectra(1, use_algebraic_closure=False)                  # needs sage.rings.function_field
            Traceback (most recent call last):
            ...
            ValueError: failed to compute the full multiplier spectra. Try use_algebraic_closure=True
            or extend the base ring of this dynamical system
        """
        PS = self.domain()
        n = Integer(n)

        if (n < 1):
            raise ValueError("period must be a positive integer")
        if not isinstance(PS, ProjectiveSpace_ring):
            raise NotImplementedError("not implemented for subschemes")

        if PS.dimension_relative() > 1:
            K = self.domain().base_ring()

            # if we are already using an algebraic closure, we move the
            # map into a finite extension and set use_algebraic_closure to True
            # in order to get a scheme defined over a finite extension
            if isinstance(K, (sage.rings.abc.AlgebraicField,
                              AlgebraicClosureFiniteField_generic)):
                f = self.reduce_base_field()
                K = f.base_ring()
                use_algebraic_closure = True
            else:
                f = self

            # in order to calculate multiplicity, we need to have a scheme defined
            # over a finite extension, not an algebraic closure
            X = f.periodic_points(n, minimal=False, formal=formal, return_scheme=True)
            if use_algebraic_closure:
                number_field = False
                finite_field = False
                if K in NumberFields():
                    number_field = True
                if K in FiniteFields():
                    finite_field = True
                if not (number_field or finite_field):
                    raise NotImplementedError('Only implemented for number fields, QQbar, finite fields, and algebraic closures of finite fields')
                Kbar = K.algebraic_closure()
                if Kbar.has_coerce_map_from(K):
                    f = f.change_ring(Kbar)
                    rat_points = X.rational_points(F=Kbar)
                else:
                    embeds = K.embeddings(Kbar)
                    if embeds:
                        X2 = X.change_ring(embeds[0])
                        rat_points = X2.rational_points()
                        f = self.change_ring(embeds[0])
                    else:
                        raise ValueError("no embeddings of base field to algebraic closure")
            else:
                rat_points = X.rational_points()
                f = self
            PS = f.domain()
            points = []
            for point in rat_points:
                if use_algebraic_closure:
                    if number_field:
                        # in order to calculate multiplicity, the point must be defined over a finite extension
                        K2, pnt_lst, _ = number_field_elements_from_algebraics(list(point))
                        # we coerce if we can
                        if K.has_coerce_map_from(K2):
                            for i in range(X.multiplicity(pnt_lst)):
                                points.append(PS(point))
                        elif K2.has_coerce_map_from(K):
                            X_k = X.change_ring(K2)
                            for i in range(X_k.multiplicity(pnt_lst)):
                                points.append(PS(point))
                        # otherwise, we need to calculate a composite field
                        else:
                            _, K_embed, K2_embed, _ = K.composite_fields(K2, both_maps=True)[0]
                            X_k = X.change_ring(K_embed)
                            pnt_lst = [K2_embed(pnt) for pnt in pnt_lst]
                            new_point = X_k.ambient_space()(pnt_lst)
                            for i in range(X_k.multiplicity(new_point)):
                                points.append(PS(point))
                    else:
                        # we find a finite extension which the current point
                        # and X coerce into
                        final_degree = K.degree()
                        new_point = []
                        for num in list(point):
                            ff_num = num.as_finite_field_element()[1]
                            new_point.append(ff_num)
                            degree = ff_num.parent().degree()
                            final_degree = final_degree.lcm(degree)
                        K_prime = GF(K.characteristic()**final_degree)
                        X_k = X.change_ring(K_prime)
                        for i in range(X_k.multiplicity(new_point)):
                            points.append(PS(point))
                else:
                    for i in range(X.multiplicity(point)):
                        points.append(point)
            if not use_algebraic_closure:
                if check:
                    # we check if we computed the full multiplier spectra
                    d = self.degree()
                    N = self.domain().ambient_space().dimension_relative()
                    if not formal:
                        expected_number = sum(d**(n*i) for i in range(N+1))
                    else:
                        expected_number = 0
                        for D in n.divisors():
                            u = moebius(n/D)
                            inner_sum = sum(d**(D*j) for j in range(N+1))
                            expected_number += u*inner_sum
                    if len(points) != expected_number:
                        raise ValueError('failed to compute the full multiplier spectra. Try use_algebraic_closure=True'
                         + ' or extend the base ring of this dynamical system')
        else:
            K = FractionField(self.codomain().base_ring())
            if use_algebraic_closure:
                Kbar = K.algebraic_closure()
                if Kbar.has_coerce_map_from(K):
                    f = self.change_ring(Kbar)
                else:
                    embeds = K.embeddings(Kbar)
                    if embeds:
                        f = self.change_ring(embeds[0])
                    else:
                        raise ValueError("no embeddings of base field to algebraic closure")
            else:
                embedding = self.field_of_definition_periodic(n, formal=formal, return_embedding=True)[1]
                f = self.change_ring(embedding)

            PS = f.domain()
            if not formal:
                G = f.nth_iterate_map(n)
                F = G[0]*PS.gens()[1] - G[1]*PS.gens()[0]
            else:
                # periodic points of formal period n are the roots of the nth dynatomic polynomial
                F = f.dynatomic_polynomial(n)

            other_roots = F.parent()(F([(f.domain().gens()[0]),1])).univariate_polynomial().roots(ring=f.base_ring())

            points = []

            minfty = min(ex[1] for ex in F.exponents()) # include the point at infinity with the right multiplicity
            for i in range(minfty):
                points.append(PS([1,0]))

            for R in other_roots:
                for i in range(R[1]):
                    points.append(PS([R[0],1])) # include copies of higher multiplicity roots

        if type == 'cycle':
            # should include one representative point per cycle, included with the right multiplicity
            newpoints = []

            while points:
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
            points = newpoints

        if PS.dimension_relative() > 1:
            multipliers = [f.multiplier(pt,n) for pt in points]
        else:
            multipliers = [f.multiplier(pt,n)[0,0] for pt in points]

        return multipliers

    def sigma_invariants(self, n, formal=False, embedding=None, type='point',
                        return_polynomial=False, chow=False, deform=False, check=True):
        r"""
        Compute the values of the elementary symmetric polynomials evaluated
        on the ``n`` multiplier spectra of this dynamical system.

        The sigma invariants are the symmetric polynomials evaluated on the
        characteristic polynomial of the multipliers. See [Hutz2019]_ for
        the full definition. Spepcifically, this function returns either
        the following polynomial or its coefficients (with signs
        appropriately adjusted):

         .. MATH::

            \prod_{P \text{ period n}} ( w - c(P,t)),

        where `c(P,t)` is the charateristic polynomial (variable `t`) of the
        multiplier at `P`. Note that in dimension 1, only the coefficients
        of the constant term is returned.

        The invariants can be computed for points of period ``n`` or
        points of formal period ``n``. The base
        ring should be a number field, number field order, or
        a finite field or a polynomial ring or function field over a
        number field, number field order, or finite field.

        The parameter ``type`` determines if the sigma are computed from
        the multipliers calculated at one per cycle (with multiplicity)
        or one per point (with multiplicity). Only implemented
        for dimension 1. Note that in the ``cycle`` case, a map with a cycle
        which collapses into multiple smaller cycles, this is still
        considered one cycle. In other words, if a 4-cycle collapses into
        a 2-cycle with multiplicity 2, there is only one multiplier used
        for the doubled 2-cycle when computing ``n=4``.

        ALGORITHM:

        In dimension 1, we use the Poisson product of the resultant of
        two polynomials:

        .. MATH::

            res(f,g) = \prod_{f(a)=0} g(a).

        In higher dimensions, we use elimination theory (Groebner bases)
        to compute the equivalent of the Poisson product. Letting `f` be
        the polynomial defining the periodic or formal
        periodic points and `g` the polynomial `w - F` for an auxilarly
        variable `w` and `F` the characteristic polynomial of the Jacobian matrix
        of `f`. Note that if `f` is a rational function, we clear
        denominators for `g`.

        To calculate the full polynomial defining the sigma invariants,
        we follow the algorithm outlined in section 4 of [Hutz2019]_. There
        are 4 cases:

        - multipliers and ``n`` periodic points all distinct -- in this case,
          we can use Proposition 4.1 of [Hutz2019]_ to compute the sigma invariants

        - ``n`` -- periodic points are all distinct, multipliers are repeated; here we
          can use Proposition 4.2 of [Hutz2019]_ to compute the sigma invariants.
          This corresponds to ``chow=True``.

        - ``n`` -- periodic points are repeated, multipliers are all distinct; to deal
          with this case, we deform the map by a formal parameter `k`. The deformation
          separates the ``n`` periodic points, making them distinct, and we can recover
          the ``n`` periodic points of the original map by specializing `k` to 0.
          This corresponds to ``deform=True``.

        - ``n`` -- periodic points are repeated, multipliers are repeated; here we
          can use both cases 2 and 3 together. This corresponds to ``deform=True``
          and ``chow=True``.

        As we do not want to check which case we are in beforehand, we throw a
        :exc:`ValueError` if the computed polynomial does not have the correct
        degree.

        INPUT:

        - ``n`` -- positive integer, the period

        - ``formal`` -- boolean (default: ``False``); ``True`` specifies
          to find the values of the elementary symmetric polynomials
          corresponding to the formal ``n`` multiplier spectra and ``False``
          specifies to instead find the values corresponding to the ``n``
          multiplier spectra, which includes the multipliers of all
          periodic points of period ``n``

        - ``embedding`` -- (default: ``None``) must be ``None``, passing an embedding
          is no longer supported, see :issue:`32205`

        - ``type`` -- (default: ``'point'``) string; either ``'point'``
          or ``'cycle'`` depending on whether you compute with one
          multiplier per point or one per cycle. Not implemented for
          dimension greater than 1.

        - ``return polynomial`` -- boolean (default: ``False``);
          ``True`` specifies returning the polynomial which generates
          the sigma invariants, see [Hutz2019]_ for the full definition.
          The polynomial is always a multivariate polynomial with variables
          ``w`` and ``t``.

        - ``chow`` -- boolean (default: ``False``); ``True`` specifies
          using the Chow algorithm from [Hutz2019]_ to compute the sigma
          invariants. While slower, the Chow algorithm does not lose
          information about multiplicities of the multipliers. In order
          to accurately compute the sigma polynomial when there is a
          repeated multiplier, ``chow`` must be ``True``.

        - ``deform`` -- boolean (default: ``False``); ``True`` specifies
          first deforming the map so that all periodic points are distinct
          and then calculating the sigma invariants. In order to accurately
          calculate the sigma polynomial when there is a periodic point with
          multiplicity, ``deform`` must be ``True``.

        - ``check`` -- boolean (default: ``True``); when ``True`` the degree of
          the sigma polynomial is checked against the expected degree. This is
          done as the sigma polynomial may drop degree if multiplicities of periodic
          points or multipliers are not correctly accounted for using ``chow`` or
          ``deform``.

        .. WARNING::

            Setting ``check`` to ``False`` can lead to mathematically incorrect
            answers.

        OUTPUT: list of elements in the base ring, unless ``return_polynomial``
        is ``True``, in which case a polynomial in ``w`` and ``t`` is returned.
        The variable ``t`` is the variable of the characteristic
        polynomials of the multipliers.

        If this map is defined over `\mathbb{P}^N`, where `N > 1`, then
        the list is the coefficients of `w` and `t`, in lexographical order with `w > t`.

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSystem_projective([x^2 + x*y + y^2, y^2 + x*y])
            sage: f.sigma_invariants(1)                                                 # needs sage.rings.number_field
            [3, 3, 1]

        If ``return_polynomial`` is ``True``, then following [Hutz2019]_
        we return a two variable polynomial in `w` and `t`::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSystem_projective([x^2 + 2*y^2, y^2])
            sage: poly = f.sigma_invariants(1, return_polynomial=True); poly
            w^3 - 3*w^2*t + 2*w^2 + 3*w*t^2 - 4*w*t + 8*w - t^3 + 2*t^2 - 8*t

        From the full polynomial, we can easily recover the one variable polynomial whose coefficients
        are symmetric functions in the multipliers, up to sign::

            sage: w, t = poly.variables()
            sage: poly.specialization({w:0}).monic()
            t^3 - 2*t^2 + 8*t
            sage: f.sigma_invariants(1)                                                 # needs sage.rings.number_field
            [2, 8, 0]

        For dynamical systems on `\mathbb{P}^N`, where `N > 1`, the full polynomial
        is needed to distinguish the conjugacy class. We can, however, still return
        a list in this case::

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: f = DynamicalSystem_projective([x^2, z^2, y^2])
            sage: f.sigma_invariants(1, chow=True)
            [1, 7, -6, -12, 21, -36, -60, 72, 48, 35, -90, -120, 352, 96, -288, -64, 35,
             -120, -120, 688, -96, -1056, 320, 384, 0, 21, -90, -60, 672, -384, -1440, 1344,
             768, -768, 0, 0, 7, -36, -12, 328, -336, -864, 1472, 384, -1536, 512, 0, 0, 0,
             1, -6, 0, 64, -96, -192, 512, 0, -768, 512, 0, 0, 0, 0, 0]

        When calculating the sigma invariants for `\mathbb{P}^N`, with `N > 1`,
        the default algorithm loses information about multiplicities. Note that
        the following call to sigma invariants returns a degree 6 polynomial in `w`::

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: f = DynamicalSystem_projective([x^2, y^2, z^2])
            sage: f.sigma_invariants(1, return_polynomial=True, check=False)
            w^6 - 6*w^5*t^2 + 8*w^5*t - 4*w^5 + 15*w^4*t^4 - 40*w^4*t^3 + 40*w^4*t^2 -
            16*w^4*t - 20*w^3*t^6 + 80*w^3*t^5 - 120*w^3*t^4 + 80*w^3*t^3 - 16*w^3*t^2 +
            15*w^2*t^8 - 80*w^2*t^7 + 160*w^2*t^6 - 144*w^2*t^5 + 48*w^2*t^4 - 6*w*t^10 +
            40*w*t^9 - 100*w*t^8 + 112*w*t^7 - 48*w*t^6 + t^12 - 8*t^11 + 24*t^10 -
            32*t^9 + 16*t^8

        Setting ``chow`` to ``True``, while much slower, accounts correctly for multiplicities.
        Note that the following returns a degree 7 polynomial in `w`::

            sage: f.sigma_invariants(1, return_polynomial=True, chow=True)
            w^7 - 7*w^6*t^2 + 10*w^6*t - 4*w^6 + 21*w^5*t^4 - 60*w^5*t^3 + 60*w^5*t^2 -
            24*w^5*t - 35*w^4*t^6 + 150*w^4*t^5 - 240*w^4*t^4 + 176*w^4*t^3 - 48*w^4*t^2 +
            35*w^3*t^8 - 200*w^3*t^7 + 440*w^3*t^6 - 464*w^3*t^5 + 224*w^3*t^4 -
            32*w^3*t^3 - 21*w^2*t^10 + 150*w^2*t^9 - 420*w^2*t^8 + 576*w^2*t^7 -
            384*w^2*t^6 + 96*w^2*t^5 + 7*w*t^12 - 60*w*t^11 + 204*w*t^10 - 344*w*t^9 +
            288*w*t^8 - 96*w*t^7 - t^14 + 10*t^13 - 40*t^12 + 80*t^11 - 80*t^10 + 32*t^9

        ::

            sage: # needs sage.rings.number_field
            sage: set_verbose(None)
            sage: z = QQ['z'].0
            sage: K = NumberField(z^4 - 4*z^2 + 1, 'z')
            sage: P.<x,y> = ProjectiveSpace(K, 1)
            sage: f = DynamicalSystem_projective([x^2 - 5/4*y^2, y^2])
            sage: f.sigma_invariants(2, formal=False, type='cycle')
            [13, 11, -25, 0]
            sage: f.sigma_invariants(2, formal=False, type='point')
            [12, -2, -36, 25, 0]

        check that infinity as part of a longer cycle is handled correctly::

            sage: # needs sage.rings.number_field
            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSystem_projective([y^2, x^2])
            sage: f.sigma_invariants(2, type='cycle')
            [12, 48, 64, 0]
            sage: f.sigma_invariants(2, type='point')
            [12, 48, 64, 0, 0]
            sage: f.sigma_invariants(2, type='cycle', formal=True)
            [0]
            sage: f.sigma_invariants(2, type='point', formal=True)
            [0, 0]

        ::

            sage: # needs sage.rings.number_field
            sage: K.<w> = QuadraticField(3)
            sage: P.<x,y> = ProjectiveSpace(K, 1)
            sage: f = DynamicalSystem_projective([x^2 - w*y^2, (1-w)*x*y])
            sage: f.sigma_invariants(2, formal=False, type='cycle')
            [6*w + 21, 78*w + 159, 210*w + 367, 90*w + 156]
            sage: f.sigma_invariants(2, formal=False, type='point')
            [6*w + 24, 96*w + 222, 444*w + 844, 720*w + 1257, 270*w + 468]

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSystem_projective([512*x^5 - 378128*x^4*y + 76594292*x^3*y^2
            ....:         - 4570550136*x^2*y^3 - 2630045017*x*y^4 + 28193217129*y^5, 512*y^5])
            sage: f.sigma_invariants(1)                                                 # needs sage.rings.number_field
            [19575526074450617/1048576, -9078122048145044298567432325/2147483648,
             -2622661114909099878224381377917540931367/1099511627776,
             -2622661107937102104196133701280271632423/549755813888,
             338523204830161116503153209450763500631714178825448006778305/72057594037927936, 0]

        ::

            sage: P.<x,y,z> = ProjectiveSpace(GF(5), 2)
            sage: f = DynamicalSystem([x^2, y^2, z^2])
            sage: f.sigma_invariants(1, chow=True, return_polynomial=True)
            w^7 - 2*w^6*t^2 + w^6 + w^5*t^4 + w^5*t + w^4*t^3 + 2*w^4*t^2 + w^3*t^5 -
            w^3*t^4 - 2*w^3*t^3 - w^2*t^10 + w^2*t^7 + w^2*t^6 + w^2*t^5 + 2*w*t^12 -
            w*t^10 + w*t^9 - 2*w*t^8 - w*t^7 - t^14 + 2*t^9

        ::

            sage: # needs sage.rings.number_field
            sage: R.<c> = QQ[]
            sage: Pc.<x,y> = ProjectiveSpace(R, 1)
            sage: f = DynamicalSystem_projective([x^2 + c*y^2, y^2])
            sage: f.sigma_invariants(1)
            [2, 4*c, 0]
            sage: f.sigma_invariants(2, formal=True, type='point')
            [8*c + 8, 16*c^2 + 32*c + 16]
            sage: f.sigma_invariants(2, formal=True, type='cycle')
            [4*c + 4]

        ::

            sage: R.<c> = QQ[]
            sage: P.<x,y> = ProjectiveSpace(R, 1)
            sage: f = DynamicalSystem([x^2 + c*y^2, y^2])
            sage: f.sigma_invariants(1, return_polynomial=True)
            w^3 + (-3)*w^2*t + 2*w^2 + 3*w*t^2 + (-4)*w*t + 4*c*w - t^3 + 2*t^2 + (-4*c)*t
            sage: f.sigma_invariants(2, chow=True, formal=True,                         # needs sage.libs.pari
            ....:                    return_polynomial=True)
            w^2 + (-2)*w*t + (8*c + 8)*w + t^2 + (-8*c - 8)*t + 16*c^2 + 32*c + 16

        ::

            sage: R.<c,d> = QQ[]
            sage: P.<x,y,z> = ProjectiveSpace(R, 2)
            sage: f = DynamicalSystem([x^2 + c*z^2, y^2 + d*z^2, z^2])
            sage: len(dict(f.sigma_invariants(1, return_polynomial=True)))
            51

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSystem([x^2 + 3*y^2, x*y])
            sage: f.sigma_invariants(1, deform=True, return_polynomial=True)            # needs sage.rings.function_field
            w^3 - 3*w^2*t + 3*w^2 + 3*w*t^2 - 6*w*t + 3*w - t^3 + 3*t^2 - 3*t + 1

        doubled fixed point::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSystem_projective([x^2 - 3/4*y^2, y^2])
            sage: f.sigma_invariants(2, formal=True)                                    # needs sage.rings.number_field
            [2, 1]

        doubled 2 cycle::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSystem_projective([x^2 - 5/4*y^2, y^2])
            sage: f.sigma_invariants(4, formal=False, type='cycle')                     # needs sage.rings.number_field
            [170, 5195, 172700, 968615, 1439066, 638125, 0]

        TESTS::

            sage: F.<t> = FunctionField(GF(5))
            sage: P.<x,y> = ProjectiveSpace(F,1)
            sage: f = DynamicalSystem_projective([x^2 + (t/(t^2+1))*y^2, y^2], P)
            sage: f.sigma_invariants(1)                                                 # needs sage.rings.number_field
            [2, 4*t/(t^2 + 1), 0]

        ::

            sage: # needs sage.rings.number_field
            sage: R.<w> = QQ[]
            sage: N.<n> = NumberField(w^2 + 1)
            sage: P.<x,y,z> = ProjectiveSpace(N, 2)
            sage: f = DynamicalSystem_projective([x^2, y^2, z^2])
            sage: f.sigma_invariants(1, chow=True) == f.change_ring(QQ).sigma_invariants(1, chow=True)
            True

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSystem([x^2 + 3*y^2, x*y])
            sage: f.sigma_invariants(1, formal=True, return_polynomial=True)
            Traceback (most recent call last):
            ..
            ValueError: sigma polynomial dropped degree, as multiplicities were not accounted
            for correctly; try setting chow=True and/or deform=True

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSystem([x^2 + 3*y^2, x*y])
            sage: f.sigma_invariants(1, return_polynomial=True)
            Traceback (most recent call last):
            ..
            ValueError: sigma polynomial dropped degree, as multiplicities were not accounted
            for correctly; try setting chow=True and/or deform=True
        """
        n = ZZ(n)

        if embedding is not None:
            raise ValueError('do not specify an embedding')
        if n < 1:
            raise ValueError("period must be a positive integer")
        dom = self.domain()
        if not isinstance(dom, ProjectiveSpace_ring):
            raise NotImplementedError("not implemented for subschemes")
        if self.degree() <= 1:
            raise TypeError("must have degree at least 2")
        if type not in ['point', 'cycle']:
            raise ValueError("type must be either point or cycle")
        if dom.dimension_relative() > 1 or return_polynomial:
            if type == 'cycle':
                raise NotImplementedError('cycle not implemented for dimension greater than 1')
            base_ring = self.base_ring()
            d = self.degree()
            N = dom.dimension_relative()
            f = copy(self)
            Fn = f.nth_iterate_map(n)
            CR = f.codomain().ambient_space().coordinate_ring()
            if deform:
                # we need a model with all affine periodic points
                new_f = f.affine_preperiodic_model(0, n)
                new_f.normalize_coordinates()
                # we now deform by a parameter t
                T = base_ring['k']
                k = T.gens()[0]
                Pt = ProjectiveSpace(N, R=T, names=[str(i) for i in CR.gens()])
                deformed_polys = [poly + k*Pt.gens()[-1]**d for poly in new_f.defining_polynomials()[:-1]]
                deformed_polys += [new_f.defining_polynomials()[-1]]
                f_deformed = DynamicalSystem(deformed_polys)
                sigma_poly = f_deformed.sigma_invariants(n, chow=chow, deform=False, return_polynomial=True, check=False)
                sigma_polynomial = sigma_poly.specialization({k:0})
                # we fix the ordering of the parent polynomial ring
                new_parent = sigma_polynomial.parent().change_ring(order='lex')
                sigma_polynomial = new_parent(sigma_polynomial)
                sigma_polynomial *= sigma_polynomial.coefficients()[0].inverse_of_unit()
            else:
                if not base_ring.is_field():
                    F = FractionField(base_ring)
                    f.normalize_coordinates()
                    X = f.periodic_points(n, minimal=False, formal=formal, return_scheme=True)
                    X = X.change_ring(F)
                else:
                    F = base_ring
                    if isinstance(base_ring, FractionField_generic):
                        if isinstance(base_ring.ring(), (MPolynomialRing_base,
                                                         PolynomialRing_generic)):
                            f.normalize_coordinates()
                            f_ring = f.change_ring(base_ring.ring())
                            X = f_ring.periodic_points(n, minimal=False, formal=formal, return_scheme=True)
                            X = X.change_ring(F)
                    else:
                        X = f.periodic_points(n, minimal=False, formal=formal, return_scheme=True)
                newR = PolynomialRing(F, 'w, t', 2, order='lex')
                if not base_ring.is_field():
                    ringR = PolynomialRing(base_ring, 'w, t', 2, order='lex')
                if chow:
                    # create full polynomial ring
                    R = PolynomialRing(F, 'v', 2*N+3, order='lex')
                    var = list(R.gens())
                    # create polynomial ring for result
                    R2 = PolynomialRing(F, var[:N] + var[-2:])
                    psi = R2.hom(N*[0]+list(newR.gens()), newR)
                    # create substitution to set extra variables to 0
                    R_zero = {R.gen(N): 1}
                    for j in range(N+1, 2*N+1):
                        R_zero[R.gen(j)] = 0
                    t = var.pop()
                    w = var.pop()
                    var = var[:N]
                else:
                    R = PolynomialRing(F, 'v', N+2, order='lex')
                    psi = R.hom(N*[0] + list(newR.gens()), newR)
                    var = list(R.gens())
                    t = var.pop()
                    w = var.pop()
                sigma_polynomial = 1
                # go through each affine patch to avoid repeating periodic points
                # setting the visited coordinates to 0 as we go
                for j in range(N,-1,-1):
                    Xa = X.affine_patch(j)
                    fa = Fn.dehomogenize(j)
                    Pa = fa.domain()
                    Ra = Pa.coordinate_ring()
                    # create the images for the Hom to the ring we will do the elimination over
                    # with done affine patch coordinates as 0
                    if chow:
                        im = [R.gen(i) for i in range(j)] + (N-j)*[0] + [R.gen(i) for i in range(N, R.ngens())]
                    else:
                        im = list(R.gens())[:j] + (N-j)*[0] + [R.gen(i) for i in range(N, R.ngens())]
                    phi = Ra.hom(R.gens()[0:len(Ra.gens())])
                    # create polymomial that evaluates to the characteristic polynomial
                    M = t*matrix.identity(R, N)
                    g = (M-jacobian([phi(F.numerator())/phi(F.denominator()) for F in fa], var)).det()
                    # create the terms of the sigma invariants prod(w-lambda)
                    g_prime = w*R(g.denominator())(im)-R(g.numerator())(im)
                    # move the defining polynomials to the polynomial ring
                    L = [phi(h)(im) for h in Xa.defining_polynomials()]
                    # add the appropriate final polynomial to compute the sigma invariant polynomial
                    # via a Poisson product in elimination
                    if chow:
                        L += [g_prime + sum(R.gen(j-1)*R.gen(N+j)*(R(g.denominator())(im)) for j in range(1,N+1))]
                    else:
                        L += [g_prime]
                    I = R.ideal(L)
                    # since R is lex ordering, this is an elimination step
                    G = I.groebner_basis()
                    # the polynomial we need is the one just in w and t
                    if chow:
                        poly = psi(G[-1].specialization(R_zero))
                        if len(list(poly)) > 0:
                            poly *= poly.coefficients()[0].inverse_of_unit()
                    else:
                        poly = psi(G[-1])
                    if not base_ring.is_field():
                        denom = lcm([coeff[0].denominator() for coeff in poly])
                        poly *= denom
                    sigma_polynomial *= poly
                if not base_ring.is_field():
                    sigma_polynomial = ringR(sigma_polynomial)
            if check:
                degree_w = sigma_polynomial.degrees()[0]
                if formal:
                    expected_degree = 0
                    for D in n.divisors():
                        u = moebius(n/D)
                        inner_sum = sum(d**(D*j) for j in range(N+1))
                        expected_degree += u*inner_sum
                else:
                    expected_degree = sum(d**(n*i) for i in range(N+1))
                if degree_w != expected_degree:
                    raise ValueError('sigma polynomial dropped degree, as multiplicities were not accounted for correctly; '
                                     'try setting chow=True and/or deform=True')
            if return_polynomial:
                return sigma_polynomial
            # if we are returning a numerical list, read off the coefficients
            # in order of degree adjusting sign appropriately
            sigmas = []
            sigma_dictionary = dict([list(reversed(i)) for i in list(sigma_polynomial)])
            degree_w = sigma_polynomial.degrees()[0]
            w, t = sigma_polynomial.variables()
            for i in range(degree_w + 1):
                for j in range(2*i, -1, -1):
                    sigmas.append((-1)**(i+j)*sigma_dictionary.pop(w**(degree_w - i)*t**(j), 0))
            return sigmas

        base_ring = dom.base_ring()
        if isinstance(base_ring, FractionField_generic):
            base_ring = base_ring.ring()
        if isinstance(base_ring, (PolynomialRing_generic,
                                  MPolynomialRing_base)):
            base_ring = base_ring.base_ring()
        elif base_ring in FunctionFields():
            base_ring = base_ring.constant_base_field()
        if not (base_ring in NumberFields() or base_ring == ZZ or isinstance(base_ring, sage.rings.abc.Order)
                or (base_ring in FiniteFields())):
            raise NotImplementedError("incompatible base field, see documentation")

        #now we find the two polynomials for the resultant
        Fn = self.nth_iterate_map(n)
        fn = Fn.dehomogenize(1)
        R = fn.domain().coordinate_ring()
        S = PolynomialRing(FractionField(self.base_ring()), 'z', 2)
        phi = R.hom([S.gen(0)], S)
        psi = dom.coordinate_ring().hom([S.gen(0), 1], S)  #dehomogenize
        dfn = fn[0].derivative(R.gen())

        #polynomial to be evaluated at the periodic points
        mult_poly = phi(dfn.denominator())*S.gen(1) - phi(dfn.numerator()) #w-f'(z)

        #polynomial defining the periodic points
        x,y = dom.gens()
        if formal:
            fix_poly = self.dynatomic_polynomial(n)  #f(z)-z
        else:
            fix_poly = Fn[0]*y - Fn[1]*x #f(z) - z

        #check infinity
        inf = dom(1,0)
        inf_per = ZZ(1)
        Q = self(inf)
        while Q != inf and inf_per <= n:
            inf_per += 1
            Q = self(Q)
        #get multiplicity
        if inf_per <= n:
            e_inf = 0
            while (y**(e_inf + 1)).divides(fix_poly):
                e_inf += 1

        if type == 'cycle':
            #now we need to deal with having the correct number of factors
            #1 multiplier for each cycle. But we need to be careful about
            #the length of the cycle and the multiplicities
            good_res = 1
            if formal:
                #then we are working with the n-th dynatomic and just need
                #to take one multiplier per cycle

                #evaluate the resultant
                fix_poly = psi(fix_poly)
                res = fix_poly.resultant(mult_poly, S.gen(0))
                #take infinity into consideration
                if inf_per.divides(n):
                    res *= (S.gen(1) - self.multiplier(inf, n)[0,0])**e_inf
                res = res.univariate_polynomial()
                #adjust multiplicities
                L = res.factor()
                for p,exp in L:
                    good_res *= p**(exp/n)
            else:
                #For each d-th dynatomic for d dividing n, take
                #one multiplier per cycle; e.g., this treats a double 2
                #cycle as a single 4 cycle for n=4
                for d in n.divisors():
                    fix_poly_d = self.dynatomic_polynomial(d)
                    resd = mult_poly.resultant(psi(fix_poly_d), S.gen(0))
                    #check infinity
                    if inf_per == d:
                        e_inf_d = 0
                        while (y**(e_inf_d + 1)).divides(fix_poly_d):
                            e_inf_d += 1
                        resd *= (S.gen(1) - self.multiplier(inf, n)[0,0])**e_inf
                    resd = resd.univariate_polynomial()
                    Ld = resd.factor()
                    for pd,ed in Ld:
                        good_res *= pd**(ed/d)
            res = good_res
        else: #type is 'point'
            #evaluate the resultant
            fix_poly = psi(fix_poly)
            res = fix_poly.resultant(mult_poly, S.gen(0))
            #take infinity into consideration
            if inf_per.divides(n):
                res *= (S.gen(1) - self.multiplier(inf, n)[0,0])**e_inf
            res = res.univariate_polynomial()

        # the sigmas are the coefficients
        # needed to fix the signs and the order
        sig = res.coefficients(sparse=False)
        den = sig.pop(-1)
        sig.reverse()
        sig = [sig[i] * (-1)**(i+1) / den for i in range(len(sig))]
        return sig

    def reduced_form(self, **kwds):
        r"""
        Return reduced form of this dynamical system.

        The reduced form is the `SL(2, \ZZ)` equivalent morphism obtained
        by applying the binary form reduction algorithm from Stoll and
        Cremona [CS2003]_ to the homogeneous polynomial defining the periodic
        points (the dynatomic polynomial). The smallest period `n` with
        enough periodic points is used and without roots of too large
        multiplicity.

        This should also minimize the size  of the coefficients,
        but this is not always the case. By default the coefficient minimizing
        algorithm in [HS2018]_ is applied.

        See :meth:`sage.rings.polynomial.multi_polynomial.reduced_form` for
        the information on binary form reduction.

        Implemented by Rebecca Lauren Miller as part of GSOC 2016.
        Minimal height added by Ben Hutz July 2018.

        INPUT: keyword arguments:

        - ``prec`` -- integer (default: 300); desired precision

        - ``return_conjuagtion`` -- boolean (default: ``True``); return
          an element of `SL(2, \ZZ)`

        - ``error_limit`` -- (default: 0.000001) a real number, sets
          the error tolerance

        - ``smallest_coeffs`` -- boolean (default: ``True``); whether to find the
          model with smallest coefficients

        - ``dynatomic`` -- boolean (default: ``True``); to use formal periodic points

        - ``start_n`` -- positive integer (default: 1); first period to try to find
          appropriate binary form

        - ``emb`` -- (optional) embedding of based field into CC

        - ``algorithm`` -- (optional) which algorithm to use to find all
          minimal models. Can be one of the following:

          * ``'BM'`` -- Bruin-Molnar algorithm [BM2012]_
          * ``'HS'`` -- Hutz-Stoll algorithm [HS2018]_

        - ``check_minimal`` -- boolean (default: ``True``); whether to check
          if this map is a minimal model

        - ``smallest_coeffs`` -- boolean (default: ``True``); whether to find the
          model with smallest coefficients

        OUTPUT:

        - a projective morphism

        - a matrix

        EXAMPLES::

            sage: PS.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSystem_projective([x^3 + x*y^2, y^3])
            sage: m = matrix(QQ, 2, 2, [-201221, -1, 1, 0])
            sage: f = f.conjugate(m)
            sage: f.reduced_form(prec=50, smallest_coeffs=False)  # this needs 2 periodic
            Traceback (most recent call last):
            ...
            ValueError: accuracy of Newton's root not within tolerance(0.00006... > 1e-06), increase precision
            sage: f.reduced_form(smallest_coeffs=False)
            (
            Dynamical System of Projective Space of dimension 1 over Rational Field
              Defn: Defined on coordinates by sending (x : y) to (x^3 + x*y^2 : y^3)
            ,
            <BLANKLINE>
            [     0     -1]
            [     1 201221]
            )

        ::

            sage: PS.<x,y> = ProjectiveSpace(ZZ, 1)
            sage: f = DynamicalSystem_projective([x^2 + x*y, y^2])  # this needs 3 periodic
            sage: m = matrix(QQ, 2, 2, [-221, -1, 1, 0])
            sage: f = f.conjugate(m)
            sage: f.reduced_form(prec=200, smallest_coeffs=False)
            (
            Dynamical System of Projective Space of dimension 1 over Integer Ring
              Defn: Defined on coordinates by sending (x : y) to (-x^2 + x*y - y^2 : -y^2)
            ,
            [  0  -1]
            [  1 220]
            )

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSystem_projective([x^3, y^3])
            sage: f.reduced_form(smallest_coeffs=False)
            (
            Dynamical System of Projective Space of dimension 1 over Rational Field
              Defn: Defined on coordinates by sending (x : y) to (x^3 : y^3)
            ,
            <BLANKLINE>
            [1 0]
            [0 1]
            )

        ::

            sage: PS.<X,Y> = ProjectiveSpace(QQ,1)
            sage: f = DynamicalSystem_projective([7365*X^4 + 12564*X^3*Y + 8046*X^2*Y^2
            ....:                                   + 2292*X*Y^3 + 245*Y^4,
            ....:                                 -12329*X^4 - 21012*X^3*Y - 13446*X^2*Y^2
            ....:                                   - 3828*X*Y^3 - 409*Y^4])
            sage: f.reduced_form(prec=30, smallest_coeffs=False)
            Traceback (most recent call last):
            ...
            ValueError: accuracy of Newton's root not within tolerance(0.00008... > 1e-06), increase precision
            sage: f.reduced_form(smallest_coeffs=False)
            (
            Dynamical System of Projective Space of dimension 1 over Rational Field
              Defn: Defined on coordinates by sending (X : Y) to
                    (-7*X^4 - 12*X^3*Y - 42*X^2*Y^2 - 12*X*Y^3 - 7*Y^4
                     : -X^4 - 4*X^3*Y - 6*X^2*Y^2 - 4*X*Y^3 - Y^4),
            <BLANKLINE>
            [-1  2]
            [ 2 -5]
            )

        ::

            sage: # needs sage.rings.real_mpfr sage.symbolic
            sage: P.<x,y> = ProjectiveSpace(RR, 1)
            sage: f = DynamicalSystem_projective([x^4, RR(sqrt(2))*y^4])
            sage: m = matrix(RR, 2, 2, [1,12,0,1])
            sage: f = f.conjugate(m)
            sage: g, m = f.reduced_form(smallest_coeffs=False); m
            [  1 -12]
            [  0   1]

        ::

            sage: # needs sage.rings.real_mpfr sage.symbolic
            sage: P.<x,y> = ProjectiveSpace(CC, 1)
            sage: f = DynamicalSystem_projective([x^4, CC(sqrt(-2))*y^4])
            sage: m = matrix(CC, 2, 2, [1,12,0,1])
            sage: f = f.conjugate(m)
            sage: g, m = f.reduced_form(smallest_coeffs=False); m
            [  1 -12]
            [  0   1]

        ::

            sage: # needs sage.rings.number_field
            sage: K.<w> = QuadraticField(2)
            sage: P.<x,y> = ProjectiveSpace(K, 1)
            sage: f = DynamicalSystem_projective([x^3, w*y^3])
            sage: m = matrix(K, 2, 2, [1,12,0,1])
            sage: f = f.conjugate(m)
            sage: f.reduced_form(smallest_coeffs=False)
            (
            Dynamical System of Projective Space of dimension 1 over Number Field in w
             with defining polynomial x^2 - 2 with w = 1.414213562373095?
              Defn: Defined on coordinates by sending (x : y) to (x^3 : w*y^3) ,
            <BLANKLINE>
            [  1 -12]
            [  0   1]
            )

        ::

            sage: # needs sage.rings.number_field
            sage: R.<x> = QQ[]
            sage: K.<w> = NumberField(x^5 + x - 3,
            ....:                     embedding=(x^5 + x - 3).roots(ring=CC)[0][0])
            sage: P.<x,y> = ProjectiveSpace(K, 1)
            sage: f = DynamicalSystem_projective([12*x^3, 2334*w*y^3])
            sage: m = matrix(K, 2, 2, [-12,1,1,0])
            sage: f = f.conjugate(m)
            sage: f.reduced_form(smallest_coeffs=False)
            (
            Dynamical System of Projective Space of dimension 1 over Number Field in w
             with defining polynomial x^5 + x - 3 with w = 1.132997565885066?
              Defn: Defined on coordinates by sending (x : y) to
                    (12*x^3 : (2334*w)*y^3) ,
            <BLANKLINE>
            [  0  -1]
            [  1 -12]
            )

        ::

            sage: P.<x,y> = QQ[]
            sage: f = DynamicalSystem([-4*y^2, 9*x^2 - 12*x*y])
            sage: f.reduced_form()
            (
            Dynamical System of Projective Space of dimension 1 over Rational Field
              Defn: Defined on coordinates by sending (x : y) to
                    (2*x^2 - 2*y^2 : -x^2 - 2*y^2)
            ,
            <BLANKLINE>
            [ 2 -2]
            [ 3  0]
            )

        ::

            sage: P.<x,y> = QQ[]
            sage: f = DynamicalSystem([-2*x^3 - 9*x^2*y - 12*x*y^2 - 6*y^3 , y^3])
            sage: f.reduced_form()
            (
            Dynamical System of Projective Space of dimension 1 over Rational Field
              Defn: Defined on coordinates by sending (x : y) to
                    (x^3 + 3*x^2*y : 3*x*y^2 + y^3)
            ,
            <BLANKLINE>
            [-1 -2]
            [ 1  1]
            )

        ::

            sage: P.<x,y> = QQ[]
            sage: f = DynamicalSystem([4*x^2 - 7*y^2, 4*y^2])
            sage: f.reduced_form(start_n=2, dynatomic=False) #long time
            (
            Dynamical System of Projective Space of dimension 1 over Rational Field
              Defn: Defined on coordinates by sending (x : y) to
                    (x^2 - x*y - y^2 : y^2)
            ,
            <BLANKLINE>
            [ 2 -1]
            [ 0  2]
            )

        ::

            sage: P.<x,y> = QQ[]
            sage: f = DynamicalSystem([4*x^2 + y^2, 4*y^2])
            sage: f.reduced_form()  # long time
            (
            Dynamical System of Projective Space of dimension 1 over Rational Field
              Defn: Defined on coordinates by sending (x : y) to
                    (x^2 - x*y + y^2 : y^2)
            ,
            <BLANKLINE>
            [ 2 -1]
            [ 0  2]
            )
        """
        if self.domain().ambient_space().dimension_relative() != 1:
            return NotImplementedError('only implemented for dimension 1')
        return_conjugation = kwds.get('return_conjugation', True)
        emb = kwds.get('emb', None)
        prec = kwds.get('prec', 300)
        start_n = kwds.get('start_n', 1)
        algorithm = kwds.get('algorithm', None)
        dynatomic = algorithm = kwds.get('dynatomic', True)
        smallest_coeffs = kwds.get('smallest_coeffs', True)
        if smallest_coeffs:
            if self.base_ring() not in [ZZ, QQ]:
                raise NotImplementedError("smallest coeff only over ZZ or QQ")
            check_min = kwds.get('check_minimal', True)
            from sage.dynamics.arithmetic_dynamics.endPN_minimal_model import smallest_dynamical
            sm_f, m = smallest_dynamical(self, dynatomic=dynatomic, start_n=start_n,
                 prec=prec, emb=emb, algorithm=algorithm, check_minimal=check_min)
        else:
            #reduce via covariant
            PS = self.domain()
            CR = PS.coordinate_ring()
            x,y = CR.gens()
            n = start_n # sometimes you get a problem later with 0,infty as roots
            pts_poly = self.dynatomic_polynomial(n)
            d = ZZ(pts_poly.degree())
            try:
                max_mult = max([ex for p,ex in pts_poly.factor()])
            except NotImplementedError: #not factorization in numerical rings
                CF = ComplexField(prec=prec)
                if pts_poly.base_ring() != CF:
                    if emb is None:
                        pts_poly_CF = pts_poly.change_ring(CF)
                    else:
                        pts_poly_CF = pts_poly.change_ring(emb)
                pp_d = pts_poly.degree()
                pts_poly_CF = pts_poly_CF.subs({pts_poly_CF.parent().gen(1):1}).univariate_polynomial()
                max_mult = max([pp_d - pts_poly_CF.degree()] + [ex for p,ex in pts_poly_CF.roots()])
            while ((d < 3) or (max_mult >= d/2) and (n < 5)):
                n = n+1
                if dynatomic:
                    pts_poly = self.dynatomic_polynomial(n)
                else:
                    gn = self.nth_iterate_map(n)
                    pts_poly = y*gn[0] - x*gn[1]
                d = ZZ(pts_poly.degree())
                try:
                    max_mult = max([ex for p,ex in pts_poly.factor()])
                except NotImplementedError: #not factorization in numerical rings
                    CF = ComplexField(prec=prec)
                    if pts_poly.base_ring() != CF:
                        if emb is None:
                            pts_poly_CF = pts_poly.change_ring(CF)
                        else:
                            pts_poly_CF = pts_poly.change_ring(emb)
                    pp_d = pts_poly.degree()
                    pts_poly_CF = pts_poly_CF.subs({pts_poly_CF.parent().gen(1):1}).univariate_polynomial()
                    max_mult = max([pp_d - pts_poly_CF.degree()] + [ex for p,ex in pts_poly_CF.roots()])
            assert (n <= 4), "n > 4, failed to find usable poly"
            G,m = pts_poly.reduced_form(prec=prec, emb=emb, smallest_coeffs=False)
            sm_f = self.conjugate(m)

        if return_conjugation:
            return (sm_f, m)
        return sm_f

    def _is_preperiodic(self, P, err=0.1, return_period=False):
        r"""
        Determine if the point is preperiodic with respect to this
        dynamical system.

        ALGORITHM:

        We know that a point is preperiodic if and only if it has
        canonical height zero. However, we can only compute the canonical
        height up to numerical precision. This function first computes
        the canonical height of the point to the given error bound. If
        it is larger than that error bound, then it must not be preperiodic.
        If it is less than the error bound, then we expect preperiodic. In
        this case we begin computing the orbit stopping if either we
        determine the orbit is finite, or the height of the point is large
        enough that it must be wandering. We can determine the height
        cutoff by computing the height difference constant, i.e., the bound
        between the height and the canonical height of a point (which
        depends only on the map and not the point itself). If the height
        of the point is larger than the difference bound, then the canonical
        height cannot be zero so the point cannot be preperiodic.

        INPUT:

        - ``P`` -- a point of this dynamical system's codomain

        kwds:

        - ``err`` -- (default: 0.1) a positive real number;
          sets the error_bound used in the canonical height computation
          and ``return_period`` a boolean which

        - ``return_period`` -- boolean (default: ``False``); controls if
          the period is returned if the point is preperiodic

        OUTPUT:

        - boolean -- ``True`` if preperiodic

        - if ``return_period`` is ``True``, then ``(0,0)`` if wandering,
          and ``(m,n)`` if preperiod ``m`` and period ``n``

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: f = DynamicalSystem_projective([x^3 - 3*x*y^2, y^3], domain=P)
            sage: Q = P(-1, 1)
            sage: f._is_preperiodic(Q)                                                  # needs sage.rings.function_field
            True

        Check that :issue:`23814` is fixed (works even if domain is not specified)::

            sage: # needs sage.rings.number_field
            sage: R.<X> = PolynomialRing(QQ)
            sage: K.<a> = NumberField(X^2 + X - 1)
            sage: P.<x,y> = ProjectiveSpace(K,1)
            sage: f = DynamicalSystem_projective([x^2 - 2*y^2, y^2])
            sage: Q = P.point([a,1])
            sage: Q.is_preperiodic(f)                                                   # needs sage.rings.function_field
            True

        ::

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: X = P.subscheme(z)
            sage: f = DynamicalSystem([x^2 - y^2, y^2, z^2], domain=X)
            sage: p = X((-1, 1, 0))
            sage: f._is_preperiodic(p, return_period=True)                              # needs sage.rings.function_field
            (0, 2)

        ::

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: X = P.subscheme(x)
            sage: f = DynamicalSystem([x^2 - y^2, y^2, z^2], domain=X)
            sage: p = X((0, 1, 0))
            sage: f._is_preperiodic(p, return_period=True)                              # needs sage.rings.function_field
            Traceback (most recent call last):
            ...
            ValueError: orbit of point leaves domain

        ::

            sage: # needs sage.rings.number_field
            sage: R.<t> = QQ[]
            sage: K.<a> = NumberField(t^2 - t - 1)
            sage: P.<x,y,z> = ProjectiveSpace(K, 2)
            sage: X = P.subscheme(z)
            sage: f = DynamicalSystem([x^2 - y^2, y^2, z^2], domain=X)
            sage: p = X((-a + 1, 1, 0))
            sage: f._is_preperiodic(p)                                                  # needs sage.rings.function_field
            True
        """
        codomain = self.codomain()
        if not isinstance(codomain, ProjectiveSpace_ring):
            # in order to calculate the canonical height, we need
            # this map to be a morphism of projective space
            ambient_space = codomain.ambient_space()
            f = DynamicalSystem(self.defining_polynomials(), domain=ambient_space)
            if not f.is_morphism():
                raise ValueError('must be a morphism of projective space')
        else:
            f = self
            if not f.is_morphism():
                raise TypeError("must be a morphism")
        if not P.codomain() == self.domain():
            raise TypeError("point must be in domain of map")

        # we calculate the canonical height without considering
        # if the domain is a subscheme
        h = f.canonical_height(P, error_bound=err)
        # we know canonical height 0 if and only if preperiodic
        # however precision issues can occur so we can only tell *not* preperiodic
        # if the value is larger than the error
        if h <= err:
            # if the canonical height is less than the
            # error, then we suspect preperiodic so check
            # either we can find the cycle or the height is
            # larger than the difference between the canonical height
            # and the height, so the canonical height cannot be 0
            B = f.height_difference_bound()
            orbit = [P]
            n = 1 # to compute period
            try:
                Q = self(P)
            except TypeError:
                raise ValueError('orbit of point leaves domain')
            H = Q.global_height()
            while Q not in orbit and H <= B:
                orbit.append(Q)
                try:
                    Q = self(Q)
                except TypeError:
                    raise ValueError('orbit of point leaves domain')
                H = Q.global_height()
                n += 1
            if H <= B: #it must have been in the cycle
                if return_period:
                    m = orbit.index(Q)
                    return (m, n - m)
                else:
                    return True
        if return_period:
            return (0,0)
        else:
            return False

    def postcritical_set(self, check=True):
        r"""
        Return the postcritical set of this dynamical system.

        Raises an error if this dynamical system is not postcritically finite.

        The postcritical set is union of points which are in the forward orbits
        of the critical points. In other words, the set of points `Q` such that
        `f^n(P) = Q` for some positive integer `n` and critical point `P`, where
        `f` is this map.

        Note that the orbit of all critical points is found, even if the
        critical points are defined in an extension of the base ring of
        this dynamical system. We extend to the field defined by
        ``f.field_of_definition_critical()``, where ``f`` is this map.

        INPUT:

        - ``check`` -- boolean (default: ``True``); whether to check
          if this dynamical system is postcritically finite or not

        OUTPUT: the set of postcritical points

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSystem([x^3 - 3/2* x*y^2, y^3])
            sage: f.postcritical_set()                                                  # needs sage.rings.number_field
            [(1/2*a : 1), (-1/2*a : 1), (1 : 0)]

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSystem([3*x^3 - 9/2* x^2*y+y^3, y^3])
            sage: f.postcritical_set(check=False)                                       # needs sage.rings.number_field
            [(1 : 1), (-1/2 : 1), (1 : 0)]

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSystem([-4*y^2, 9*x^2 - 12*x*y])
            sage: f.postcritical_set()                                                  # needs sage.rings.number_field
            [(1 : 1), (4/3 : 1), (1 : 0), (0 : 1)]

        ::

            sage: # needs sage.rings.number_field
            sage: K.<v> = QuadraticField(2)
            sage: P.<x,y> = ProjectiveSpace(K,1)
            sage: f = DynamicalSystem([x^2 + (-2)*y^2, y^2])
            sage: m = matrix(K, 2, 2, [v, 1, 0, 1])
            sage: g = f.conjugate(m)
            sage: g.postcritical_set()
            [(-3/2*a : 1), (1/2*a : 1), (1 : 0)]

        ::

            sage: # needs sage.rings.finite_rings
            sage: F.<z> = FiniteField(9)
            sage: P.<x,y> = ProjectiveSpace(F, 1)
            sage: f = DynamicalSystem([x^2 + (-2)*y^2, y^2])
            sage: m = matrix(F, 2, 2, [z, 1, 0, 1])
            sage: g = f.conjugate(m)
            sage: g.postcritical_set()
            [(1 : 0), (0 : 1), (a + 2 : 1)]
        """
        if not isinstance(self.domain(), ProjectiveSpace_ring):
            raise ValueError('must be a dynamical system on projective space')
        if self.domain().dimension_relative() != 1:
            raise ValueError('must be defined on projective space of dimension 1')
        if check:
            if not self.is_postcritically_finite():
                raise ValueError('map must be postcritically finite')
        new_base_ring = self.field_of_definition_critical(return_embedding=True)[1]
        f = self.change_ring(new_base_ring)
        critical_points = f.critical_points()
        post_critical_list = []
        for point in critical_points:
            next_point = f(point)
            while next_point not in post_critical_list:
                post_critical_list.append(next_point)
                next_point = f(next_point)
        return post_critical_list

    def is_chebyshev(self):
        r"""
        Check if ``self`` is a Chebyshev polynomial.

        OUTPUT: ``True`` if ``self`` is Chebyshev, ``False`` otherwise

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: F = DynamicalSystem_projective([x^4, y^4])
            sage: F.is_chebyshev()
            False

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: F = DynamicalSystem_projective([x^2 + y^2, y^2])
            sage: F.is_chebyshev()
            False

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: F = DynamicalSystem_projective([2*x^2 - y^2, y^2])
            sage: F.is_chebyshev()
            True

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: F = DynamicalSystem_projective([x^3, 4*y^3 - 3*x^2*y])
            sage: F.is_chebyshev()
            True

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: F = DynamicalSystem_projective([2*x^2 - y^2, y^2])
            sage: L.<i> = CyclotomicField(4)
            sage: M = Matrix([[0,i],[-i,0]])
            sage: F.conjugate(M)
            Dynamical System of Projective Space of dimension 1 over
             Cyclotomic Field of order 4 and degree 2
              Defn: Defined on coordinates by sending (x : y) to
                    ((-i)*x^2 : (-i)*x^2 + (2*i)*y^2)
            sage: F.is_chebyshev()
            True

        REFERENCES:

        - [Mil2006]_
        """
        # We need `F` to be defined over a number field for
        # the function `is_postcrtically_finite` to work
        if self.base_ring() not in NumberFields():
            raise NotImplementedError("Base ring must be a number field")

        if self.domain().dimension() != 1:
            return False

        # All Chebyshev polynomials are postcritically finite
        if not self.is_postcritically_finite():
            return False

        # Get field of definition for critical points and change base field
        critical_field, phi = self.field_of_definition_critical(return_embedding=True)
        F_crit = self.change_ring(phi)

        # Get the critical points and post-critical set
        crit_set = set(F_crit.critical_points())
        post_crit_set = set()
        images_needed = copy(crit_set)
        while len(images_needed) != 0:
            Q = images_needed.pop()
            Q2 = F_crit(Q)
            if Q2 not in post_crit_set:
                post_crit_set.add(Q2)
                images_needed.add(Q2)

        crit_fixed_pts = set()
        for crit in crit_set:
            if F_crit.nth_iterate(crit, 1) == crit:
                crit_fixed_pts.add(crit)

        # All Chebyshev maps have 3 post-critical values
        if (len(post_crit_set) != 3) or (len(crit_fixed_pts) != 1):
            return False

        f = F_crit.dehomogenize(1)[0]
        x = f.parent().gen()
        ram_points = {}
        for crit in crit_set:
            g = f
            new_crit = crit

            # Check if critical point is infinity
            if crit[1] == 0:
                g = g.subs(x=1/x)
                new_crit = F_crit.domain()([0, 1])

            # Check if output is infinity
            if F_crit.nth_iterate(crit, 1)[1] == 0:
                g = 1/g

            new_crit = new_crit.dehomogenize(1)[0]
            e = 1
            g = g.derivative(x)

            while g(new_crit) == 0:
                e += 1
                g = g.derivative(x)

            ram_points[crit] = e

        r = {}

        # Set r value to 0 to represent infinity
        r[crit_fixed_pts.pop()] = 0

        # Get non-fixed tail points in the post-critical portrait
        crit_tails = crit_set.difference(post_crit_set)
        for crit in crit_tails:
            # Each critical tail point has r value 1
            r[crit] = 1
            point = crit

            # Assign r values to every point in the orbit of crit.
            # If we find a point in the orbit which already has an r value assigned,
            # check that we get a consistent r value.
            while F_crit.nth_iterate(point, 1) not in r.keys():
                if point not in ram_points.keys():
                    r[F_crit.nth_iterate(point,1)] = r[point]
                else:
                    r[F_crit.nth_iterate(point,1)] = r[point] * ram_points[point]

                point = F_crit.nth_iterate(point,1)

            # Once we get here, the image of point has an assigned r value
            # We check that this value is consistent
            if point not in ram_points.keys():
                if r[F_crit.nth_iterate(point,1)] != r[point]:
                    return False
            elif r[F_crit.nth_iterate(point,1)] != r[point] * ram_points[point]:
                return False

        # The non-one r values must be one of the following in order for F to be Chebyshev
        r_vals = sorted([val for val in r.values() if val != 1])
        return r_vals == [0, 2, 2]

    def is_Lattes(self):
        r"""
        Check if ``self`` is a Lattes map.

        OUTPUT: ``True`` if ``self`` is Lattes, ``False`` otherwise

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: F = DynamicalSystem_projective([x^3, y^3])
            sage: F.is_Lattes()
            False

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: F = DynamicalSystem_projective([x^2 - 2*y^2, y^2])
            sage: F.is_Lattes()
            False

        ::

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: F = DynamicalSystem_projective([x^2 + y^2 + z^2, y^2, z^2])
            sage: F.is_Lattes()
            False

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: F = DynamicalSystem_projective([(x + y)*(x - y)^3, y*(2*x + y)^3])
            sage: F.is_Lattes()
            True

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: F = DynamicalSystem_projective([(x + y)^4, 16*x*y*(x - y)^2])
            sage: F.is_Lattes()
            True

        ::

            sage: f = P.Lattes_map(EllipticCurve([0, 0, 0, 0, 2]),2)
            sage: f.is_Lattes()
            True

        ::

            sage: f = P.Lattes_map(EllipticCurve([0, 0, 0, 0, 2]), 2)
            sage: L.<i> = CyclotomicField(4)
            sage: M = Matrix([[i, 0], [0, -i]])
            sage: f.conjugate(M)
            Dynamical System of Projective Space of dimension 1 over
             Cyclotomic Field of order 4 and degree 2
              Defn: Defined on coordinates by sending (x : y) to
                    ((-1/4*i)*x^4 + (-4*i)*x*y^3 : (-i)*x^3*y + (2*i)*y^4)
            sage: f.is_Lattes()
            True

        REFERENCES:

        - [Mil2006]_
        """
        # We need `f` to be defined over a number field for
        # the function `is_postcrtically_finite` to work
        if self.base_ring() is not QQbar:
            if self.base_ring() not in NumberFields():
                    raise NotImplementedError("Base ring must be a number field")

        if self.domain().dimension() != 1:
            return False

        # All Lattes maps are postcritically finite
        if not self.is_postcritically_finite():
            return False

        # Get field of definition for critical points and change basefield
        critical_field, phi = self.field_of_definition_critical(return_embedding=True)
        F_crit = self.change_ring(phi)

        # Get the critical points and post-critical set
        crit = F_crit.critical_points()
        post_crit = set()
        images_needed = copy(crit)
        while len(images_needed) != 0:
            Q = images_needed.pop()
            Q2 = F_crit(Q)
            if Q2 not in post_crit:
                post_crit.add(Q2)
                images_needed.append(Q2)

        (crit_set, post_crit_set) = crit, list(post_crit)

        # All Lattes maps have 3 or 4 post critical values
        if len(post_crit_set) not in [3, 4]:
            return False

        f = F_crit.dehomogenize(1)[0]
        x = f.parent().gen()
        ram_points = {}
        for crit in crit_set:
            g = f
            new_crit = crit

            # Check if critical point is infinity
            if crit[1] == 0:
                g = g.subs(x=1/x)
                new_crit = F_crit.domain()([0, 1])

            # Check if output is infinity
            if F_crit.nth_iterate(crit, 1)[1] == 0:
                g = 1/g

            new_crit = new_crit.dehomogenize(1)[0]

            e = 1
            g = g.derivative(x)
            while g(new_crit) == 0:
                e += 1
                g = g.derivative(x)
            ram_points[crit] = e

        r = {}

        # Get tail points in the post-critical portrait
        crit_tails = set(crit_set).difference(set(post_crit_set))
        for crit in crit_tails:
            # Each critical tail point has r value 1
            r[crit] = 1
            point = crit

            # Assign r values to every point in the orbit of crit
            # If we find a point in the orbit which already has an r value assigned,
            # check that we get a consistent r value
            while F_crit.nth_iterate(point, 1) not in r.keys():
                if point not in ram_points.keys():
                    r[F_crit.nth_iterate(point, 1)] = r[point]
                else:
                    r[F_crit.nth_iterate(point, 1)] = r[point] * ram_points[point]

                point = F_crit.nth_iterate(point,1)

            # Once we get here the image of point has an assigned r value
            # We check that this value is consistent.
            if point not in ram_points.keys():
                if r[F_crit.nth_iterate(point, 1)] != r[point]:
                    return False
            else:
                if r[F_crit.nth_iterate(point, 1)] != r[point] * ram_points[point]:
                    return False

        # The non-one r values must be one of the following in order for F to be Lattes
        r_lattes_cases = [[2, 2, 2, 2], [3, 3, 3], [2, 4, 4], [2, 3, 6]]
        r_vals = sorted([val for val in r.values() if val != 1])
        return r_vals in r_lattes_cases

    def Lattes_to_curve(self, return_conjugation=False, check_lattes=False):
        r"""
        Finds a Short Weierstrass Model Elliptic curve of self
        self assumed to be Lattes map and not in charateristic 2 or 3

        INPUT:

        ``return_conjugation`` -- (default: ``False``) if ``True``, then
        return the conjugation that moves self to a map that comes from a
        Short Weierstrass Model Elliptic curve
        ``check_lattes``.-.(default:.``False``) if ``True``, then  will ValueError if not Lattes

        OUTPUT: a Short Weierstrass Model Elliptic curve which is isogenous to
        the Elliptic curve of 'self',
        If ``return_conjugation`` is ``True``
        then also returns conjugation of 'self' to short form as a matrix

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f = P.Lattes_map(EllipticCurve([0, 0, 0, 10, 2]), 2)
            sage: f.Lattes_to_curve()
            Elliptic Curve defined by y^2 = x^3 + 10*x + 2 over Rational Field

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: M = matrix(QQ,2,2,[[1,2],[-1,2]])
            sage: f = P.Lattes_map(EllipticCurve([1, 1, 1, 1, 2]), 2)
            sage: f = f.conjugate(M)
            sage: f.Lattes_to_curve(return_conjugation = True)
            (
            [  -7/36*a^2 + 7/12*a + 7/3 -17/18*a^2 + 17/6*a + 34/3]
            [    -1/8*a^2 + 1/4*a + 3/2        1/4*a^2 - 1/2*a - 3],
            Elliptic Curve defined by y^2 = x^3 + (-94/27*a^2+94/9*a+376/9)*x +
            12232/243 over Number Field in a with defining polynomial y^3 - 18*y - 30
            )

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: f = P.Lattes_map(EllipticCurve([1, 1, 1, 2, 2]), 2)
            sage: L.<i> = CyclotomicField(4)
            sage: M = Matrix([[1+i,2*i], [0, -i]])
            sage: f = f.conjugate(M)
            sage: f.Lattes_to_curve(return_conjugation = True)
            (
            [              1 19/24*a + 19/24]
            [              0               1],
            Elliptic Curve defined by y^2 = x^3 + 95/96*a*x + (-1169/3456*a+1169/3456)
            over Number Field in a with defining polynomial y^2 + 1
            )

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: M = matrix(QQ,2,2,[[1,3],[2,1]])
            sage: E = EllipticCurve([1, 1, 1, 2, 3])
            sage: f = P.Lattes_map(E, 2)
            sage: f = f.conjugate(M)
            sage: f.Lattes_to_curve(return_conjugation = True)
            (
            [11/1602*a^2 41/3204*a^2]
            [     -2/5*a      -1/5*a],
            Elliptic Curve defined by y^2 = x^3 + 2375/3421872*a^2*x + (-254125/61593696)
            over Number Field in a with defining polynomial y^3 - 267
            )

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ , 1)
            sage: M = matrix(QQ,2,2,[[1 , 3],[2 , 1]])
            sage: E = EllipticCurve([1, 1, 1, 2, 3])
            sage: f = P.Lattes_map(E , 2)
            sage: f = f.conjugate(M)
            sage: m,H = f.Lattes_to_curve(true)
            sage: J.<x,y> = ProjectiveSpace(H.base_ring(), 1)
            sage: K = J.Lattes_map(H,2)
            sage: K = K.conjugate(m)
            sage: K.scale_by(f[0].lc()/K[0].lc())
            sage: K == f.change_ring(K.base_ring())
            True

        ::

            sage: P.<x,y> = ProjectiveSpace(RR, 1)
            sage: F = DynamicalSystem_projective([x^4, y^4])
            sage: F.Lattes_to_curve(check_lattes=True)
            Traceback (most recent call last):
            ...
            NotImplementedError: Base ring must be a number field

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: F = DynamicalSystem_projective([x^4, y^4])
            sage: F.Lattes_to_curve(check_lattes=True)
            Traceback (most recent call last):
            ...
            ValueError: Map is not Lattes

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: F = DynamicalSystem_projective([x^4, y^4])
            sage: F.Lattes_to_curve()
            Traceback (most recent call last):
            ...
            ValueError: No Solutions found. Check if map is Lattes

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: F = DynamicalSystem_projective([x^3, y^3])
            sage: F.Lattes_to_curve(check_lattes=True)
            Traceback (most recent call last):
            ...
            NotImplementedError: Map is not Lattes or is Complex Lattes

        ::

            sage: K.<x>=QuadraticField(2)
            sage: P.<a,y>=ProjectiveSpace(K, 1)
            sage: E=EllipticCurve([1, x])
            sage: f=P.Lattes_map(E, 2)
            sage: f.Lattes_to_curve()
            Elliptic Curve defined by y^2 = x^3 + x + a
            over Number Field in a with defining polynomial y^2 - 2

        ::

            sage: P.<x,y>=ProjectiveSpace(QQbar, 1)
            sage: E=EllipticCurve([1, 2])
            sage: f=P.Lattes_map(E, 2)
            sage: f.Lattes_to_curve(check_lattes=true)
            Elliptic Curve defined by y^2 = x^3 + x + 2 over Rational Field

        """
        if self.base_ring() is not QQbar:
            if self.base_ring() not in NumberFields():
                    raise NotImplementedError("Base ring must be a number field")
    #The Complex case is hard to implement and needs to be done later
        if sqrt(self.degree()) != int(sqrt(self.degree())):
            raise NotImplementedError("Map is not Lattes or is Complex Lattes")
        if check_lattes:
            V = self.is_Lattes()
            if not V:
                raise ValueError("Map is not Lattes")
        n = int(sqrt(self.degree()))
    #Creating a Symbolic Lattes map f_sym from a short Elliptic curve
        R = PolynomialRing(self.base_ring(), 6, "avar, bvar, uvar, vvar, wvar, tvar")
        a, b, u, v, w, t = R.gens()
        P = ProjectiveSpace(R, 1, self.domain().gens())
        E_sym = EllipticCurve([a, b])
        f_sym = P.Lattes_map(E_sym, n)
    # Conjugating f_sym map to have the right form so we can solve for the conjugating matrix later
        m = matrix(R, 2, [u, v, t, w])
        f_sym = f_sym.conjugate(m)
        f_sym.scale_by(u*w - v*t)
        F_sym = f_sym.dehomogenize(1)
    #extracting the base variables to do term by term matching
        self.scale_by(1/self[0].lc())
        F = self.dehomogenize(1)
    #Creating a set of equations, eq, from term by term matching
        eq = [u*w - v*t-1]
        for j in range(2):
            if j == 0:
                g = F[0].numerator()
                g_sym = F_sym[0].numerator()
            else:
                g = F[0].denominator()
                g_sym = F_sym[0].denominator()
            eq += (g - g_sym).coefficients()
    #Solving the equations
        phi = QQbar.coerce_map_from(R.base_ring())
        if phi is None:
            phi = R.base_ring().embeddings(QQbar)[0]
        eq = [poly.numerator().change_ring(phi) for poly in eq]
        I = eq[0].parent().ideal(eq)
        pts = I.variety()
        if len(pts) == 0:
            raise ValueError("No Solutions found. Check if map is Lattes")
        a = pts[0]['avar']
        b = pts[0]['bvar']
        u = pts[0]['uvar']
        v = pts[0]['vvar']
        t = pts[0]['tvar']
        w = pts[0]['wvar']
        K, [a, b, u, v, t, w], phi = number_field_elements_from_algebraics([a, b, u, v, t, w])
    #creating our end products
        E = EllipticCurve([a, b])
        if return_conjugation:
            M = matrix(K, 2, 2, [u, v, t, w])
            return (M, E)
        return E


class DynamicalSystem_projective_field(DynamicalSystem_projective,
                                       SchemeMorphism_polynomial_projective_space_field):

    def lift_to_rational_periodic(self, points_modp, B=None):
        r"""
        Given a list of points in projective space over `\GF{p}`,
        determine if they lift to `\QQ`-rational periodic points.

        The map must be an endomorphism of projective space defined
        over `\QQ`.

        ALGORITHM:

        Use Hensel lifting to find a `p`-adic approximation for that
        rational point. The accuracy needed is determined by the height
        bound ``B``. Then apply the LLL algorithm to determine if the
        lift corresponds to a rational point.

        If the point is a point of high multiplicity (multiplier 1), the
        procedure can be very slow.

        INPUT:

        - ``points_modp`` -- list or tuple of pairs containing a point
          in projective space over `\GF{p}` and the possible period

        - ``B`` -- (optional) a positive integer; the height bound for
          a rational preperiodic point

        OUTPUT: list of projective points

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: f = DynamicalSystem_projective([x^2 - y^2, y^2])
            sage: f.lift_to_rational_periodic([[P(0,1).change_ring(GF(7)), 4]])         # needs sage.symbolic
            [[(0 : 1), 2]]

        There may be multiple points in the lift. ::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: f = DynamicalSystem_projective([-5*x^2 + 4*y^2, 4*x*y])
            sage: f.lift_to_rational_periodic([[P(1,0).change_ring(GF(3)), 1]])  # long time
            [[(1 : 0), 1], [(2/3 : 1), 1], [(-2/3 : 1), 1]]

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: f = DynamicalSystem_projective([16*x^2 - 29*y^2, 16*y^2])
            sage: f.lift_to_rational_periodic([[P(3,1).change_ring(GF(13)), 3]])        # needs sage.symbolic
            [[(-1/4 : 1), 3]]

        ::

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: f = DynamicalSystem_projective([76*x^2 - 180*x*y + 45*y^2
            ....:                                   + 14*x*z + 45*y*z - 90*z^2,
            ....:                                 67*x^2 - 180*x*y - 157*x*z + 90*y*z,
            ....:                                 -90*z^2])
            sage: f.lift_to_rational_periodic([[P(14,19,1).change_ring(GF(23)), 9]])    # long time
            [[(-9 : -4 : 1), 9]]
        """
        if not points_modp:
            return []

        if B is None:
            from sage.symbolic.constants import e
            B = e ** self.height_difference_bound()

        p = points_modp[0][0].codomain().base_ring().characteristic()
        if p == 0:
            raise TypeError("must be positive characteristic")
        PS = self.domain()
        N = PS.dimension_relative()
        R = RealField()
        #compute the maximum p-adic precision needed to conclusively determine
        #if the rational point exists
        L = R((R(2 ** (N/2 + 1) * sqrt(N+1) * B**2).log()) / R(p).log() + 1).trunc()

        points = []
        for i in range(len(points_modp)):
            #[point mod p, period, current p-adic precision]
            points.append([points_modp[i][0].change_ring(QQ, check=False), points_modp[i][1], 1])
        good_points = []
        #shifts is used in non-Hensel lifting
        shifts = None
        #While there are still points to consider try to lift to next precision
        while points:
            q = points.pop()
            qindex = N
            #Find the last nonzero coordinate to use for normalizations
            while q[0][qindex] % p == 0:
                qindex -= 1
            T = q[0]
            n = q[1]
            k = q[2]
            T.scale_by(1 / T[qindex]) #normalize
            bad = 0
            #stop where we reach the needed precision or the point is bad
            while k < L and bad == 0:
                l = self._multipliermod(T, n, p, 2*k)
                l -= l.parent().one() #f^n(x) - x
                lp = l.change_ring(Zmod(p**k))
                ldet = lp.determinant()
                # if the matrix is invertible then we can Hensel lift
                if ldet % p != 0:
                    RQ = ZZ.quo(p**(2*k))
                    T.clear_denominators()
                    newT = T.change_ring(RQ, check=False)
                    fp = self.change_ring(RQ, check=False)
                    S = fp.nth_iterate(newT, n, normalize=False).change_ring(QQ, check=False)
                    T.scale_by(1 / T[qindex])
                    S.scale_by(1 / S[qindex])
                    newS = list(S)
                    for i in range(N + 1):
                        newS[i] = S[i] - T[i]
                        if newS[i] % (p**k) != 0 and i != N:
                            bad = 1
                            break
                    if bad == 1:
                        break
                    S = PS.point(newS, False)
                    S.scale_by(-1 / p**k)
                    vecs = [Zmod(p**k)(S._coords[iS]) for iS in range(N + 1)]
                    vecs.pop(qindex)
                    newvecs = list((lp.inverse()) * vector(vecs)) #l.inverse should be mod p^k!!
                    newS = []
                    [newS.append(QQ(newvecs[i])) for i in range(qindex)]
                    newS.append(0)
                    [newS.append(QQ(newvecs[i])) for i in range(qindex, N)]
                    for i in range(N + 1):
                        newS[i] = newS[i] % (p**k)
                    S = PS.point(newS, False) #don't check for [0,...,0]
                    newT = list(T)
                    for i in range(N + 1):
                        newT[i] += S[i] * (p**k)
                    T = PS.point(newT, False)
                    T.normalize_coordinates()
                    #Hensel gives us 2k for the newprecision
                    k = min(2*k, L)
                else:
                    #we are unable to Hensel Lift so must try all possible lifts
                    #to the next precision (k+1)
                    first = 0
                    newq = []
                    RQ = Zmod(p**(k+1))
                    fp = self.change_ring(RQ, check=False)
                    if shifts is None:
                        shifts = xmrange([p for i in range(N)])
                    for shift in shifts:
                        newT = [RQ(t) for t in T]  #T.change_ring(RQ, check = False)
                        shiftindex = 0
                        for i in range(N + 1):
                            if i != qindex:
                                newT[i] = newT[i] + shift[shiftindex] * p**k
                                shiftindex += 1
                        newT = fp.domain().point(newT, check=False)
                        TT = fp.nth_iterate(newT, n, normalize=False)
                        if TT == newT:
                            if first == 0:
                                newq.append(newT.change_ring(QQ, check=False))
                                newq.append(n)
                                newq.append(k + 1)
                                first = 1
                            else:
                                points.append([newT.change_ring(QQ, check=False), n, k+1])
                    if not newq:
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
                    M[i+1, i] = p**L
                M[N+1, N] = p**L
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
                    while not done and k <= n:
                        newP = self(newP)
                        if newP == P:
                            if [P, k] not in good_points:
                                good_points.append([newP, k])
                            done = True
                        k += 1

        return good_points

    def all_periodic_points(self, **kwds):
        r"""
        Determine the set of rational periodic points
        for this dynamical system.

        The map must be defined over `\QQ` and be an endomorphism of
        projective space. If the map is a polynomial endomorphism of
        `\mathbb{P}^1`, i.e. has a totally ramified fixed point, then
        the base ring can be an absolute number field.
        This is done by passing to the Weil restriction.

        The default parameter values are typically good choices for
        `\mathbb{P}^1`. If you are having trouble getting a particular
        map to finish, try first computing the possible periods, then
        try various different ``lifting_prime`` values.

        ALGORITHM:

        Modulo each prime of good reduction `p` determine the set of
        periodic points modulo `p`. For each cycle modulo `p` compute
        the set of possible periods (`mrp^e`). Take the intersection
        of the list of possible periods modulo several primes of good
        reduction to get a possible list of minimal periods of rational
        periodic points. Take each point modulo `p` associated to each
        of these possible periods and try to lift it to a rational point
        with a combination of `p`-adic approximation and the LLL basis
        reduction algorithm.

        See [Hutz2015]_.

        INPUT: keyword arguments:

        - ``R`` -- (default: domain of dynamical system) the base ring
          over which the periodic points of the dynamical system are found

        - ``prime_bound`` -- (default: ``[1,20]``) a pair (list or tuple)
          of positive integers that represent the limits of primes to use
          in the reduction step or an integer that represents the upper bound

        - ``lifting_prime`` -- (default: 23) a prime integer; argument that
          specifies modulo which prime to try and perform the lifting

        - ``period_degree_bounds`` -- (default: ``[4,4]``) a pair of positive integers
          (max period, max degree) for which the dynatomic polynomial should be solved for

        - ``algorithm`` -- (optional) specifies which algorithm to use;
          current options are `dynatomic` and `lifting`; defaults to solving the
          dynatomic for low periods and degrees and lifts for everything else

        - ``periods`` -- (optional) a list of positive integers that is
          the list of possible periods

        - ``bad_primes`` -- (optional) a list or tuple of integer primes;
          the primes of bad reduction

        - ``ncpus`` -- (default: all cpus) number of cpus to use in parallel

        OUTPUT: list of rational points in projective space

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: f = DynamicalSystem_projective([x^2 - 3/4*y^2, y^2])
            sage: sorted(f.all_periodic_points(prime_bound=20, lifting_prime=7))  # long time
            [(-1/2 : 1), (1 : 0), (3/2 : 1)]

        ::

            sage: P.<x,y,z> = ProjectiveSpace(QQ,2)
            sage: f = DynamicalSystem_projective([2*x^3 - 50*x*z^2 + 24*z^3,
            ....:                                 5*y^3 - 53*y*z^2 + 24*z^3, 24*z^3])
            sage: sorted(f.all_periodic_points(prime_bound=[1,20])) # long time
            [(-3 : -1 : 1), (-3 : 0 : 1), (-3 : 1 : 1), (-3 : 3 : 1), (-1 : -1 : 1),
             (-1 : 0 : 1), (-1 : 1 : 1), (-1 : 3 : 1), (0 : 1 : 0), (1 : -1 : 1),
             (1 : 0 : 0), (1 : 0 : 1), (1 : 1 : 1), (1 : 3 : 1), (3 : -1 : 1),
             (3 : 0 : 1), (3 : 1 : 1), (3 : 3 : 1), (5 : -1 : 1), (5 : 0 : 1),
             (5 : 1 : 1), (5 : 3 : 1)]

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: f = DynamicalSystem_projective([-5*x^2 + 4*y^2, 4*x*y])
            sage: sorted(f.all_periodic_points())  # long time
            [(-2 : 1), (-2/3 : 1), (2/3 : 1), (1 : 0), (2 : 1)]

        ::

            sage: # needs sage.rings.number_field
            sage: R.<x> = QQ[]
            sage: K.<w> = NumberField(x^2 - x + 1)
            sage: P.<u,v> = ProjectiveSpace(K,1)
            sage: f = DynamicalSystem_projective([u^2 + v^2, v^2])
            sage: sorted(f.all_periodic_points())                                       # needs sage.rings.function_field
            [(-w + 1 : 1), (w : 1), (1 : 0)]

        ::

            sage: # needs sage.rings.number_field
            sage: R.<x> = QQ[]
            sage: K.<w> = NumberField(x^2 - x + 1)
            sage: P.<u,v> = ProjectiveSpace(K,1)
            sage: f = DynamicalSystem_projective([u^2 + v^2, u*v])
            sage: f.all_periodic_points()                                               # needs sage.rings.function_field
            Traceback (most recent call last):
            ...
            NotImplementedError: rational periodic points for number fields
            only implemented for polynomials

        ::

            sage: # needs sage.rings.number_field
            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: K.<v> = QuadraticField(5)
            sage: phi = QQ.embeddings(K)[0]
            sage: f = DynamicalSystem_projective([x^2 - y^2, y^2])
            sage: sorted(f.all_periodic_points(R=phi))
            [(-1 : 1), (-1/2*v + 1/2 : 1), (0 : 1), (1 : 0), (1/2*v + 1/2 : 1)]

        ::

            sage: P.<x,y,z,w> = ProjectiveSpace(QQ, 3)
            sage: f = DynamicalSystem_projective([x^2 - (3/4)*w^2, y^2 - 3/4*w^2,
            ....:                                 z^2 - 3/4*w^2, w^2])
            sage: sorted(f.all_periodic_points(algorithm='dynatomic'))                  # needs sage.rings.function_field
            [(-1/2 : -1/2 : -1/2 : 1),
             (-1/2 : -1/2 : 3/2 : 1),
             (-1/2 : 3/2 : -1/2 : 1),
             (-1/2 : 3/2 : 3/2 : 1),
             (0 : 0 : 1 : 0),
             (0 : 1 : 0 : 0),
             (0 : 1 : 1 : 0),
             (1 : 0 : 0 : 0),
             (1 : 0 : 1 : 0),
             (1 : 1 : 0 : 0),
             (1 : 1 : 1 : 0),
             (3/2 : -1/2 : -1/2 : 1),
             (3/2 : -1/2 : 3/2 : 1),
             (3/2 : 3/2 : -1/2 : 1),
             (3/2 : 3/2 : 3/2 : 1)]

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSystem_projective([x^2 - 3/4*y^2, y^2])
            sage: sorted(f.all_periodic_points(period_degree_bounds=[2,2]))             # needs sage.rings.function_field
            [(-1/2 : 1), (1 : 0), (3/2 : 1)]

        TESTS::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSystem([x^2 + y^2, x*y])
            sage: f.all_periodic_points(algorithm='banana')
            Traceback (most recent call last):
            ...
            ValueError: algorithm must be 'dynatomic' or 'lifting'
        """
        ring = kwds.pop("R", None)
        if ring is not None:
            # changes to the new ring
            DS = self.change_ring(ring)
            # ensures that the correct method is run, in case user switches to a finite field
            return DS.all_periodic_points(**kwds)
        else:
            DS = self
        PS = DS.domain()
        K = PS.base_ring()
        if K in NumberFields():
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
                g = DS.dehomogenize(1)
                inf = PS([1,0])
                k = 1
                if isinstance(g[0], FractionFieldElement):
                    g = DS.dehomogenize(0)
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
                Fper = F.all_periodic_points(**kwds)
                for P in Fper:
                    #take the 'good' points in the weil restriction and find the
                    #associated number field points.
                    if P[d] == 1:
                        pt = [sum([P[i]*w**i for i in range(d)])]
                        pt.insert(k,1)
                        Q = PS(pt)
                        # for each periodic point get the entire cycle
                        if Q not in periodic_points:
                            #check periodic not preperiodic and add all points in cycle
                            orb = set([Q])
                            Q2 = DS(Q)
                            while Q2 not in orb:
                                orb.add(Q2)
                                Q2 = DS(Q2)
                            if Q2 == Q:
                                periodic_points = periodic_points.union(orb)
                return list(periodic_points)
            else:
                primebound = kwds.pop("prime_bound", [1, 20])
                p = kwds.pop("lifting_prime", 23)
                pd_bounds = kwds.pop("period_degree_bounds", [4,4])
                alg = kwds.pop("algorithm", None)
                periods = kwds.pop("periods", None)
                badprimes = kwds.pop("bad_primes", None)
                num_cpus = kwds.pop("ncpus", ncpus())
                if alg is not None and alg not in ['dynatomic', 'lifting']:
                    raise ValueError("algorithm must be 'dynatomic' or 'lifting'")

                if not isinstance(primebound, (list, tuple)):
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
                    badprimes = DS.primes_of_bad_reduction()
                if periods is None:
                    periods = DS.possible_periods(prime_bound=primebound, bad_primes=badprimes, ncpus=num_cpus)
                PS = DS.domain()
                periodic = set()
                N = PS.ambient_space().dimension_relative()

                if alg != 'lifting':
                    for i in periods[:]:
                        if (alg == 'dynatomic') or ((N == 1)
                                and i <= pd_bounds[0] and DS.degree() <= pd_bounds[1]):
                            periodic.update(DS.periodic_points(i))
                            periods.remove(i)
                    if not periods:
                        return list(periodic)
                while p in badprimes:
                    p = next_prime(p + 1)

                from sage.symbolic.constants import e

                B = e ** DS.height_difference_bound()
                f = DS.change_ring(GF(p))
                all_points = f.possible_periods(True) # return the list of points and their periods.
                pos_points = []
                # check period, remove duplicates
                for i in range(len(all_points)):
                    if all_points[i][1] in periods and all_points[i] not in pos_points:
                        pos_points.append(all_points[i])
                periodic_points = DS.lift_to_rational_periodic(pos_points,B)
                for p,n in periodic_points:
                    for k in range(n):
                        p.normalize_coordinates()
                        periodic.add(p)
                        p = DS(p)
                return list(periodic)
        else:
            raise TypeError("base field must be an absolute number field")

    def all_rational_preimages(self, points):
        r"""
        Given a set of rational points in the domain of this
        dynamical system, return all the rational preimages of those points.

        In others words, all the rational points which have some
        iterate in the set points. This function repeatedly calls
        :meth:`rational_preimages`. If the degree is at least two,
        by Northocott, this is always a finite set. The map must be
        defined over number fields and be an endomorphism.

        INPUT:

        - ``points`` -- list of rational points in the domain of this map

        OUTPUT: list of rational points in the domain of this map

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: f = DynamicalSystem_projective([16*x^2 - 29*y^2, 16*y^2])
            sage: sorted(f.all_rational_preimages([P(-1,4)]))                           # needs sage.rings.function_field
            [(-7/4 : 1), (-5/4 : 1), (-3/4 : 1), (-1/4 : 1), (1/4 : 1), (3/4 : 1),
             (5/4 : 1), (7/4 : 1)]

        ::

            sage: P.<x,y,z> = ProjectiveSpace(QQ,2)
            sage: f = DynamicalSystem_projective([76*x^2 - 180*x*y + 45*y^2 + 14*x*z
            ....:                                   + 45*y*z - 90*z^2,
            ....:                                 67*x^2 - 180*x*y - 157*x*z + 90*y*z,
            ....:                                 -90*z^2])
            sage: sorted(f.all_rational_preimages([P(-9,-4,1)]))                        # needs sage.rings.function_field
            [(-9 : -4 : 1), (0 : -1 : 1), (0 : 0 : 1), (0 : 1 : 1), (0 : 4 : 1),
             (1 : 0 : 1), (1 : 1 : 1), (1 : 2 : 1), (1 : 3 : 1)]

        A non-periodic example ::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: f = DynamicalSystem_projective([x^2 + y^2, 2*x*y])
            sage: sorted(f.all_rational_preimages([P(17,15)]))                          # needs sage.rings.function_field
            [(1/3 : 1), (3/5 : 1), (5/3 : 1), (3 : 1)]

        A number field example::

            sage: # needs sage.rings.number_field
            sage: z = QQ['z'].0
            sage: K.<w> = NumberField(z^3 + (z^2)/4 - (41/16)*z + 23/64)
            sage: P.<x,y> = ProjectiveSpace(K,1)
            sage: f = DynamicalSystem_projective([16*x^2 - 29*y^2, 16*y^2])
            sage: sorted(f.all_rational_preimages([P(16*w^2 - 29, 16)]), key=str)
            [(-w - 1/2 : 1),
             (-w : 1),
             (-w^2 + 21/16 : 1),
             (-w^2 + 29/16 : 1),
             (-w^2 - w + 25/16 : 1),
             (-w^2 - w + 33/16 : 1),
             (w + 1/2 : 1),
             (w : 1),
             (w^2 + w - 25/16 : 1),
             (w^2 + w - 33/16 : 1),
             (w^2 - 21/16 : 1),
             (w^2 - 29/16 : 1)]

        ::

            sage: # needs sage.rings.number_field
            sage: K.<w> = QuadraticField(3)
            sage: P.<u,v> = ProjectiveSpace(K,1)
            sage: f = DynamicalSystem_projective([u^2 + v^2, v^2])
            sage: f.all_rational_preimages(P(4))                                        # needs sage.rings.function_field
            [(-w : 1), (w : 1)]
        """
        if self.domain().base_ring() not in NumberFields():
            raise TypeError("field won't return finite list of elements")
        if not isinstance(points, (list, tuple)):
            points = [points]

        preperiodic = set()
        while points != []:
            P = points.pop()
            preimages = self.rational_preimages(P)
            for i in range(len(preimages)):
                if preimages[i] not in preperiodic:
                    points.append(preimages[i])
                    preperiodic.add(preimages[i])
        return list(preperiodic)

    def all_preperiodic_points(self, **kwds):
        r"""
        Determine the set of rational preperiodic points for
        this dynamical system.

        The map must be defined over `\QQ` and be an endomorphism of
        projective space. If the map is a polynomial endomorphism of
        `\mathbb{P}^1`, i.e. has a totally ramified fixed point, then
        the base ring can be an absolute number field.
        This is done by passing to the Weil restriction.

        The default parameter values are typically good choices for
        `\mathbb{P}^1`. If you are having trouble getting a particular
        map to finish, try first computing the possible periods, then
        try various different values for ``lifting_prime``.

        ALGORITHM:

        - Determines the list of possible periods.

        - Determines the rational periodic points from the possible periods.

        - Determines the rational preperiodic points from the rational
          periodic points by determining rational preimages.

        INPUT: keyword arguments:

        - ``R`` -- (default: domain of dynamical system) the base ring
          over which the periodic points of the dynamical system are found

        - ``prime_bound`` -- (default: ``[1, 20]``) a pair (list or tuple)
          of positive integers that represent the limits of primes to use
          in the reduction step or an integer that represents the upper bound

        - ``lifting_prime`` -- (default: 23) a prime integer; specifies
          modulo which prime to try and perform the lifting

        - ``periods`` -- (optional) a list of positive integers that is
          the list of possible periods

        - ``bad_primes`` -- (optional) a list or tuple of integer primes;
          the primes of bad reduction

        - ``ncpus`` -- (default: all cpus) number of cpus to use in parallel

        - ``period_degree_bounds`` -- (default: ``[4,4]``) a pair of positive integers
          (max period, max degree) for which the dynatomic polynomial should be solved
          for when in dimension 1

        - ``algorithm`` -- (optional) specifies which algorithm to use;
          current options are `dynatomic` and `lifting`; defaults to solving the
          dynatomic for low periods and degrees and lifts for everything else

        OUTPUT: list of rational points in projective space

        EXAMPLES::

            sage: PS.<x,y> = ProjectiveSpace(1,QQ)
            sage: f = DynamicalSystem_projective([x^2 - y^2, 3*x*y])
            sage: sorted(f.all_preperiodic_points())                                    # needs sage.rings.function_field
            [(-2 : 1), (-1 : 1), (-1/2 : 1), (0 : 1), (1/2 : 1), (1 : 0), (1 : 1), (2 : 1)]

        ::

            sage: PS.<x,y> = ProjectiveSpace(1,QQ)
            sage: f = DynamicalSystem_projective([5*x^3 - 53*x*y^2 + 24*y^3, 24*y^3])
            sage: sorted(f.all_preperiodic_points(prime_bound=10))                      # needs sage.rings.function_field
            [(-1 : 1), (0 : 1), (1 : 0), (1 : 1), (3 : 1)]

        ::

            sage: PS.<x,y,z> = ProjectiveSpace(2,QQ)
            sage: f = DynamicalSystem_projective([x^2 - 21/16*z^2, y^2 - 2*z^2, z^2])
            sage: sorted(f.all_preperiodic_points(prime_bound=[1,8],              # long time
            ....:                                 lifting_prime=7, periods=[2]))
            [(-5/4 : -2 : 1), (-5/4 : -1 : 1), (-5/4 : 0 : 1), (-5/4 : 1 : 1), (-5/4 : 2 : 1),
             (-1/4 : -2 : 1), (-1/4 : -1 : 1), (-1/4 : 0 : 1), (-1/4 : 1 : 1), (-1/4 : 2 : 1),
             (1/4 : -2 : 1), (1/4 : -1 : 1), (1/4 : 0 : 1), (1/4 : 1 : 1), (1/4 : 2 : 1),
             (5/4 : -2 : 1), (5/4 : -1 : 1), (5/4 : 0 : 1), (5/4 : 1 : 1), (5/4 : 2 : 1)]

        ::

            sage: # needs sage.rings.number_field
            sage: K.<w> = QuadraticField(33)
            sage: PS.<x,y> = ProjectiveSpace(K,1)
            sage: f = DynamicalSystem_projective([x^2 - 71/48*y^2, y^2])
            sage: sorted(f.all_preperiodic_points())    # long time
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
        ring = kwds.pop("R",None)
        if ring is not None:
            DS = self.change_ring(ring)
        else:
            DS = self
        PS = DS.domain()
        K = PS.base_ring()
        if K not in NumberFields() or not K.is_absolute():
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
            g = DS.dehomogenize(1)
            inf = PS([1,0])
            k = 1
            if isinstance(g[0], FractionFieldElement):
                g = DS.dehomogenize(0)
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
            Fpre = F.all_preperiodic_points(**kwds)
            for P in Fpre:
                #take the 'good' points in the weil restriction and find the
                #associated number field points.
                if P[d] == 1:
                    pt = [sum([P[i]*w**i for i in range(d)])]
                    pt.insert(k,1)
                    Q = PS(pt)
                    # for each preperiodic point get the entire connected component
                    if Q not in preper:
                        for t in DS.connected_rational_component(Q):
                            preper.add(t)
            preper = list(preper)
        else:
            #input error checking done in possible_periods and rational_periodic_points
            badprimes = kwds.pop("bad_primes", None)
            periods = kwds.pop("periods", None)
            primebound = kwds.pop("prime_bound", [1, 20])
            num_cpus = kwds.pop("ncpus", ncpus())
            if badprimes is None:
                badprimes = DS.primes_of_bad_reduction()
            if periods is None:
                #determine the set of possible periods
                periods = DS.possible_periods(prime_bound=primebound,
                                                bad_primes=badprimes, ncpus=num_cpus)
            if periods == []:
                return []  #no rational preperiodic points
            else:
                p = kwds.pop("lifting_prime", 23)
                #find the rational preperiodic points
                T = DS.all_periodic_points(prime_bound=primebound, lifting_prime=p,
                                                  periods=periods, bad_primes=badprimes,
                                                  ncpus=num_cpus, **kwds)
                preper = DS.all_rational_preimages(T) #find the preperiodic points
                preper = list(preper)
        return preper

    def rational_preperiodic_graph(self, **kwds):
        r"""
        Determine the directed graph of the rational preperiodic points
        for this dynamical system.

        The map must be defined over `\QQ` and be an endomorphism of
        projective space. If this map is a polynomial endomorphism of
        `\mathbb{P}^1`, i.e. has a totally ramified fixed point, then
        the base ring can be an absolute number field.
        This is done by passing to the Weil restriction.

        ALGORITHM:

        - Determines the list of possible periods.

        - Determines the rational periodic points from the possible periods.

        - Determines the rational preperiodic points from the rational
          periodic points by determining rational preimages.

        INPUT: keyword arguments:

        - ``prime_bound`` -- (default: ``[1, 20]``) a pair (list or tuple)
          of positive integers that represent the limits of primes to use
          in the reduction step or an integer that represents the upper bound

        - ``lifting_prime`` -- (default: 23) a prime integer; specifies
          modulo which prime to try and perform the lifting

        - ``periods`` -- (optional) a list of positive integers that is
          the list of possible periods

        - ``bad_primes`` -- (optional) a list or tuple of integer primes;
          the primes of bad reduction

        - ``ncpus`` -- (default: all cpus) number of cpus to use in parallel

        OUTPUT:

        A digraph representing the orbits of the rational preperiodic
        points in projective space.

        EXAMPLES::

            sage: PS.<x,y> = ProjectiveSpace(1,QQ)
            sage: f = DynamicalSystem_projective([7*x^2 - 28*y^2, 24*x*y])
            sage: f.rational_preperiodic_graph()                                        # needs sage.rings.function_field
            Looped digraph on 12 vertices

        ::

            sage: PS.<x,y> = ProjectiveSpace(1,QQ)
            sage: f = DynamicalSystem_projective([-3/2*x^3 + 19/6*x*y^2, y^3])
            sage: f.rational_preperiodic_graph(prime_bound=[1,8])                       # needs sage.rings.function_field
            Looped digraph on 12 vertices

        ::

            sage: PS.<x,y,z> = ProjectiveSpace(2,QQ)
            sage: f = DynamicalSystem_projective([2*x^3 - 50*x*z^2 + 24*z^3,
            ....:                                 5*y^3 - 53*y*z^2 + 24*z^3, 24*z^3])
            sage: f.rational_preperiodic_graph(prime_bound=[1,11],  # long time
            ....:                              lifting_prime=13)
            Looped digraph on 30 vertices

        ::

            sage: # needs sage.rings.number_field
            sage: K.<w> = QuadraticField(-3)
            sage: P.<x,y> = ProjectiveSpace(K,1)
            sage: f = DynamicalSystem_projective([x^2 + y^2, y^2])
            sage: f.rational_preperiodic_graph()        # long time
            Looped digraph on 5 vertices
        """
        #input checking done in .rational_preperiodic_points()
        preper = self.all_preperiodic_points(**kwds)
        g = self._preperiodic_points_to_cyclegraph(preper)
        return g

    def connected_rational_component(self, P, n=0):
        r"""
        Compute the connected component of a rational preperiodic
        point ``P`` by this dynamical system.

        Will work for non-preperiodic points if ``n`` is positive.
        Otherwise this will not terminate.

        INPUT:

        - ``P`` -- a rational preperiodic point of this map

        - ``n`` -- (default: 0) integer; maximum distance from ``P`` to
          branch out; a value of 0 indicates no bound

        OUTPUT: list of points connected to ``P`` up to the specified distance

        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: R.<x> = PolynomialRing(QQ)
            sage: K.<w> = NumberField(x^3 + 1/4*x^2 - 41/16*x + 23/64)
            sage: PS.<x,y> = ProjectiveSpace(1,K)
            sage: f = DynamicalSystem_projective([x^2 - 29/16*y^2, y^2])
            sage: P = PS([w,1])
            sage: sorted(f.connected_rational_component(P), key=str)
            [(-w - 1/2 : 1),
             (-w : 1),
             (-w^2 + 21/16 : 1),
             (-w^2 + 29/16 : 1),
             (-w^2 - w + 25/16 : 1),
             (-w^2 - w + 33/16 : 1),
             (w + 1/2 : 1),
             (w : 1),
             (w^2 + w - 25/16 : 1),
             (w^2 + w - 33/16 : 1),
             (w^2 - 21/16 : 1),
             (w^2 - 29/16 : 1)]

        ::

            sage: PS.<x,y,z> = ProjectiveSpace(2,QQ)
            sage: f = DynamicalSystem_projective([x^2 - 21/16*z^2, y^2 - 2*z^2, z^2])
            sage: P = PS([17/16, 7/4, 1])
            sage: f.connected_rational_component(P, 3)                                  # needs sage.rings.function_field
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

    def conjugating_set(self, other, R=None, num_cpus=2):
        r"""
        Return the set of elements in PGL over the base ring
        that conjugates one dynamical system to the other.

        Given two nonconstant rational functions of equal degree,
        determine if there is a rational element of PGL that
        conjugates one rational function to another.

        The optional argument `R` specifies the field of definition
        of the PGL elements. The set is determined
        by taking the fixed points of one map and mapping
        them to permutations of the fixed points of the other map.
        As conjugacy preserves the multipliers as a set, fixed points
        are only mapped to fixed points with the same multiplier.
        If there are not enough fixed points the
        function compares the mapping between rational preimages of
        fixed points and the rational preimages of the preimages of
        fixed points until there are enough points; such that there
        are `n+2` points with all `n+1` subsets linearly independent.

        .. WARNING::

           For degree 1 maps that are conjugate, there is a positive dimensional
           set of conjugations. This function returns only one such element.

        ALGORITHM:

        Implementing invariant set algorithm from the paper [FMV2014]_.
        Uses the set of  `n` th preimages of fixed points, as this set is
        invariant under conjugation to find all elements of PGL that
        take one set to another. Additionally, keeps track of multiplier
        information to reduce the necessary combinatorics.

        INPUT:

        - ``other`` -- a rational function of same degree
          as this map

        - ``R`` -- a field or embedding

        - ``num_cpus`` -- (default: 2) the number of threads to run in parallel;
          increasing ``num_cpus`` can potentially greatly speed up this function

        OUTPUT: set of conjugating `n+1` by `n+1` matrices

        AUTHORS:

        - Original algorithm written by Xander Faber, Michelle Manes,
          Bianca Viray [FMV2014]_.

        - Implemented by Rebecca Lauren Miller as part of GSOC 2016.

        - Algorithmic improvement by Alexander Galarraga as part of GSOC 2021.

        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSystem_projective([x^2 - 2*y^2, y^2])
            sage: m = matrix(QQbar, 2, 2, [-1, 3, 2, 1])
            sage: g = f.conjugate(m)
            sage: f.conjugating_set(g)
            [
            [-1  3]
            [ 2  1]
            ]

        Increasing ``num_cpus`` can speed up computation::

            sage: P.<x,y,z,w> = ProjectiveSpace(QQ, 3)
            sage: f = DynamicalSystem_projective([x^2, y^2, z^2, w^2])
            sage: len(f.conjugating_set(f, num_cpus=3))                                 # needs sage.rings.function_field
            24

        ::

            sage: # needs sage.rings.number_field
            sage: K.<w> = QuadraticField(-1)
            sage: P.<x,y> = ProjectiveSpace(K, 1)
            sage: f = DynamicalSystem_projective([x^2 + y^2, x*y])
            sage: m = matrix(K, 2, 2, [1, 1, 2, 1])
            sage: g = f.conjugate(m)
            sage: sorted(f.conjugating_set(g))
            [
            [-1 -1]  [1 1]
            [ 2  1], [2 1]
            ]

        ::

            sage: # needs sage.rings.number_field
            sage: K.<i> = QuadraticField(-1)
            sage: P.<x,y> = ProjectiveSpace(K, 1)
            sage: D8 = DynamicalSystem_projective([y^3, x^3])
            sage: sorted(D8.conjugating_set(D8))
            [
            [-1  0]  [-i  0]  [ 0 -1]  [ 0 -i]  [0 i]  [0 1]  [i 0]  [1 0]
            [ 0  1], [ 0  1], [ 1  0], [ 1  0], [1 0], [1 0], [0 1], [0 1]
            ]

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: D8 = DynamicalSystem_projective([y^2, x^2])
            sage: D8.conjugating_set(D8)                                                # needs sage.rings.function_field
            Traceback (most recent call last):
            ...
            ValueError: no more rational preimages;
            try extending the base field and trying again

        ::

            sage: P.<x,y> = ProjectiveSpace(GF(7), 1)
            sage: D6 = DynamicalSystem_projective([y^2, x^2])
            sage: sorted(D6.conjugating_set(D6))                                        # needs sage.rings.function_field
            [
            [0 1]  [0 2]  [0 4]  [1 0]  [2 0]  [4 0]
            [1 0], [1 0], [1 0], [0 1], [0 1], [0 1]
            ]

        ::

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: f = DynamicalSystem_projective([x^2 + x*z, y^2, z^2])
            sage: f.conjugating_set(f)                                                  # needs sage.rings.function_field
            [
            [1 0 0]
            [0 1 0]
            [0 0 1]
            ]

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: R = P.coordinate_ring()
            sage: f = DynamicalSystem_projective([R(3), R(4)])
            sage: g = DynamicalSystem_projective([R(5), R(2)])
            sage: m = f.conjugating_set(g)[0]
            sage: f.conjugate(m) == g
            True

        ::

            sage: # needs sage.rings.number_field
            sage: P.<x,y> = ProjectiveSpace(QQbar, 1)
            sage: f = DynamicalSystem_projective([7*x + 12*y, 8*x])
            sage: g = DynamicalSystem_projective([1645*x - 318*y, 8473*x - 1638*y])
            sage: m = f.conjugating_set(g)[0]
            sage: f.conjugate(m) == g
            True

        note that only one possible conjugation is returned::

            sage: P.<x,y,z> = ProjectiveSpace(GF(11), 2)
            sage: f = DynamicalSystem_projective([2*x + 12*y, 11*y + 2*z, x + z])
            sage: m1 = matrix(GF(11), 3, 3, [1,4,1,0,2,1,1,1,1])
            sage: g = f.conjugate(m1)
            sage: f.conjugating_set(g)
            [
            [ 1  0  0]
            [ 9  1  4]
            [ 4 10  8]
            ]

        ::

            sage: # needs sage.rings.number_field
            sage: L.<v> = CyclotomicField(8)
            sage: P.<x,y,z> = ProjectiveSpace(L, 2)
            sage: f = DynamicalSystem_projective([2*x + 12*y, 11*y + 2*z, x + z])
            sage: m1 = matrix(L, 3, 3, [1,4,v^2,0,2,1,1,1,1])
            sage: g = f.conjugate(m1)
            sage: m = f.conjugating_set(g)[0]   # long time
            sage: f.conjugate(m) == g           # long time
            True

        TESTS:

        Make sure the caching problem is fixed, see #28070 ::

            sage: # needs sage.rings.number_field
            sage: K.<i> = QuadraticField(-1)
            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSystem_projective([x^2 - 2*y^2, y^2])
            sage: m = matrix(QQ, 2, 2, [-1, 3, 2, 1])
            sage: g = f.conjugate(m)
            sage: f.conjugating_set(g)                                                  # needs sage.rings.function_field
            [
            [-1  3]
            [ 2  1]
            ]
            sage: f = f.change_ring(K)
            sage: g = g.change_ring(K)
            sage: f.conjugating_set(g)                                                  # needs sage.rings.function_field
            [
            [-1  3]
            [ 2  1]
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
        gens = f[0].parent().gens()
        M = len(gens)
        base = f.base_ring()
        if f.degree() == 0: # all constant maps are conjugate
            zer = [0 for i in range(M-1)]
            m = []
            for i in range(M):
                m1 = copy(zer)
                m1.insert(i, f[i]/g[i])
                m += m1
            return [matrix(base, M, M, m)]
        if f.degree() == 1: # for degree 1 maps, check if matrix representations are similar
            # make matrix forms of f1 and f2
            m1 = matrix(base,M,M,[F.coefficient(var) for F in f for var in gens])
            m2 = matrix(base,M,M,[F.coefficient(var) for F in g for var in gens])
            # Note: det_ratio will be nonzero for invertible f1, f2
            if m1.det() != m2.det():
                det_ratio = m1.det()/m2.det()
                try:
                    det_root = det_ratio.nth_root(M)
                except ValueError: #no root in field
                    return []
                #matrices must have same determinant to be similar, but were in PGL
                #so we can scale so the determinants are equal
                m1 = (1/det_root)*m1
            bol,m = m2.is_similar(m1, transformation=True)
            if bol:
                if m.base_ring() == base:
                    return [m]
                #else is_similar went to algebraic closure
                if base in NumberFields():
                    from sage.rings.qqbar import number_field_elements_from_algebraics
                    K,mK,phi = number_field_elements_from_algebraics([u for t in list(m) for u in t],
                                minimal=True)
                    if K == base:
                        return [matrix(K, M, M, mK)]
                    else: #may be a subfield
                        embeds = K.embeddings(base)
                        if len(embeds) == 0:
                            #not a subfield
                            return []
                        else:
                            for emb in embeds:
                                m_emb = matrix(base, M,M, [emb(u) for u in mK])
                                #check that it is the right embedding
                                if f.conjugate(m_emb) == g:
                                    return [m_emb]
                else: #finite field case
                    #always comes from prime field so can coerce
                    m = matrix(base, M, M, [base(u.as_finite_field_element()[1]) for t in list(m) for u in t])
                    return [m]
            #not similar
            return []
        # sigma invariants are invariant under conjugacy but are only fast in dim 1
        n = f.domain().dimension_relative()
        if (n == 1) and (R in NumberFields() or R in FiniteFields())\
            and (f.sigma_invariants(1) != g.sigma_invariants(1)):
            return []
        tup = conjugating_set_initializer(f, g)
        if tup == []:
            return []
        source, possible_targets = tup
        return conjugating_set_helper(f, g, num_cpus, source, possible_targets)

    def is_conjugate(self, other, R=None, num_cpus=2):
        r"""
        Return whether two dynamical systems are conjugate over their
        base ring (by default) or over the ring `R` entered as an
        optional parameter.

        ALGORITHM:

        Implementing invariant set algorithm from the paper [FMV2014]_.
        Uses the set of  `n` th preimages of fixed points, as this set is
        invariant under conjugation to find all elements of PGL that
        take one set to another. Additionally, keeps track of multiplier
        information to reduce the necessary combinatorics.

        INPUT:

        - ``other`` -- a nonconstant rational function of the same
          degree as this map

        - ``R`` -- a field or embedding

        - ``num_cpus`` -- (default: 2) the number of threads to run in parallel;
          increasing ``num_cpus`` can potentially greatly speed up this function

        OUTPUT: boolean

        AUTHORS:

        - Original algorithm written by Xander Faber, Michelle Manes,
          Bianca Viray [FMV2014]_.

        - Implemented by Rebecca Lauren Miller as part of GSOC 2016.

        - Algorithmic improvement by Alexander Galarraga as part of GSOC 2021.

        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: K.<w> = CyclotomicField(3)
            sage: P.<x,y> = ProjectiveSpace(K, 1)
            sage: D8 = DynamicalSystem_projective([y^2, x^2])
            sage: D8.is_conjugate(D8)
            True

        We can speed up computation by increasing ``num_cpus``::

            sage: P.<x,y,z,w> = ProjectiveSpace(QQ,3)
            sage: f = DynamicalSystem_projective([x^2, y^2, z^2, w^2])
            sage: f.is_conjugate(f, num_cpus=2)                                         # needs sage.rings.function_field
            True

        ::

            sage: # needs sage.rings.number_field
            sage: set_verbose(None)
            sage: P.<x,y> = ProjectiveSpace(QQbar, 1)
            sage: f = DynamicalSystem_projective([x^2 + x*y, y^2])
            sage: m = matrix(QQbar, 2, 2, [1, 1, 2, 1])
            sage: g = f.conjugate(m)
            sage: f.is_conjugate(g)
            True

        ::

            sage: P.<x,y> = ProjectiveSpace(GF(5), 1)
            sage: f = DynamicalSystem_projective([x^3 + x*y^2, y^3])
            sage: m = matrix(GF(5), 2, 2, [1, 3, 2, 9])
            sage: g = f.conjugate(m)
            sage: f.is_conjugate(g)                                                     # needs sage.rings.number_field
            True

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSystem_projective([x^2 + x*y, y^2])
            sage: g = DynamicalSystem_projective([x^3 + x^2*y, y^3])
            sage: f.is_conjugate(g)
            False

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSystem_projective([x^2 + x*y, y^2])
            sage: g = DynamicalSystem_projective([x^2 - 2*y^2, y^2])
            sage: f.is_conjugate(g)                                                     # needs sage.rings.number_field
            False

        ::

            sage: # needs sage.rings.number_field
            sage: P.<x,y> = ProjectiveSpace(QQbar, 1)
            sage: f = DynamicalSystem_projective([7*x + 12*y, 8*x])
            sage: g = DynamicalSystem_projective([1645*x - 318*y, 8473*x - 1638*y])
            sage: f.is_conjugate(g)
            True

        Conjugation is only checked over the base field by default::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSystem_projective([-3*y^2, 3*x^2])
            sage: g = DynamicalSystem_projective([-x^2 - 2*x*y, 2*x*y + y^2])
            sage: f.is_conjugate(g), f.is_conjugate(g, R=QQbar)                         # needs sage.rings.number_field
            (False, True)

        ::

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: f = DynamicalSystem_projective([7*x + 12*y, 8*y + 2*z, x + z])
            sage: m1 = matrix(QQ, 3, 3, [1,4,1,0,2,1,1,1,1])
            sage: g = f.conjugate(m1)
            sage: f.is_conjugate(g)
            True

        ::

            sage: P.<x,y,z> = ProjectiveSpace(GF(7), 2)
            sage: f = DynamicalSystem_projective([2*x + 12*y, 11*y + 2*z, x + z])
            sage: m1 = matrix(GF(7), 3, 3, [1,4,1,0,2,1,1,1,1])
            sage: g = f.conjugate(m1)
            sage: f.is_conjugate(g)
            True

        ::

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: f = DynamicalSystem_projective([2*x^2 + 12*y*x, 11*y*x + 2*y^2, x^2 + z^2])
            sage: m1 = matrix(QQ, 3, 3, [1,4,1,0,2,1,1,1,1])
            sage: g = f.conjugate(m1)
            sage: f.is_conjugate(g)
            True

        TESTS:

        Make sure the caching problem is fixed, see #28070 ::

            sage: # needs sage.rings.number_field
            sage: K.<i> = QuadraticField(5)
            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: f = DynamicalSystem_projective([x^2 - 2*y^2, y^2])
            sage: m = matrix(QQ, 2, 2, [-1, 3, 2, 1])
            sage: g = f.conjugate(m)
            sage: f.is_conjugate(g)
            True
            sage: f = f.change_ring(K)
            sage: g = g.change_ring(K)
            sage: f.is_conjugate(g)
            True
        """
        f = copy(self)
        g = copy(other)
        if R is None:
            R = f.base_ring()
        else:
            f = self.change_ring(R)
            g = other.change_ring(R)
        if not (R in NumberFields() or isinstance(R, sage.rings.abc.AlgebraicField) or R in FiniteFields()):
            raise NotImplementedError("ring must be a number field or finite field")
        try:
            f.normalize_coordinates()
            g.normalize_coordinates()
        except (ValueError):
            pass
        if f.degree() != g.degree(): # checks that maps are of equal degree
            return False
        if f.degree() == 0: # all constant maps are conjugate
            return True
        if f.degree() == 1: # for degree 1 maps, check if matrix representations are similar
            # make matrix forms of f1 and f2
            gens = f[0].parent().gens()
            M = len(gens)
            m1 = matrix(f.base_ring(),M,M,[F.coefficient(var) for F in f for var in gens])
            m2 = matrix(f.base_ring(),M,M,[F.coefficient(var) for F in g for var in gens])
            # Note: det_ratio will be nonzero for invertible f1, f2
            if m1.det() != m2.det():
                det_ratio = m1.det()/m2.det()
                try:
                    det_root = det_ratio.nth_root(M)
                except ValueError: #no root in field
                    return False
                # matrices must have same determinant to be similar, but were in PGL
                # so we can scale to have the determinants equal
                m1 = (1/det_root)*m1
            return m1.is_similar(m2)
        # sigma invariants are invariant under conjugacy but are only fast in dim 1
        n = f.domain().dimension_relative()
        if (n == 1) and (R in NumberFields() or R in FiniteFields())\
          and (f.sigma_invariants(1) != g.sigma_invariants(1)):
            return False
        tup = conjugating_set_initializer(f, g)
        if tup == []:
            return False
        source, possible_targets = tup
        return is_conjugate_helper(f, g, num_cpus, source, possible_targets)

    def is_polynomial(self):
        r"""
        Check to see if the dynamical system has a totally ramified
        fixed point.

        The function must be defined over an absolute number field or a
        finite field.

        OUTPUT: boolean

        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: R.<x> = QQ[]
            sage: K.<w> = QuadraticField(7)
            sage: P.<x,y> = ProjectiveSpace(K, 1)
            sage: f = DynamicalSystem_projective([x**2 + 2*x*y - 5*y**2, 2*x*y])
            sage: f.is_polynomial()
            False

        ::

            sage: # needs sage.rings.number_field
            sage: R.<x> = QQ[]
            sage: K.<w> = QuadraticField(7)
            sage: P.<x,y> = ProjectiveSpace(K, 1)
            sage: f = DynamicalSystem_projective([x**2 - 7*x*y, 2*y**2])
            sage: m = matrix(K, 2, 2, [w, 1, 0, 1])
            sage: f = f.conjugate(m)
            sage: f.is_polynomial()
            True

        ::

            sage: # needs sage.rings.number_field
            sage: K.<w> = QuadraticField(4/27)
            sage: P.<x,y> = ProjectiveSpace(K,1)
            sage: f = DynamicalSystem_projective([x**3 + w*y^3, x*y**2])
            sage: f.is_polynomial()
            False

        ::

            sage: # needs sage.rings.finite_rings
            sage: K = GF(3**2, prefix='w')
            sage: P.<x,y> = ProjectiveSpace(K,1)
            sage: f = DynamicalSystem_projective([x**2 + K.gen()*y**2, x*y])
            sage: f.is_polynomial()
            False

        ::

            sage: PS.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSystem_projective([6*x^2 + 12*x*y + 7*y^2, 12*x*y + 42*y^2])
            sage: f.is_polynomial()
            False

        TESTS:

        See :issue:`25242`::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: F = DynamicalSystem([x^2 + y^2, x*y])
            sage: F2 = F.conjugate(matrix(QQ,2,2, [1,2,3,5]))
            sage: F2.is_polynomial()                                                    # needs sage.libs.pari
            False
        """
        if self.codomain().dimension_relative() != 1:
            raise NotImplementedError("space must have dimension equal to 1")
        K = self.base_ring()
        if K not in FiniteFields() and (K not in NumberFields() or not K.is_absolute()):
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
            for p,exp in Y:
                if p.degree() == 1:
                    if len((g[0]*p[1] + g[1]*p[0]).factor()) == 1:
                        return True
                    G = R(G/(p**exp)) # we already checked this root
                else:
                    u = p #need to extend to get these roots
            if G.degree() != 0:
                #create the next extension
                if K == QQ:
                    from sage.rings.number_field.number_field import NumberField
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
        Return a normal form in the moduli space of dynamical systems.

        Currently implemented only for polynomials. The totally ramified
        fixed point is moved to infinity and the map is conjugated to the form
        `x^n + a_{n-2} x^{n-2} + \cdots + a_{0}`. Note that for finite fields
        we can only remove the `(n-1)`-st term when the characteristic
        does not divide `n`.

        INPUT:

        - ``return_conjugation`` -- boolean (default: ``False``); if ``True``,
          then return the conjugation element of PGL along with the embedding
          into the new field

        OUTPUT:

        - :class:`SchemeMorphism_polynomial`

        - (optional) an element of PGL as a matrix

        - (optional) the field embedding

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSystem_projective([x^2 + 2*x*y - 5*x^2, 2*x*y])
            sage: f.normal_form()
            Traceback (most recent call last):
            ...
            NotImplementedError: map is not a polynomial

        ::

            sage: # needs sage.rings.number_field
            sage: R.<x> = QQ[]
            sage: K.<w> = NumberField(x^2 - 5)
            sage: P.<x,y> = ProjectiveSpace(K,1)
            sage: f = DynamicalSystem_projective([x^2 + w*x*y, y^2])
            sage: g,m,psi = f.normal_form(return_conjugation=True); m
            [     1 -1/2*w]
            [     0      1]
            sage: f.change_ring(psi).conjugate(m) == g
            True

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: f = DynamicalSystem_projective([13*x^2 + 4*x*y + 3*y^2, 5*y^2])
            sage: f.normal_form()                                                       # needs sage.libs.pari
            Dynamical System of Projective Space of dimension 1 over Rational Field
              Defn: Defined on coordinates by sending (x : y) to
                    (5*x^2 + 9*y^2 : 5*y^2)

        ::

            sage: # needs sage.rings.finite_rings
            sage: K = GF(3^3, prefix='w')
            sage: P.<x,y> = ProjectiveSpace(K,1)
            sage: f = DynamicalSystem_projective([x^3 + 2*x^2*y + 2*x*y^2 + K.gen()*y^3, y^3])
            sage: f.normal_form()
            Dynamical System of Projective Space of dimension 1
             over Finite Field in w3 of size 3^3
              Defn: Defined on coordinates by sending (x : y) to
                    (x^3 + x^2*y + x*y^2 + (-w3)*y^3 : y^3)

        ::

            sage: P.<x,y> = ProjectiveSpace(GF(3),1)
            sage: f = DynamicalSystem_projective([2*x**3 + x**2*y, y**3])
            sage: g,m,psi = f.normal_form(return_conjugation=True); psi                 # needs sage.rings.finite_rings
            Ring morphism:
              From: Finite Field of size 3
              To:   Finite Field in z2 of size 3^2
              Defn: 1 |--> 1

        Fixes :issue:`38012` by not forcing univariate polynomial to be univariate::

            sage: R.<z> = PolynomialRing(QQ)
            sage: f = DynamicalSystem_affine(z^2 + z + 1).homogenize(1)
            sage: f.normal_form()
            Dynamical System of Projective Space of dimension 1 over Rational Field
             Defn: Defined on coordinates by sending (x0 : x1) to
                   (x0^2 + 5/4*x1^2 : x1^2)
        """
        # defines the field of fixed points
        if self.codomain().dimension_relative() != 1:
            raise NotImplementedError("space must have dimension equal to 1")
        K = self.base_ring()
        if K not in FiniteFields() and (K not in NumberFields() or not K.is_absolute()):
            raise NotImplementedError("must be over an absolute number field or finite field")
        if K in FiniteFields():
            q = K.characteristic()
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
            if isinstance(G.parent(), MPolynomialRing_base):
                G = G.polynomial(G.variable(0))
        else:
            #no other fixed points
            raise NotImplementedError("map is not a polynomial")
        #check other fixed points
        while not done:
            Y = G.factor()
            R = G.parent()
            done = True
            for p,exp in Y:
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
                    from sage.rings.number_field.number_field import NumberField
                    L = NumberField(u, 't'+str(i))
                    i += 1
                    phi = K.embeddings(L)[0]
                    K = L
                elif K in FiniteFields():
                    K, phi = K.extension(G.degree(), map=True)
                else:
                    L = K.extension(u, 't'+str(i))
                    i += 1
                    phi1 = K.embeddings(L)[0]
                    K = L
                    L = K.absolute_field('t'+str(i))
                    i += 1
                    phi = K.embeddings(L)[0]*phi1
                    K = L
                psi = phi * psi
                #switch to the new field
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
                M, phi = N.extension(d-1, map=True)
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
            mc2 = matrix(M, 2, 2, [1, M((-gcc[0].coefficient([d-1, 1])
                / (d*gcc[1].coefficient([0, d]))).constant_coefficient()), 0, 1])
        else:
            mc2 = mc.parent().one()
        gccc = gcc.conjugate(mc2)
        if return_conjugation:
            return gccc, m * mc * mc2, psi
        return gccc

    def potential_good_reduction(self, prime, return_conjugation=False):
        r"""
        Return ``True`` if this dynamical system has potential good reduction at ``prime``.

        A dynamical system has good reduction at ``prime`` if after the coefficients
        are reduced modulo ``prime`` the degree remains the same. A dynamical system
        `f` has `\textit{potential}` good reduction if there exists
        `\phi \in PGL(n,\overline{K})` such that `\phi^{-1} \circ f \circ \phi`
        has good reduction.

        If this dynamical system `f` has potential good reduction at ``prime``,
        a dynamical system `g = \phi^{-1} \circ f \circ \phi` which has good
        reduction at ``prime`` is returned.

        This dynamical system must have as its domain `\mathbb{P}^1(K)`, where
        `K` is a number field.

        INPUT:

        - ``prime`` -- a prime ideal of the field of definition of the fixed
          points of the map, or a prime number in `\QQ` if the field of definition
          of the fixed points is `\QQ`.

        - ``return_conjugation`` -- boolean (default: ``False``); if set to ``True``,
          the `PGL_2` map used to achieve good reduction will be returned

        OUTPUT: a tuple:

        - The first element is:
          - ``False`` if this dynamical system does not have potential good reduction.
          - ``True`` if this dynamical system does have potential good reduction.

        - The second element is:
          - ``None`` if this dynamical system does not have potential good reduction.
          - A dynamical system with good reduction at ``prime`` otherwise.

        - If ``return_conjugation`` is ``True``, then the tuple will have a third element, which is:
          - ``None`` if this dynamical system does not have potential good reduction.
          - The `PGL_2` map used to achieve good reduction otherwise.

        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: system = DynamicalSystem_projective([x^2 - y^2, 2*x*y])
            sage: prime = system.field_of_definition_periodic(1).prime_above(2)
            sage: new_system = system.potential_good_reduction(prime)[1]
            sage: new_system
            Dynamical System of Projective Space of dimension 1 over Number Field
             in a with defining polynomial x^2 + 1
              Defn: Defined on coordinates by sending (x : y) to
                    ((-1/2*a)*x^2 + (-5/2*a)*y^2 : (-a)*x*y + y^2)

        Note that this map has good reduction at 2::

            sage: new_system.resultant()                                                # needs sage.rings.number_field
            1

        Using ``return_conjugation``, we can get the conjugation that achieves good reduction::

            sage: conj = system.potential_good_reduction(prime, True)[2]; conj          # needs sage.rings.number_field
            [-1/2*a    1/2]
            [     0      1]

        We can check that this conjugation achieves good reduction::

            sage: system.conjugate(conj).resultant()                                    # needs sage.rings.number_field
            1

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: system = DynamicalSystem_projective([3^4*x^3 + 3*x*y^2 + y^3, 3^6*y^3])
            sage: prime = system.field_of_definition_periodic(1).prime_above(3)         # needs sage.rings.number_field
            sage: system.potential_good_reduction(prime)                                # needs sage.rings.number_field
            (False, None)

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: system = DynamicalSystem_projective([x^5 - x*y^4, 5*y^5])
            sage: prime = system.field_of_definition_periodic(1).prime_above(5)         # needs sage.rings.number_field
            sage: system.potential_good_reduction(prime)                                # needs sage.rings.number_field
            (False, None)

        TESTS::

            sage: # needs sage.rings.number_field
            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: R.<z> = QQ[]
            sage: A.<a> = NumberField(z^2 + 1)
            sage: prime = A.prime_above(2)
            sage: system = DynamicalSystem_projective([x^2 - y^2, 2*x*y])
            sage: system.potential_good_reduction(prime)
            (True,
             Dynamical System of Projective Space of dimension 1 over
              Number Field in a with defining polynomial x^2 + 1
               Defn: Defined on coordinates by sending (x : y) to
                     ((-1/2*a)*x^2 + (-5/2*a)*y^2 : (-a)*x*y + y^2))

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: system = DynamicalSystem_projective([3^5*x^3 + x^2*y - 3^5*x*y^2,
            ....:                                      -3^5*x^2*y + x*y^2 + 3^5*y^3])
            sage: system.potential_good_reduction(3, return_conjugation=True)  # long time
            (False, None, None)

        ::

            sage: # needs sage.rings.number_field
            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: system = DynamicalSystem_projective([x**5 - 11*y**5, x**4*y])
            sage: B, new_sys, conj = system.potential_good_reduction(11, True)
            sage: system.conjugate(conj).resultant() == 1
            True
            sage: system.conjugate(conj) == new_sys
            True

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: system = DynamicalSystem_projective([3*x^2 + x*y + y^2, 9*y^2])
            sage: prime = system.field_of_definition_periodic(1).prime_above(3)         # needs sage.rings.number_field
            sage: system.potential_good_reduction(prime)                                # needs sage.rings.number_field
            (False, None)
        """
        if self.domain().base_ring() not in NumberFields():
            raise ValueError('dynamical system must be defined over number field')

        field_of_definition_periodic = self.field_of_definition_periodic(1)

        if not (isinstance(prime, NumberFieldFractionalIdeal) or prime in QQ):
            raise TypeError('prime must be an ideal of a number field or an element of QQ')
        if prime not in QQ:
            if prime.number_field() != field_of_definition_periodic:
                K = prime.number_field()
                old_parent = K.defining_polynomial().parent()
                new_parent = field_of_definition_periodic.defining_polynomial().parent()
                hom = old_parent.hom([new_parent.gens()[0]])
                L = field_of_definition_periodic
                if hom(K.defining_polynomial()) != L.defining_polynomial():
                    raise ValueError('prime ideal of %s ' % K +
                        'but field of definition of fixed points is %s. ' % L +
                        'see documentation for examples')
                embedding = K.embeddings(field_of_definition_periodic)[0]
                prime = embedding(prime)
        else:
            if field_of_definition_periodic is not QQ:
                raise ValueError('field of definition of fixed ' +
                    'points is %s but prime is in QQ. ' % field_of_definition_periodic)

        system = self.change_ring(field_of_definition_periodic)
        fixed_points = system.periodic_points(1)
        multipliers = [system.multiplier(i, 1)[0][0] for i in fixed_points]
        indifferent_point = None
        for mult in multipliers:
            if field_of_definition_periodic is not QQ:
                valuation = mult.valuation(prime) / prime.absolute_ramification_index()
            else:
                valuation = mult.valuation(prime)
            if valuation < 0:
                if return_conjugation:
                    return (False, None, None)
                else:
                    return (False, None)
            elif valuation == 0:
                indifferent_point = fixed_points[multipliers.index(mult)]
        if indifferent_point is not None:
            point = indifferent_point
            field_of_definition, embedding_preimage = system.field_of_definition_preimage(point, 2, True)
            system = system.change_ring(embedding_preimage)
            point = point.change_ring(embedding_preimage)
            preimages = [point]
            for i in [1,2]:
                preimages_of_point = system.rational_preimages(point, 1)
                for preimage in preimages_of_point:
                    if preimage != point:
                        preimages.append(preimage)
                        point = preimage
                        break
        else:
            preimages = [fixed_points[0], fixed_points[1], fixed_points[2]]
            field_of_definition = field_of_definition_periodic
        P = ProjectiveSpace(field_of_definition,1)
        preimages = [P(i) for i in preimages]
        conjugation = P.point_transformation_matrix(preimages,[P(0),P(1),P([1,0])])
        new_system = system.change_ring(field_of_definition)
        new_system = new_system.conjugate(conjugation)
        res = new_system.resultant()
        if 'embedding_preimage' in locals():
            check_value = res.valuation(field_of_definition.prime_above(embedding_preimage(prime)))
        else:
            check_value = res.valuation(field_of_definition.prime_above(prime))
        if check_value != 0:
            if return_conjugation:
                return (False, None, None)
            else:
                return (False, None)
        if return_conjugation:
            return (True, new_system, conjugation)
        return (True, new_system)

    def reduce_base_field(self):
        """
        Return this map defined over the field of definition of the coefficients.

        The base field of the map could be strictly larger than
        the field where all of the coefficients are defined. This function
        reduces the base field to the minimal possible. This can be done when
        the base ring is a number field, QQbar, a finite field, or algebraic
        closure of a finite field.

        OUTPUT: a dynamical system

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<t> = GF(2^3)
            sage: P.<x,y,z> = ProjectiveSpace(K, 2)
            sage: f = DynamicalSystem_projective([x^2 + y^2, y^2, z^2+z*y])
            sage: f.reduce_base_field()
            Dynamical System of Projective Space of dimension 2 over Finite Field of size 2
              Defn: Defined on coordinates by sending (x : y : z) to
                    (x^2 + y^2 : y^2 : y*z + z^2)

        ::

            sage: # needs sage.rings.number_field sage.symbolic
            sage: P.<x,y,z> = ProjectiveSpace(QQbar, 2)
            sage: f = DynamicalSystem_projective([x^2 + QQbar(sqrt(3))*y^2,
            ....:                                 y^2, QQbar(sqrt(2))*z^2])
            sage: f.reduce_base_field()
            Dynamical System of Projective Space of dimension 2 over Number Field in a
             with defining polynomial y^4 - 4*y^2 + 1 with a = -0.5176380902050415?
              Defn: Defined on coordinates by sending (x : y : z) to
                    (x^2 + (-a^2 + 2)*y^2 : y^2 : (a^3 - 3*a)*z^2)

        ::

            sage: # needs sage.rings.number_field
            sage: R.<x> = QQ[]
            sage: K.<v> = NumberField(x^3 - 2, embedding=(x^3 - 2).roots(ring=CC)[0][0])
            sage: R.<x> = QQ[]
            sage: L.<w> = NumberField(x^6 + 9*x^4 - 4*x^3 + 27*x^2 + 36*x + 31,
            ....:                     embedding=(x^6 + 9*x^4 - 4*x^3
            ....:                                 + 27*x^2 + 36*x + 31).roots(ring=CC)[0][0])
            sage: P.<x,y> = ProjectiveSpace(L,1)
            sage: f = DynamicalSystem([L(v)*x^2 + y^2, x*y])
            sage: f.reduce_base_field().base_ring().is_isomorphic(K)
            True

        ::

            sage: # needs sage.rings.number_field
            sage: K.<v> = CyclotomicField(5)
            sage: A.<x,y> = ProjectiveSpace(K, 1)
            sage: f = DynamicalSystem_projective([3*x^2 + y^2, x*y])
            sage: f.reduce_base_field()
            Dynamical System of Projective Space of dimension 1 over Rational Field
              Defn: Defined on coordinates by sending (x : y) to
                    (3*x^2 + y^2 : x*y)
        """
        return self.as_scheme_morphism().reduce_base_field().as_dynamical_system()

    def is_newton(self, return_conjugation=False):
        r"""
        Return whether ``self`` is a Newton map.

        A map `g` is *Newton* if it is conjugate to a map of the form
        `f(z) = z - \frac{p(z)}{p'(z)}` after dehomogenization,
        where `p(z)` is a squarefree polynomial.

        INPUT:

        - ``return_conjugation`` -- boolean (default: ``False``); if the map is Newton
          and ``True``, then return the conjugation that moves this map to
          the above form

        OUTPUT:

        A boolean. If ``return_conjugation`` is ``True``, then this also
        returns the conjugation as a matrix if ``self`` is Newton or ``None``
        otherwise.

        The conjugation may be defined over an extension if the map has
        fixed points not defined over the base field.

        EXAMPLES::

            sage: A.<z> = AffineSpace(QQ, 1)
            sage: f = DynamicalSystem_affine([z - (z^2 + 1)/(2*z)])
            sage: F = f.homogenize(1)
            sage: F.is_newton(return_conjugation=True)                                  # needs sage.rings.number_field
            (
                  [1 0]
            True, [0 1]
            )

        ::

            sage: A.<z> = AffineSpace(QQ, 1)
            sage: f = DynamicalSystem_affine([z^2 + 1])
            sage: F = f.homogenize(1)
            sage: F.is_newton()                                                         # needs sage.rings.number_field
            False
            sage: F.is_newton(return_conjugation=True)                                  # needs sage.rings.number_field
            (False, None)

        ::

            sage: PP.<x,y> = ProjectiveSpace(QQ, 1)
            sage: F = DynamicalSystem_projective([-4*x^3 - 3*x*y^2, -2*y^3])
            sage: F.is_newton(return_conjugation=True)[1]                               # needs sage.rings.number_field
            [   0    1]
            [-4*a  2*a]

        ::

            sage: # needs sage.rings.number_field
            sage: K.<zeta> = CyclotomicField(2*4)
            sage: A.<z> = AffineSpace(K, 1)
            sage: f = DynamicalSystem_affine(z-(z^3+zeta*z)/(3*z^2+zeta))
            sage: F = f.homogenize(1)
            sage: F.is_newton()
            True
        """
        if self.degree() == 1:
            raise NotImplementedError("degree one Newton maps are trivial")
        if self.base_ring() not in NumberFields():
            raise NotImplementedError("only implemented over number fields")
        # check if Newton map
        sigma_1 = self.sigma_invariants(1)
        d = ZZ(self.degree())
        Newton_sigma = [d/(d-1)] + [0] * d  # almost Newton
        if sigma_1 != Newton_sigma:
            if return_conjugation:
                return False, None
            else:
                return False
        from sage.rings.qqbar import QQbar
        Fbar = self.change_ring(QQbar)
        Pbar = Fbar.domain()
        fixed = Fbar.periodic_points(1)
        for Q in fixed:
            if Fbar.multiplier(Q, 1) != 0:
                inf = Q
                break
        if inf != Pbar([1,0]):
            # need to move to inf to infinity
            fixed.remove(inf)
            source = [inf] + fixed[:2]
            target = [Pbar([1, 0]), Pbar([0, 1]), Pbar([1, 1])]
            M = Pbar.point_transformation_matrix(source, target)
            M = M.inverse()
            Newton = Fbar.conjugate(M)
            K, el, psi = number_field_elements_from_algebraics([t for r in M for t in r])
            M = matrix(M.nrows(), M.ncols(), el)
            Newton = Newton._number_field_from_algebraics()
        else:
            Newton = self
            M = matrix(QQ, 2, 2, [1,0,0,1])
        N_aff = Newton.dehomogenize(1)
        z = N_aff.domain().gen(0)
        Npoly = (z - N_aff[0]).numerator()
        if return_conjugation:
            if Npoly.derivative(z) == (z - N_aff[0]).denominator():
                return True, M
            else:
                return False, None
        else:
            return Npoly.derivative(z) == (z - N_aff[0]).denominator()


class DynamicalSystem_projective_finite_field(DynamicalSystem_projective_field,
                                              SchemeMorphism_polynomial_projective_space_finite_field):

    def is_postcritically_finite(self, **kwds):
        r"""
        Every point is postcritically finite in a finite field.

        INPUT: None. ``kwds`` is to parallel the overridden function

        OUTPUT: the boolean ``True``

        EXAMPLES::

            sage: P.<x,y,z> = ProjectiveSpace(GF(5),2)
            sage: f = DynamicalSystem_projective([x^2 + y^2, y^2, z^2 + y*z], domain=P)
            sage: f.is_postcritically_finite()
            True

        ::

            sage: P.<x,y> = ProjectiveSpace(GF(13),1)
            sage: f = DynamicalSystem_projective([x^4 - x^2*y^2 + y^4, y^4])
            sage: f.is_postcritically_finite(use_algebraic_closure=False)
            True
        """
        return True

    def _is_preperiodic(self, P, **kwds):
        r"""
        Every point in a finite field is preperiodic.

        INPUT:

        - ``P`` -- a point in the domain of this map

        keywords:

        - ``return_period`` -- boolean (default: ``False``); controls if
          the period is returned

        OUTPUT: the boolean ``True`` or a tuple ``(m,n)`` of integers

        EXAMPLES::

            sage: P.<x,y,z> = ProjectiveSpace(GF(5),2)
            sage: f = DynamicalSystem_projective([x^2 + y^2, y^2, z^2 + y*z], domain=P)
            sage: f._is_preperiodic(P(2,1,2))
            True

        ::

            sage: P.<x,y,z> = ProjectiveSpace(GF(5),2)
            sage: f = DynamicalSystem_projective([x^2 + y^2, y^2, z^2 + y*z], domain=P)
            sage: f._is_preperiodic(P(2,1,2), return_period=True)
            (0, 6)
        """
        return_period = kwds.pop('return_period', False)
        if return_period:
            return self.orbit_structure(P)
        else:
            return True

    def orbit_structure(self, P):
        r"""
        Return the pair ``(m,n)``, where ``m`` is the preperiod and ``n``
        is the period of the point ``P`` by this dynamical system.

        Every point is preperiodic over a finite field so every point
        will be preperiodic.

        INPUT:

        - ``P`` -- a point in the domain of this map

        OUTPUT: a tuple ``(m,n)`` of integers

        EXAMPLES::

            sage: P.<x,y,z> = ProjectiveSpace(GF(5),2)
            sage: f = DynamicalSystem_projective([x^2 + y^2, y^2, z^2 + y*z], domain=P)
            sage: f.orbit_structure(P(2,1,2))
            (0, 6)

        ::

            sage: P.<x,y,z> = ProjectiveSpace(GF(7),2)
            sage: X = P.subscheme(x^2 - y^2)
            sage: f = DynamicalSystem_projective([x^2, y^2, z^2], domain=X)
            sage: f.orbit_structure(X(1,1,2))                                           # needs sage.rings.function_field
            (0, 2)

        ::

            sage: P.<x,y> = ProjectiveSpace(GF(13),1)
            sage: f = DynamicalSystem_projective([x^2 - y^2, y^2], domain=P)
            sage: f.orbit_structure(P(3,4))
            (2, 3)

        ::

            sage: # needs sage.rings.finite_rings
            sage: R.<t> = GF(13^3)
            sage: P.<x,y> = ProjectiveSpace(R,1)
            sage: f = DynamicalSystem_projective([x^2 - y^2, y^2], domain=P)
            sage: f.orbit_structure(P(t, 4))
            (11, 6)
        """
        orbit = []
        index = 1
        Q = copy(P)
        Q.normalize_coordinates()
        F = copy(self)
        F.normalize_coordinates()
        while Q not in orbit:
            orbit.append(Q)
            Q = F(Q)
            Q.normalize_coordinates()
            index += 1
        I = orbit.index(Q)
        return (I, index-I-1)

    def cyclegraph(self):
        r"""
        Return the digraph of all orbits of this dynamical system.

        Over a finite field this is a finite graph. For subscheme domains, only points
        on the subscheme whose image are also on the subscheme are in the digraph.

        OUTPUT: a digraph

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(GF(13),1)
            sage: f = DynamicalSystem_projective([x^2 - y^2, y^2])
            sage: f.cyclegraph()                                                        # needs sage.graphs
            Looped digraph on 14 vertices

        ::

            sage: # needs sage.rings.finite_rings
            sage: P.<x,y,z> = ProjectiveSpace(GF(3^2,'t'),2)
            sage: f = DynamicalSystem_projective([x^2 + y^2, y^2, z^2 + y*z])
            sage: f.cyclegraph()                                                        # needs sage.graphs
            Looped digraph on 91 vertices

        ::

            sage: P.<x,y,z> = ProjectiveSpace(GF(7),2)
            sage: X = P.subscheme(x^2 - y^2)
            sage: f = DynamicalSystem_projective([x^2, y^2, z^2], domain=X)
            sage: f.cyclegraph()                                                        # needs sage.graphs
            Looped digraph on 15 vertices

        ::

            sage: P.<x,y,z> = ProjectiveSpace(GF(3),2)
            sage: f = DynamicalSystem_projective([x*z - y^2, x^2 - y^2, y^2 - z^2])
            sage: f.cyclegraph()                                                        # needs sage.graphs
            Looped digraph on 13 vertices

        ::

            sage: P.<x,y,z> = ProjectiveSpace(GF(3),2)
            sage: X = P.subscheme([x - y])
            sage: f = DynamicalSystem_projective([x^2 - y^2, x^2 - y^2, y^2 - z^2], domain=X)
            sage: f.cyclegraph()                                                        # needs sage.graphs
            Looped digraph on 4 vertices
        """
        V = []
        E = []
        if isinstance(self.domain(), ProjectiveSpace_ring):
            for P in self.domain():
                V.append(P)
                try:
                    Q = self(P)
                    Q.normalize_coordinates()
                    E.append([Q])
                except ValueError: #indeterminacy
                    E.append([])
        else:
            X = self.domain()
            for P in X.ambient_space():
                try:
                    XP = X.point(P)
                    V.append(XP)
                    try:
                        Q = self(XP)
                        Q.normalize_coordinates()
                        E.append([Q])
                    except ValueError: #indeterminacy
                        E.append([])
                except TypeError:  # not a point on the scheme
                    pass
        from sage.graphs.digraph import DiGraph
        g = DiGraph(dict(zip(V, E)), loops=True)
        return g

    def possible_periods(self, return_points=False):
        r"""
        Return the list of possible minimal periods of a periodic point
        over `\QQ` and (optionally) a point in each cycle.

        ALGORITHM:

        See [Hutz2009]_.

        INPUT:

        - ``return_points`` -- boolean (default: ``False``); if ``True``,
          then return the points as well as the possible periods

        OUTPUT:

        A list of positive integers, or a list of pairs of projective
        points and periods if ``return_points`` is ``True``.

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(GF(23),1)
            sage: f = DynamicalSystem_projective([x^2 - 2*y^2, y^2])
            sage: f.possible_periods()                                                  # needs sage.libs.pari
            [1, 5, 11, 22, 110]

        ::

            sage: P.<x,y> = ProjectiveSpace(GF(13),1)
            sage: f = DynamicalSystem_projective([x^2 - y^2, y^2])
            sage: sorted(f.possible_periods(True))                                      # needs sage.libs.pari
            [[(0 : 1), 2], [(1 : 0), 1], [(3 : 1), 3], [(3 : 1), 36]]

        ::

            sage: PS.<x,y,z> = ProjectiveSpace(2,GF(7))
            sage: f = DynamicalSystem_projective([-360*x^3 + 760*x*z^2,
            ....:                                 y^3 - 604*y*z^2 + 240*z^3, 240*z^3])
            sage: f.possible_periods()                                                  # needs sage.libs.pari
            [1, 2, 4, 6, 12, 14, 28, 42, 84]

        .. TODO::

            - do not return duplicate points

            - improve hash to reduce memory of point-table
        """
        return _fast_possible_periods(self, return_points)

    def automorphism_group(self, **kwds):
        r"""
        Return the subgroup of `PGL2` that is the automorphism group of this
        dynamical system.

        The automorphism group is the set of `PGL2` elements that fixed the map under conjugation.

        For dimension 1, see [FMV2014]_ for the algorithm.

        For dimension greater than 1, we compute the conjugating set of this
        dynamical system with itself.

        INPUT:

        The following keywords are used when the dimension of the domain
        is greater than 1:

        - ``num_cpus`` -- (default: 2) the number of threads to use; setting to a
          larger number can greatly speed up this function

        The following keywords are used when the dimension of the domain is 1:

        - ``absolute`` -- boolean (default: ``False``); if ``True``, then
          return the absolute automorphism group and a field of definition

        - ``iso_type`` -- boolean (default: ``False``); if ``True``, then
          return the isomorphism type of the automorphism group

        - ``return_functions`` -- boolean (default: ``False``); ``True``
          returns elements as linear fractional transformations and
          ``False`` returns elements as `PGL2` matrices

        OUTPUT: list of elements of the automorphism group

        AUTHORS:

        - Original algorithm written by Xander Faber, Michelle Manes,
          Bianca Viray

        - Modified by Joao Alberto de Faria, Ben Hutz, Bianca Thompson

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: R.<x,y> = ProjectiveSpace(GF(7^3,'t'),1)
            sage: f = DynamicalSystem_projective([x^2 - y^2, x*y])
            sage: f.automorphism_group()
            [
            [1 0]  [6 0]
            [0 1], [0 1]
            ]

        ::

            sage: # needs sage.rings.finite_rings
            sage: R.<x,y> = ProjectiveSpace(GF(3^2,'t'),1)
            sage: f = DynamicalSystem_projective([x^3, y^3])
            sage: lst, label = f.automorphism_group(return_functions=True,      # long time
            ....:                                   iso_type=True)
            sage: sorted(lst, key=str), label                                   # long time
            ([(2*x + 1)/(x + 1),
              (2*x + 1)/x,
              (2*x + 2)/(x + 2),
              (2*x + 2)/x,
              (x + 1)/(x + 2),
              (x + 1)/x,
              (x + 2)/(x + 1),
              (x + 2)/x,
              1/(x + 1),
              1/(x + 2),
              1/x,
              2*x,
              2*x + 1,
              2*x + 2,
              2*x/(x + 1),
              2*x/(x + 2),
              2/(x + 1),
              2/(x + 2),
              2/x,
              x,
              x + 1,
              x + 2,
              x/(x + 1),
              x/(x + 2)],
             'PGL(2,3)')

        ::

            sage: # needs sage.rings.finite_rings
            sage: R.<x,y> = ProjectiveSpace(GF(2^5,'t'),1)
            sage: f = DynamicalSystem_projective([x^5, y^5])
            sage: f.automorphism_group(return_functions=True, iso_type=True)
            ([x, 1/x], 'Cyclic of order 2')

        ::

            sage: # needs sage.rings.finite_rings
            sage: R.<x,y> = ProjectiveSpace(GF(3^4,'t'),1)
            sage: f = DynamicalSystem_projective([x^2 + 25*x*y + y^2, x*y + 3*y^2])
            sage: f.automorphism_group(absolute=True)
            [Univariate Polynomial Ring in w over Finite Field in b of size 3^4,
             [
             [1 0]
             [0 1]
             ]]

        ::

            sage: R.<x,y,z> = ProjectiveSpace(GF(5), 2)
            sage: f = DynamicalSystem_projective([x^3 + x*z^2, y^3 + y*z^2, z^3])
            sage: all([f.conjugate(m) == f for m in f.automorphism_group()])            # needs sage.rings.function_field
            True
        """
        absolute = kwds.get('absolute', False)
        iso_type = kwds.get('iso_type', False)
        return_functions = kwds.get('return_functions', False)
        num_cpus = kwds.get('num_cpus', 2)
        if self.domain().dimension_relative() != 1:
            return self.conjugating_set(self, num_cpus)
        else:
            f = self.dehomogenize(1)
            z = f[0].parent().gen()
        self.normalize_coordinates()
        if (self.degree() == 1) or (self.degree() == 0):
            raise NotImplementedError("rational function of degree 1 not implemented")
        if f[0].denominator() != 1:
            F = f[0].numerator().polynomial(z) / f[0].denominator().polynomial(z)
        else:
            F = f[0].numerator().polynomial(z)
        return automorphism_group_FF(F, absolute, iso_type, return_functions)

    def all_periodic_points(self, **kwds):
        r"""
        Return a list of all periodic points over a finite field.

        INPUT: keyword arguments:

        - ``R`` -- (default: base ring of dynamical system) the base ring
          over which the periodic points of the dynamical system are found

        OUTPUT: list of elements which are periodic

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: P.<x,y> = ProjectiveSpace(GF(5^2),1)
            sage: f = DynamicalSystem_projective([x^2 + y^2, x*y])
            sage: f.all_periodic_points()
            [(1 : 0), (z2 + 2 : 1), (4*z2 + 3 : 1)]

        ::

            sage: P.<x,y,z> = ProjectiveSpace(GF(5),2)
            sage: f = DynamicalSystem_projective([x^2 + y^2 + z^2, x*y + x*z, z^2])
            sage: f.all_periodic_points()
            [(1 : 0 : 0),
             (0 : 0 : 1),
             (1 : 0 : 1),
             (2 : 1 : 1),
             (1 : 4 : 1),
             (3 : 0 : 1),
             (0 : 3 : 1)]

        ::

            sage: P.<x,y> = ProjectiveSpace(GF(3), 1)
            sage: f = DynamicalSystem_projective([x^2 - y^2, y^2])
            sage: f.all_periodic_points(R=GF(3^2, 't'))                                 # needs sage.rings.finite_rings
            [(1 : 0), (0 : 1), (2 : 1), (t : 1), (2*t + 1 : 1)]
        """
        R = kwds.pop("R", None)
        if R is None:
            DS = self
        else:
            DS = self.change_ring(R)
            return DS.all_periodic_points(**kwds)  #ensures that the correct method is run, in case user switches to infinite fields
        return _all_periodic_points(DS)
